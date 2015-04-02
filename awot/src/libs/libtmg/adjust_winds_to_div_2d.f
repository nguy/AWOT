      SUBROUTINE ADJUST_WINDS_TO_DIV_2D(U,SD_U,V,SD_V,DIV,MAXX_IN,
     $MAXY_IN,DELX,DELY,NX,NY,BADDATA,MAXX_OUT,MAXY_OUT,U_ADJ,V_ADJ)

C  Thomas Matejka NOAA/NSSL 9 February 1995

C  This subroutine variationally adjusts two velocity component fields
C  so that they agree with a specified two-dimensional divergence field.
C  The solution is formulated so that the adjusted velocity component
C  fields are as close as possible to the input velocity component
C  fields.

C  Input:

C  U is a two-dimensional real array.  U(I,J) specifies the value, at
C  the Ith grid point in the first dimension and the Jth grid point in
C  the second dimension, of the velocity component field in the first
C  dimension to be adjusted.  If it is missing, it should equal BADDATA.

C  SD_U is a two-dimensional real array.  SD_U(I,J) specifies the
C  standard deviation of U(I,J).  If it is missing, it should equal
C  BADDATA.  SD_U should not be missing for U(I,J) to be usable.

C  V is a two-dimensional real array.  V(I,J) specifies the value, at
C  the Ith grid point in the first dimension and the Jth grid point in
C  the second dimension, of the velocity component field in the second
C  dimension to be adjusted.  If it is missing, it should equal BADDATA.

C  SD_V is a two-dimensional real array.  SD_V(I,J) specifies the
C  standard deviation of V(I,J).  If it is missing, it should equal
C  BADDATA.  SD_V should not be missing for V(I,J) to be usable.

C  DIV is a two-dimensional real array.  DIV(I,J) specifies the
C  two-dimensional divergence of the wind at the Ith grid point in the
C  first dimension and the Jth grid point in the second dimension.  If
C  it is missing, it should equal BADDATA.

C  MAXX_IN is an integer variable that specifies the first dimension of
C  U, SD_U, V, SD_V, and DIV in the calling program.

C  MAXY_IN is an integer variable that specifies the second dimension of
C  U, SD_U, V, SD_V, and DIV in the calling program.

C  NX is an integer variable that specifies the number of grid points in
C  the first dimension for U, SD_U, V, SD_V, DIV, U_ADJ, and V_ADJ.

C  NY is an integer variable that specifies the number of grid points in
C  the second dimension for U, SD_U, V, SD_V, DIV, U_ADJ, and V_ADJ.

C  DELX is a real variable that specifies the grid spacing in the first
C  dimension for U, SD_U, V, SD_V, DIV, U_ADJ, and V_ADJ.

C  DELY is a real variable that specifies the grid spacing in the second
C  dimension for U, SD_U, V, SD_V, DIV, U_ADJ, and V_ADJ.

C  BADDATA is a real variable that indicates a missing value as
C  described.

C  MAXX_OUT is an integer variable that specifies the first dimension of
C  U_ADJ and V_ADJ in the calling program.

C  MAXY_OUT is an integer variable that specifies the second dimension
C  of U_ADJ and V_ADJ in the calling program.

C  Output:

C  U_ADJ is a two-dimensional real array.  U_ADJ(I,J) returns the value,
C  at the Ith grid point in the first dimension and the Jth grid point
C  in the second dimension, of the adjusted wind component in the first
C  dimension.  If it is missing, it is returned as BADDATA.  U_ADJ may
C  be identical to U, in which case U will be overwritten.

C  V_ADJ is a two-dimensional real array.  V_ADJ(I,J) returns the value,
C  at the Ith grid point in the first dimension and the Jth grid point
C  in the second dimension, of the adjusted wind component in the second
C  dimension.  If it is missing, it is returned as BADDATA.  V_ADJ may
C  be identical to V, in which case V will be overwritten.

      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      INTEGER MAXX_IN,MAXY_IN,MAXX_OUT,MAXY_OUT,NX,NY,IX,IY,IR,
     $N_ITERATIONS,IL,IK,IS,IM,IN
      REAL DELX,DELY,BADDATA,RELAX,MAX_ADJUST,SAVE,A,B,C
      REAL U(MAXX_IN,MAXY_IN),SD_U(MAXX_IN,MAXY_IN),V(MAXX_IN,MAXY_IN),
     $SD_V(MAXX_IN,MAXY_IN),DIV(MAXX_IN,MAXY_IN)
      REAL U_ADJ(MAXX_OUT,MAXY_OUT),V_ADJ(MAXX_OUT,MAXY_OUT)
      REAL LAMBDA(NX,NY),RHS(NX,NY),U_TEMP(NX,NY),V_TEMP(NX,NY)
      REAL WT_U(NX,NY,-1:1),WT_V(NX,NY,-1:1)

C  Copy the input fields in case they are to be overwritten.
      DO 1 IY=1,NY
         DO 2 IX=1,NX
            U_TEMP(IX,IY)=U(IX,IY)
            V_TEMP(IX,IY)=V(IX,IY)
2        CONTINUE
1     CONTINUE

C  Set the overrelaxation coefficient.
      RELAX=1.55

C  Initialize the weights to zero.
      DO 3 IY=1,NY
         DO 4 IX=1,NX
            DO 5 IR=-1,1
               WT_U(IX,IY,IR)=0.
               WT_V(IX,IY,IR)=0.
5           CONTINUE
4        CONTINUE
3     CONTINUE

C  Calculate the U weights for the plane.
      DO 6 IY=1,NY
         DO 7 IX=1,NX
            IF(DIV(IX,IY).NE.BADDATA)THEN
               IF(IX.GT.1.AND.
     $         U_TEMP(IX-1,IY).NE.BADDATA)THEN
                  IF(IX.LT.NX.AND.
     $            U_TEMP(IX+1,IY).NE.BADDATA)THEN
                     WT_U(IX,IY,-1)=-1./(2.*DELX)
                     WT_U(IX,IY,1)=1./(2.*DELX)
                  ELSEIF(U_TEMP(IX,IY).NE.BADDATA)THEN
                     WT_U(IX,IY,-1)=-1./(DELX)
                     WT_U(IX,IY,0)=1./(DELX)
                  ENDIF
               ELSEIF(IX.LT.NX.AND.
     $         U_TEMP(IX+1,IY).NE.BADDATA.AND.
     $         U_TEMP(IX,IY).NE.BADDATA)THEN
                  WT_U(IX,IY,0)=-1./(DELX)
                  WT_U(IX,IY,1)=1./(DELX)
               ENDIF
            ENDIF
7        CONTINUE
6     CONTINUE

C  Calculate the V weights for the plane.
      DO 8 IY=1,NY
         DO 9 IX=1,NX
            IF(DIV(IX,IY).NE.BADDATA)THEN
               IF(IY.GT.1.AND.
     $         V_TEMP(IX,IY-1).NE.BADDATA)THEN
                  IF(IY.LT.NY.AND.
     $            V_TEMP(IX,IY+1).NE.BADDATA)THEN
                     WT_V(IX,IY,-1)=-1./(2.*DELY)
                     WT_V(IX,IY,1)=1./(2.*DELY)
                  ELSEIF(V_TEMP(IX,IY).NE.BADDATA)THEN
                     WT_V(IX,IY,-1)=-1./(DELY)
                     WT_V(IX,IY,0)=1./(DELY)
                  ENDIF
               ELSEIF(IY.LT.NY.AND.
     $         V_TEMP(IX,IY+1).NE.BADDATA.AND.
     $         V_TEMP(IX,IY).NE.BADDATA)THEN
                  WT_V(IX,IY,0)=-1./(DELY)
                  WT_V(IX,IY,1)=1./(DELY)
               ENDIF
            ENDIF
9        CONTINUE
8     CONTINUE

C  Calculate the right hand side for the plane.
      DO 10 IL=1,NY
         DO 11 IK=1,NX
            IF(DIV(IK,IL).NE.BADDATA)THEN
               RHS(IK,IL)=0.
               DO 12 IR=-1,1
                  IF(IK+IR.GE.1.AND.
     $            IK+IR.LE.NX)THEN
                     RHS(IK,IL)=RHS(IK,IL)+2.*WT_U(IK,IL,IR)*
     $               U_TEMP(IK+IR,IL)
                  ENDIF
12             CONTINUE
               DO 13 IS=-1,1
                  IF(IL+IS.GE.1.AND.
     $            IL+IS.LE.NY)THEN
                     RHS(IK,IL)=RHS(IK,IL)+2.*WT_V(IK,IL,IS)*
     $               V_TEMP(IK,IL+IS)
                  ENDIF
13             CONTINUE
               RHS(IK,IL)=RHS(IK,IL)-2*DIV(IK,IL)
            ELSE
               RHS(IK,IL)=BADDATA
            ENDIF
11       CONTINUE
10    CONTINUE

C  Initialize the output wind component fields to the input wind
C  component fields and the Lagrange multipliers to zero.
      DO 14 IY=1,NY
         DO 15 IX=1,NX
            U_ADJ(IX,IY)=U_TEMP(IX,IY)
            V_ADJ(IX,IY)=V_TEMP(IX,IY)
            LAMBDA(IX,IY)=0.
15       CONTINUE
14    CONTINUE

C  Iterate to a solution.
      N_ITERATIONS=0
      MAX_ADJUST=2.*ADJUST_WINDS_TO_DIV_THRESH
      DOWHILE(N_ITERATIONS.LT.ADJUST_WINDS_TO_DIV_MAX_ITERATIONS.AND.
     $MAX_ADJUST.GT.ADJUST_WINDS_TO_DIV_THRESH)   
         N_ITERATIONS=N_ITERATIONS+1
         MAX_ADJUST=0.
         DO 16 IL=1,NY
            DO 17 IK=1,NX
               IF(U_TEMP(IK,IL).NE.BADDATA.AND.
     $         SD_U(IK,IL).NE.BADDATA)THEN
                  SAVE=U_ADJ(IK,IL)
                  C=0.
                  DO 18 IM=-1,1
                     IF(IK-IM.GE.1.AND.
     $               IK-IM.LE.NX)THEN
                        C=C+LAMBDA(IK-IM,IL)*WT_U(IK-IM,IL,IM)
                     ENDIF
18                CONTINUE
                  U_ADJ(IK,IL)=U_TEMP(IK,IL)-C*SD_U(IK,IL)**2/2.
                  U_ADJ(IK,IL)=SAVE+(U_ADJ(IK,IL)-SAVE)*RELAX
                  MAX_ADJUST=AMAX1(MAX_ADJUST,ABS(U_ADJ(IK,IL)-SAVE))
               ENDIF
               IF(V_TEMP(IK,IL).NE.BADDATA.AND.
     $         SD_V(IK,IL).NE.BADDATA)THEN
                  SAVE=V_ADJ(IK,IL)
                  C=0.
                  DO 19 IN=-1,1
                     IF(IL-IN.GE.1.AND.
     $               IL-IN.LE.NY)THEN
                        C=C+LAMBDA(IK,IL-IN)*WT_V(IK,IL-IN,IN)
                     ENDIF
19                CONTINUE
                  V_ADJ(IK,IL)=V_TEMP(IK,IL)-C*SD_V(IK,IL)**2/2.
                  V_ADJ(IK,IL)=SAVE+(V_ADJ(IK,IL)-SAVE)*RELAX
                  MAX_ADJUST=AMAX1(MAX_ADJUST,ABS(V_ADJ(IK,IL)-SAVE))
               ENDIF
               IF(DIV(IK,IL).NE.BADDATA)THEN
                  A=0.
                  B=0.
                  DO 20 IR=-1,1
                     DO 21 IM=-1,1
                        IF(IK+IR.GE.1.AND.
     $                  IK+IR.LE.NX.AND.
     $                  IK+IR-IM.GE.1.AND.
     $                  IK+IR-IM.LE.NX)THEN
                           IF(IR-IM.EQ.0)THEN
                              A=A+WT_U(IK,IL,IR)*WT_U(IK+IR-IM,IL,IM)*
     $                        SD_U(IK+IR,IL)**2
                           ELSE
                              B=B+LAMBDA(IK+IR-IM,IL)*WT_U(IK,IL,IR)*
     $                        WT_U(IK+IR-IM,IL,IM)*SD_U(IK+IR,IL)**2
                           ENDIF
                        ENDIF
21                   CONTINUE
20                CONTINUE
                  DO 22 IS=-1,1
                     DO 23 IN=-1,1
                        IF(IL+IS.GE.1.AND.
     $                  IL+IS.LE.NY.AND.
     $                  IL+IS-IN.GE.1.AND.
     $                  IL+IS-IN.LE.NY)THEN
                           IF(IS-IN.EQ.0)THEN
                              A=A+WT_V(IK,IL,IS)*WT_V(IK,IL+IS-IN,IN)*
     $                        SD_V(IK,IL+IS)**2
                           ELSE
                              B=B+LAMBDA(IK,IL+IS-IN)*WT_V(IK,IL,IS)*
     $                        WT_V(IK,IL+IS-IN,IN)*SD_V(IK,IL+IS)**2
                           ENDIF
                        ENDIF
23                   CONTINUE
22                CONTINUE
                  IF(A.NE.0.)THEN
                     LAMBDA(IK,IL)=LAMBDA(IK,IL)+
     $               ((RHS(IK,IL)-B)/A-LAMBDA(IK,IL))*RELAX
                  ENDIF
               ENDIF
17          CONTINUE
16       CONTINUE
      ENDDO

C  Report if the solution did not converge within the specified number
C  of iterations.
      IF(N_ITERATIONS.GE.ADJUST_WINDS_TO_DIV_MAX_ITERATIONS)THEN
         WRITE(TMMLIB_MESSAGE_UNIT,*)'ADJUST_WINDS_TO_DIV_2D:  ',
     $   'SOLUTION DID NOT CONVERGE AFTER ',
     $   ADJUST_WINDS_TO_DIV_MAX_ITERATIONS,' ITERATIONS.  MAX_ADJUST ',
     $   '= ',MAX_ADJUST
      ENDIF

C  Done.
      RETURN
      END
