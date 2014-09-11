      SUBROUTINE ADJUST_WINDS_2D(U,SD_U,V,SD_V,DIV,
     $MAXX_IN,MAXY_IN,DELX,DELY,NX,NY,
     $BADDATA,
     $ADJUST_THRESH,MAX_ITERATIONS,
     $MAXX_OUT,MAXY_OUT,
     $U_ADJ,V_ADJ,MAX_ADJUST,N_ITERATIONS,SUCCESS)

C  Thomas Matejka NOAA/NSSL 8 April 1998

C  This subroutine variationally adjusts two velocity component fields
C  so that they agree with a specified two-dimensional divergence field.
C  The solution is formulated so that the adjusted velocity component
C  fields are as close as possible to the input velocity component
C  fields.

C  Input:

C  U (2d real array 1:MAXX_IN,1:MAXY_IN).  U(I,J) specifies the value,
C  at the Ith grid point in the first dimension and the Jth grid point
C  in the second dimension, of the velocity component field in the first
C  dimension to be adjusted.  If it is missing, it should equal BADDATA.

C  SD_U (2d real array 1:MAXX_IN,1:MAXY_IN).  SD_U(I,J) specifies the
C  rms error of U(I,J).  If it is missing, it should equal BADDATA.
C  SD_U should not be missing for U(I,J) to be usable.

C  V (2d real array 1:MAXX_IN,1:MAXY_IN).  V(I,J) specifies the value,
C  at the Ith grid point in the first dimension and the Jth grid point
C  in the second dimension, of the velocity component field in the
C  second dimension to be adjusted.  If it is missing, it should equal
C  BADDATA.

C  SD_V (2d real array 1:MAXX_IN,1:MAXY_IN).  SD_V(I,J) specifies the
C  rms error of V(I,J).  If it is missing, it should equal BADDATA.
C  SD_V should not be missing for V(I,J) to be usable.

C  DIV (2d real array 1:MAXX_IN,1:MAXY_IN).  DIV(I,J) specifies the
C  two-dimensional divergence of the wind at the Ith grid point in the
C  first dimension and the Jth grid point in the second dimension.  If
C  it is missing, it should equal BADDATA.

C  MAXX_IN (integer) specifies the first dimension of U, SD_U, V, SD_V,
C  and DIV in the calling program.

C  MAXY_IN (integer) specifies the second dimension of U, SD_U, V, SD_V,
C  and DIV in the calling program.

C  NX (integer) specifies the number of grid points in the first
C  dimension for U, SD_U, V, SD_V, DIV, U_ADJ, and V_ADJ.

C  NY (integer) specifies the number of grid points in the second
C  dimension for U, SD_U, V, SD_V, DIV, U_ADJ, and V_ADJ.

C  DELX (real) specifies the grid spacing in the first dimension for U,
C  SD_U, V, SD_V, DIV, U_ADJ, and V_ADJ.

C  DELY (real) specifies the grid spacing in the second dimension for U,
C  SD_U, V, SD_V, DIV, U_ADJ, and V_ADJ.

C  BADDATA (real) indicates a missing value as described.

C  ADJUST_THRESH (real) controls whether the solution has converged.  If
C  the maximum change in a velocity component from one iteration to the
C  next is less than ADJUST_THRESH, the solution has converged.

C  MAX_ITERATIONS (integer) specifies the maximum number of iterations
C  to a solution allowed.

C  MAXX_OUT (integer) specifies the first dimension of U_ADJ and V_ADJ
C  in the calling program.

C  MAXY_OUT (integer) specifies the second dimension of U_ADJ and V_ADJ
C  in the calling program.

C  Output:

C  U_ADJ (2d real array 1:MAXX_OUT,1:MAXY_OUT).  U_ADJ(I,J) returns the
C  value, at the Ith grid point in the first dimension and the Jth grid
C  point in the second dimension, of the adjusted wind component in the
C  first dimension.  If it is missing, it is returned as BADDATA.  U_ADJ
C  may be identical to U, in which case U will be overwritten.

C  V_ADJ (2d real array 1:MAXX_OUT,1:MAXY_OUT).  V_ADJ(I,J) returns the
C  value, at the Ith grid point in the first dimension and the Jth grid
C  point in the second dimension, of the adjusted wind component in the
C  second dimension.  If it is missing, it is returned as BADDATA.
C  V_ADJ may be identical to V, in which case V will be overwritten.

      IMPLICIT NONE
      REAL,PARAMETER::RELAX=1.55
      LOGICAL::SUCCESS
      INTEGER::MAXX_IN,MAXY_IN,MAXX_OUT,MAXY_OUT,NX,NY,IX,IY,IR,
     $N_ITERATIONS,IL,IK,IS,IM,IN,MAX_ITERATIONS
      REAL::DELX,DELY,BADDATA,MAX_ADJUST,SAVE,A,B,C,ADJUST_THRESH
      REAL,DIMENSION(MAXX_IN,MAXY_IN)::U,SD_U,V,SD_V,DIV
      REAL,DIMENSION(MAXX_OUT,MAXY_OUT)::U_ADJ,V_ADJ
      REAL,DIMENSION(NX,NY)::LAMBDA,RHS,U_TEMP,V_TEMP
      REAL,DIMENSION(NX,NY,-1:1)::WT_U,WT_V

C  Copy the input fields in case they are to be overwritten.
      DO IY=1,NY
         DO IX=1,NX
            U_TEMP(IX,IY)=U(IX,IY)
            V_TEMP(IX,IY)=V(IX,IY)
         ENDDO
      ENDDO

C  Initialize the weights to zero.
      DO IY=1,NY
         DO IX=1,NX
            DO IR=-1,1
               WT_U(IX,IY,IR)=0.
               WT_V(IX,IY,IR)=0.
            ENDDO
         ENDDO
      ENDDO

C  Calculate the U weights for the plane.
      DO IY=1,NY
         DO IX=1,NX
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
         ENDDO
      ENDDO

C  Calculate the V weights for the plane.
      DO IY=1,NY
         DO IX=1,NX
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
         ENDDO
      ENDDO

C  Calculate the right hand side for the plane.
      DO IL=1,NY
         DO IK=1,NX
            IF(DIV(IK,IL).NE.BADDATA)THEN
               RHS(IK,IL)=0.
               DO IR=-1,1
                  IF(IK+IR.GE.1.AND.
     $            IK+IR.LE.NX)THEN
                     RHS(IK,IL)=RHS(IK,IL)+2.*WT_U(IK,IL,IR)*
     $               U_TEMP(IK+IR,IL)
                  ENDIF
               ENDDO
               DO IS=-1,1
                  IF(IL+IS.GE.1.AND.
     $            IL+IS.LE.NY)THEN
                     RHS(IK,IL)=RHS(IK,IL)+2.*WT_V(IK,IL,IS)*
     $               V_TEMP(IK,IL+IS)
                  ENDIF
               ENDDO
               RHS(IK,IL)=RHS(IK,IL)-2*DIV(IK,IL)
            ELSE
               RHS(IK,IL)=BADDATA
            ENDIF
         ENDDO
      ENDDO

C  Initialize the output wind component fields to the input wind
C  component fields and the Lagrange multipliers to zero.
      DO IY=1,NY
         DO IX=1,NX
            U_ADJ(IX,IY)=U_TEMP(IX,IY)
            V_ADJ(IX,IY)=V_TEMP(IX,IY)
            LAMBDA(IX,IY)=0.
         ENDDO
      ENDDO

C  Iterate to a solution.
      N_ITERATIONS=0
      MAX_ADJUST=2.*ADJUST_THRESH
      DO
         IF(N_ITERATIONS.GE.MAX_ITERATIONS.OR.
     $   MAX_ADJUST.LE.ADJUST_THRESH)THEN
            EXIT
         ENDIF
         N_ITERATIONS=N_ITERATIONS+1
         MAX_ADJUST=0.
         DO IL=1,NY
            DO IK=1,NX
               IF(U_TEMP(IK,IL).NE.BADDATA.AND.
     $         SD_U(IK,IL).NE.BADDATA)THEN
                  SAVE=U_ADJ(IK,IL)
                  C=0.
                  DO IM=-1,1
                     IF(IK-IM.GE.1.AND.
     $               IK-IM.LE.NX)THEN
                        C=C+LAMBDA(IK-IM,IL)*WT_U(IK-IM,IL,IM)
                     ENDIF
                  ENDDO
                  U_ADJ(IK,IL)=U_TEMP(IK,IL)-
     $            C*SD_U(IK,IL)*SD_U(IK,IL)/2.
                  U_ADJ(IK,IL)=SAVE+(U_ADJ(IK,IL)-SAVE)*RELAX
                  MAX_ADJUST=AMAX1(MAX_ADJUST,ABS(U_ADJ(IK,IL)-SAVE))
               ENDIF
               IF(V_TEMP(IK,IL).NE.BADDATA.AND.
     $         SD_V(IK,IL).NE.BADDATA)THEN
                  SAVE=V_ADJ(IK,IL)
                  C=0.
                  DO IN=-1,1
                     IF(IL-IN.GE.1.AND.
     $               IL-IN.LE.NY)THEN
                        C=C+LAMBDA(IK,IL-IN)*WT_V(IK,IL-IN,IN)
                     ENDIF
                  ENDDO
                  V_ADJ(IK,IL)=V_TEMP(IK,IL)-
     $            C*SD_V(IK,IL)*SD_V(IK,IL)/2.
                  V_ADJ(IK,IL)=SAVE+(V_ADJ(IK,IL)-SAVE)*RELAX
                  MAX_ADJUST=AMAX1(MAX_ADJUST,ABS(V_ADJ(IK,IL)-SAVE))
               ENDIF
               IF(DIV(IK,IL).NE.BADDATA)THEN
                  A=0.
                  B=0.
                  DO IR=-1,1
                     DO IM=-1,1
                        IF(IK+IR.GE.1.AND.
     $                  IK+IR.LE.NX.AND.
     $                  IK+IR-IM.GE.1.AND.
     $                  IK+IR-IM.LE.NX)THEN
                           IF(IR-IM.EQ.0)THEN
                              A=A+WT_U(IK,IL,IR)*WT_U(IK+IR-IM,IL,IM)*
     $                        SD_U(IK+IR,IL)*SD_U(IK+IR,IL)
                           ELSE
                              B=B+LAMBDA(IK+IR-IM,IL)*WT_U(IK,IL,IR)*
     $                        WT_U(IK+IR-IM,IL,IM)*SD_U(IK+IR,IL)*
     $                        SD_U(IK+IR,IL)
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDDO
                  DO IS=-1,1
                     DO IN=-1,1
                        IF(IL+IS.GE.1.AND.
     $                  IL+IS.LE.NY.AND.
     $                  IL+IS-IN.GE.1.AND.
     $                  IL+IS-IN.LE.NY)THEN
                           IF(IS-IN.EQ.0)THEN
                              A=A+WT_V(IK,IL,IS)*WT_V(IK,IL+IS-IN,IN)*
     $                        SD_V(IK,IL+IS)*SD_V(IK,IL+IS)
                           ELSE
                              B=B+LAMBDA(IK,IL+IS-IN)*WT_V(IK,IL,IS)*
     $                        WT_V(IK,IL+IS-IN,IN)*SD_V(IK,IL+IS)*
     $                        SD_V(IK,IL+IS)
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDDO
                  IF(A.NE.0.)THEN
                     LAMBDA(IK,IL)=LAMBDA(IK,IL)+
     $               ((RHS(IK,IL)-B)/A-LAMBDA(IK,IL))*RELAX
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO

C  Check whether the solution converged within the specified number of
C  iterations.
      IF(N_ITERATIONS.GE.MAX_ITERATIONS)THEN
         SUCCESS=.FALSE.
      ELSE
         SUCCESS=.TRUE.
      ENDIF

      END SUBROUTINE ADJUST_WINDS_2D

