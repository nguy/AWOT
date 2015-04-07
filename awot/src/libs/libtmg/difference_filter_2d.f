      SUBROUTINE DIFFERENCE_FILTER_2D(A,SD_A,MAXX,MAXY,NX,NY,BADDATA,
     $WT_XY,MAX_REACH_X,MAX_REACH_Y,REACH_X,REACH_Y,IX1,IY1,IX2,IY2,
     $A_DIFF,SD_A_DIFF)

C  Thomas Matejka NOAA/NSSL 26 April 1994

C  This subroutine calculates the difference between the values at two
C  specified points in a two-dimensional data field.  The data field is
C  filtered in two dimensions during the computation of the difference
C  but is not changed.  A field of the standard deviation of the
C  difference is also produced under the assumption that errors in the
C  input data field are uncorrelated and taking into account the effect
C  of the filter.

C  Input:

C  A is a two-dimensional real array.  A(I,J) specifies the value, at
C  the Ith grid point in the first dimension and the Jth grid point in
C  the second dimension, of the data field for which a difference is
C  sought.  If it is missing, it should equal BADDATA.  Errors in the
C  values in A are assumed to be uncorrelated.

C  SD_A is a two-dimensional real array.  SD_A(I,J) specifies the
C  standard deviation of A(I,J).  If it is missing, it should equal
C  BADDATA.  SD_A(I,J) must not be missing for A(I,J) to be usable.

C  MAXX is an integer variable that specifies the first dimension of A
C  and SD_A in the calling program.

C  MAXY is an integer variable that specifies the second dimension of A
C  and SD_A in the calling program.

C  NX is an integer variable that specifies the number of grid points in
C  the first dimension for A and SD_A.

C  NY is an integer variable that specifies the number of grid points in
C  the second dimension for A and SD_A.

C  BADDATA is a real variable that indicates a missing value as
C  described.

C  WT_XY is a two-dimensional real array.  WT_XY(I,J) specifies the
C  weight at the Ith grid point in the first dimension and the Jth grid
C  point in the second dimension relative to the point of application
C  for the two-dimensional filter.  The weights in WT_XY can be
C  specified to within a constant factor, since the filter weights are
C  normalized in the subroutine.

C  MAX_REACH_X is an integer variable that specifies that WT_XY is
C  dimensioned from -MAX_REACH_X to MAX_REACH_X in the first dimension
C  in the calling program.

C  MAX_REACH_Y is an integer variable that specifies that WT_YY is
C  dimensioned from -MAX_REACH_Y to MAX_REACH_Y in the second dimension
C  in the calling program.

C  REACH_X is an integer variable that specifies that weights in WT_XY
C  run from -REACH_X to REACH_X grid points in the first dimension
C  relative to the point of application.

C  REACH_Y is an integer variable that specifies that weights in WT_XY
C  run from -REACH_Y to REACH_Y grid points in the second dimension
C  relative to the point of application.

C  IX1 is an integer variable that specifies the grid point number in
C  the first dimension of the point involved in the subtraction
C  A(IX1,IY1) - A(IX2,IY2).

C  IY1 is an integer variable that specifies the grid point number in
C  the second dimension of the point involved in the subtraction
C  A(IX1,IY1) - A(IX2,IY2).

C  IX2 is an integer variable that specifies the grid point number in
C  the first dimension of the point involved in the subtraction
C  A(IX1,IY1) - A(IX2,IY2).

C  IY2 is an integer variable that specifies the grid point number in
C  the second dimension of the point involved in the subtraction
C  A(IX1,IY1) - A(IX2,IY2).

C  Output:

C  A_DIFF is a real variable that returns the difference A(IX1,IY1) -
C  A(IX2,IY2).  If it is missing, it is returned as BADDATA.

C  SD_A_DIFF is a real variable that returns the standard deviation of
C  A_DIFF.  If it is missing, it is returned as BADDATA.

      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      INTEGER MAXX,MAXY,NX,NY,REACH_X,REACH_Y,MAX_REACH_X,MAX_REACH_Y,P,
     $Q,PMIN,PMAX,QMIN,QMAX,IX1,IY1,IX2,IY2,JY,JX,KX1,KY1,KX2,KY2
      REAL BADDATA,A_DIFF,SD_A_DIFF,SUM1,SUM2,SUM_VAR1,SUM_VAR2,SUM_COV,
     $SUM_WT1,SUM_WT2,VAR_A_DIFF
      REAL WT_XY(-MAX_REACH_X:MAX_REACH_X,-MAX_REACH_Y:MAX_REACH_Y)
      REAL A(MAXX,MAXY),SD_A(MAXX,MAXY)

C  Check whether values exist at the points of interest.
      IF(A(IX1,IY1).NE.BADDATA.AND.
     $SD_A(IX1,IY1).NE.BADDATA.AND.
     $A(IX2,IY2).NE.BADDATA.AND.
     $SD_A(IX2,IY2).NE.BADDATA)THEN

C  Loop through the points within the reach of the filter and calculate
C  sums.
         SUM1=0.
         SUM2=0.
         SUM_VAR1=0.
         SUM_VAR2=0.
         SUM_WT1=0.
         SUM_WT2=0.
         DO 1 JY=-REACH_Y,REACH_Y
            KY1=IY1+JY
            KY2=IY2+JY
            DO 2 JX=-REACH_X,REACH_X
               KX1=IX1+JX
               KX2=IX2+JX
               IF(KX1.GE.1.AND.
     $         KX1.LE.NX.AND.
     $         KY1.GE.1.AND.
     $         KY1.LE.NY)THEN
                  IF(A(KX1,KY1).NE.BADDATA.AND.
     $            SD_A(KX1,KY1).NE.BADDATA)THEN
                     SUM1=SUM1+A(KX1,KY1)*WT_XY(JX,JY)
                     SUM_VAR1=SUM_VAR1+(SD_A(KX1,KY1)*WT_XY(JX,JY))**2
                     SUM_WT1=SUM_WT1+WT_XY(JX,JY)
                  ENDIF
               ENDIF
               IF(KX2.GE.1.AND.
     $         KX2.LE.NX.AND.
     $         KY2.GE.1.AND.
     $         KY2.LE.NY)THEN
                  IF(A(KX2,KY2).NE.BADDATA.AND.
     $            SD_A(KX2,KY2).NE.BADDATA)THEN
                     SUM2=SUM2+A(KX2,KY2)*WT_XY(JX,JY)
                     SUM_VAR2=SUM_VAR2+(SD_A(KX2,KY2)*WT_XY(JX,JY))**2
                     SUM_WT2=SUM_WT2+WT_XY(JX,JY)
                  ENDIF
               ENDIF
2           CONTINUE
1        CONTINUE

C  Loop through points common to the filtered values at each point and
C  calculate a sum for the covariance.
         SUM_COV=0.
         PMIN=MAX0(IX1-REACH_X,IX2-REACH_X,1)
         PMAX=MIN0(IX1+REACH_X,IX2+REACH_X,NX)
         QMIN=MAX0(IY1-REACH_Y,IY2-REACH_Y,1)
         QMAX=MIN0(IY1+REACH_Y,IY2+REACH_Y,NY)
         IF(PMAX.GT.PMIN.AND.
     $   QMAX.GT.QMIN)THEN
            DO 3 P=PMIN,PMAX
               DO 4 Q=QMIN,QMAX
                  IF(A(P,Q).NE.BADDATA.AND.
     $            SD_A(P,Q).NE.BADDATA)THEN
                     SUM_COV=SUM_COV+SD_A(P,Q)**2*WT_XY(P-IX1,Q-IY1)*
     $               WT_XY(P-IX2,Q-IY2)
                  ENDIF
4              CONTINUE
3           CONTINUE
         ENDIF

C  Calculate the filtered difference and its standard deviation.
         IF(SUM_WT1.NE.0..AND.
     $   SUM_WT2.NE.0.)THEN
            A_DIFF=SUM1/SUM_WT1-SUM2/SUM_WT2
            VAR_A_DIFF=SUM_VAR1/SUM_WT1**2+SUM_VAR2/SUM_WT2**2-
     $      2.*SUM_COV/(SUM_WT1*SUM_WT2)
            IF(VAR_A_DIFF.GT.0.)THEN
               SD_A_DIFF=SQRT(VAR_A_DIFF)
            ELSE
               SD_A_DIFF=0.
            ENDIF
         ELSE
            A_DIFF=BADDATA
            SD_A_DIFF=BADDATA
         ENDIF
      ELSE
         A_DIFF=BADDATA
         SD_A_DIFF=BADDATA
      ENDIF

C  Done.
      RETURN
      END
