      SUBROUTINE DDY_FILTER_2D(A,SD_A,MAXX_IN,MAXY_IN,NX,NY,DELY,
     $BADDATA,WT_X,MAX_REACH,REACH_X,WT_Y,REACH_Y,MAXX_OUT,MAXY_OUT,
     $DDY_A,SD_DDY_A)

C  Thomas Matejka NOAA/NSSL 26 April 1994

C  This subroutine calculates the partial derivative with respect to the
C  second dimension of a two-dimensional data field.  A centered
C  difference over two grid spaces is used where possible.  At data
C  edges, a two-point difference over one grid space is used where
C  possible.  The data field is filtered in two dimensions during the
C  computation of the derivative field but is not changed.  A field of
C  the standard deviation of the derivative field is also produced under
C  the assumption that errors in the input data field are uncorrelated
C  and taking into account the effect of the filter.

C  Input:

C  A is a two-dimensional real array.  A(I,J) specifies the value, at
C  the Ith grid point in the first dimension and the Jth grid point in
C  the second dimension, of the data field whose derivative is sought.
C  If it is missing, it should equal BADDATA.  Errors in the values in A
C  are assumed to be uncorrelated.

C  SD_A is a two-dimensional real array.  SD_A(I,J) specifies the
C  standard deviation of A(I,J).  If it is missing, it should equal
C  BADDATA.  SD_A(I,J) must not be missing for A(I,J) to be usable.

C  MAXX_IN is an integer variable that specifies the first dimension of
C  A and SD_A in the calling program.

C  MAXY_IN is an integer variable that specifies the second dimension of
C  A and SD_A in the calling program.

C  NX is an integer variable that specifies the number of grid points in
C  the first dimension for A, SD_A, DDY_A, and SD_DDY_A.

C  NY is an integer variable that specifies the number of grid points in
C  the second dimension for A, SD_A, DDY_A, and SD_DDY_A.

C  DELY is a real variable that specifies the grid spacing in the second
C  dimension for A, SD_A, DDY_A, and SD_DDY_A.

C  BADDATA is a real variable that indicates a missing value as
C  described.

C  WT_X is a one-dimensional real array.  WT_X(I) specifies the weight
C  at the Ith grid point relative to the point of application for the
C  one-dimensional filter to be applied in the first dimension.  The
C  weights in WT_X and WT_Y can be specified to within a constant
C  factor, since the filter weights are normalized in the subroutine.

C  MAX_REACH is an integer variable that specifies that WT_X and WT_Y
C  are dimensioned from -MAX_REACH to MAX_REACH in the calling program.

C  REACH_X is an integer variable that specifies that weights in WT_X
C  run from -REACH_X to REACH_X grid points relative to the point of
C  application.

C  WT_Y is a one-dimensional real array.  WT_Y(I) specifies the weight
C  at the Ith grid point relative to the point of application for the
C  one-dimensional filter to be applied in the second dimension.  The
C  weights in WT_X and WT_Y can be specified to within a constant
C  factor, since the filter weights are normalized in the subroutine.

C  REACH_Y is an integer variable that specifies that weights in WT_Y
C  run from -REACH_Y to REACH_Y grid points relative to the point of
C  application.

C  MAXX_OUT is an integer variable that specifies the first dimension of
C  DDY_A and SD_DDY_A in the calling program.

C  MAXY_OUT is an integer variable that specifies the second dimension
C  of DDY_A and SD_DDY_A in the calling program.

C  Output:

C  DDY_A is a two-dimensional real array.  DDY_A(I,J) returns the value,
C  at the Ith grid point in the first dimension and the Jth grid point
C  in the second dimension, of the partial derivative of A with respect
C  to the second dimension.  If it is missing, it is returned as
C  BADDATA.

C  SD_DDY_A is a two-dimensional real array.  SD_DDY_A(I,J) returns the
C  standard deviation of DDY_A(I,J).  If it is missing, it is returned
C  as BADDATA.

      IMPLICIT NONE
      INTEGER MAXX_IN,MAXY_IN,MAXX_OUT,MAXY_OUT,NX,NY,REACH_X,REACH_Y,
     $MAX_REACH,IX,IY,JX,JY
      REAL DELY,BADDATA,A_DIFF,SD_A_DIFF
      REAL WT_X(-MAX_REACH:MAX_REACH),WT_Y(-MAX_REACH:MAX_REACH)
      REAL WT_XY(-REACH_X:REACH_X,-REACH_Y:REACH_Y)
      REAL A(MAXX_IN,MAXY_IN),SD_A(MAXX_IN,MAXY_IN)
      REAL DDY_A(MAXX_OUT,MAXY_OUT),SD_DDY_A(MAXX_OUT,MAXY_OUT)

C  Calculate the weights for the two-dimensional filter.
      DO 1 JY=-REACH_Y,REACH_Y
         DO 2 JX=-REACH_X,REACH_X
            WT_XY(JX,JY)=WT_X(JX)*WT_Y(JY)
2        CONTINUE
1     CONTINUE

C  Loop through the points on the plane and calculate the derivative of
C  the filtered data field.
      DO 3 IY=1,NY
         DO 4 IX=1,NX
            IF(IY.GT.1.AND.
     $      A(IX,IY-1).NE.BADDATA.AND.
     $      SD_A(IX,IY-1).NE.BADDATA)THEN
               IF(IY.LT.NY.AND.
     $         A(IX,IY+1).NE.BADDATA.AND.
     $         SD_A(IX,IY+1).NE.BADDATA)THEN
                  CALL DIFFERENCE_FILTER_2D(A,SD_A,MAXX_IN,MAXY_IN,NX,
     $            NY,BADDATA,WT_XY,REACH_X,REACH_Y,REACH_X,REACH_Y,IX,
     $            IY+1,IX,IY-1,A_DIFF,SD_A_DIFF)
                  DDY_A(IX,IY)=A_DIFF/(2.*DELY)
                  SD_DDY_A(IX,IY)=SD_A_DIFF/(2.*DELY)
               ELSEIF(A(IX,IY).NE.BADDATA.AND.
     $         SD_A(IX,IY).NE.BADDATA)THEN
                  CALL DIFFERENCE_FILTER_2D(A,SD_A,MAXX_IN,MAXY_IN,NX,
     $            NY,BADDATA,WT_XY,REACH_X,REACH_Y,REACH_X,REACH_Y,IX,
     $            IY,IX,IY-1,A_DIFF,SD_A_DIFF)
                  DDY_A(IX,IY)=A_DIFF/DELY
                  SD_DDY_A(IX,IY)=SD_A_DIFF/DELY
               ELSE
                  DDY_A(IX,IY)=BADDATA
                  SD_DDY_A(IX,IY)=BADDATA
               ENDIF
            ELSEIF(IY.LT.NY.AND.
     $      A(IX,IY+1).NE.BADDATA.AND.
     $      SD_A(IX,IY+1).NE.BADDATA.AND.
     $      A(IX,IY).NE.BADDATA.AND.
     $      SD_A(IX,IY).NE.BADDATA)THEN
               CALL DIFFERENCE_FILTER_2D(A,SD_A,MAXX_IN,MAXY_IN,NX,NY,
     $         BADDATA,WT_XY,REACH_X,REACH_Y,REACH_X,REACH_Y,IX,IY+1,IX,
     $         IY,A_DIFF,SD_A_DIFF)
               DDY_A(IX,IY)=A_DIFF/DELY
               SD_DDY_A(IX,IY)=SD_A_DIFF/DELY
            ELSE
               DDY_A(IX,IY)=BADDATA
               SD_DDY_A(IX,IY)=BADDATA
            ENDIF
4        CONTINUE
3     CONTINUE

C  Done.
      RETURN
      END
