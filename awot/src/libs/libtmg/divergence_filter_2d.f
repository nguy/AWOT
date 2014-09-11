      SUBROUTINE DIVERGENCE_FILTER_2D(U,SD_U,V,SD_V,MAXX_IN,MAXY_IN,NX,
     $NY,DELX,DELY,BADDATA,WT_X,MAX_REACH,REACH_X,WT_Y,REACH_Y,MAXX_OUT,
     $MAXY_OUT,DIV,SD_DIV)

C  Thomas Matejka NOAA/NSSL 22 April 1994

C  This subroutine calculates a two-dimensional divergence field from
C  velocity component fields in each dimension.  The velocity component
C  fields are filtered in two dimensions during the computation of the
C  divergence field but are not changed.  A field of the standard
C  deviation of the divergence field is also produced under the
C  assumption that errors in the input velocity component fields are
C  uncorrelated and taking into account the effect of the filter.

C  Input:

C  U is a two-dimensional real array.  U(I,J) specifies the value, at
C  the Ith grid point in the first dimension and the Jth grid point in
C  the second dimension, of the velocity component field in the first
C  dimension.  If it is missing, it should equal BADDATA.  Errors in the
C  values in U are assumed to be uncorrelated with each other and with
C  values in V.

C  SD_U is a two-dimensional real array.  SD_U(I,J) specifies the
C  standard deviation of U(I,J).  If it is missing, it should equal
C  BADDATA.  SD_U(I,J) must not be missing for U(I,J) to be usable.

C  V is a two-dimensional real array.  V(I,J) specifies the value, at
C  the Ith grid point in the first dimension and the Jth grid point in
C  the second dimension, of the velocity component field in the second
C  dimension.  If it is missing, it should equal BADDATA.  Errors in the
C  values in V are assumed to be uncorrelated with each other and with
C  values in U.

C  SD_V is a two-dimensional real array.  SD_V(I,J) specifies the
C  standard deviation of V(I,J).  If it is missing, it should equal
C  BADDATA.  SD_V(I,J) must not be missing for V(I,J) to be usable.

C  MAXX_IN is an integer variable that specifies the first dimension of
C  U, SD_U, V, and SD_V in the calling program.

C  MAXY_IN is an integer variable that specifies the second dimension of
C  U, SD_U, V, and SD_V in the calling program.

C  NX is an integer variable that specifies the number of grid points in
C  the first dimension for U, SD_U, V, SD_V, DIV, and SD_DIV.

C  NY is an integer variable that specifies the number of grid points in
C  the second dimension for U, SD_U, V, SD_V, DIV, and SD_DIV.

C  DELX is a real variable that specifies the grid spacing in the first
C  dimension for U, SD_U, V, SD_V, DIV, and SD_DIV.

C  DELY is a real variable that specifies the grid spacing in the second
C  dimension for U, SD_U, V, SD_V, DIV, and SD_DIV.

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
C  DIV and SD_DIV in the calling program.

C  MAXY_OUT is an integer variable that specifies the second dimension
C  of DIV and SD_DIV in the calling program.

C  Output:

C  DIV is a two-dimensional real array.  DIV(I,J) returns the value, at
C  the Ith grid point in the first dimension and the Jth grid point in
C  the second dimension, of the two-dimensional divergence.  If it is
C  missing, it is returned as BADDATA.

C  SD_DIV is a two-dimensional real array.  SD_DIV(I,J) returns the
C  standard deviation of DIV(I,J).  If it is missing, it is returned as
C  BADDATA.

      IMPLICIT NONE
      INTEGER MAXX_IN,MAXY_IN,MAXX_OUT,MAXY_OUT,NX,NY,REACH_X,REACH_Y,
     $MAX_REACH,IX,IY
      REAL DELX,DELY,BADDATA
      REAL WT_X(-MAX_REACH:MAX_REACH),WT_Y(-MAX_REACH:MAX_REACH)
      REAL U(MAXX_IN,MAXY_IN),SD_U(MAXX_IN,MAXY_IN),
     $V(MAXX_IN,MAXY_IN),SD_V(MAXX_IN,MAXY_IN)
      REAL DIV(MAXX_OUT,MAXY_OUT),SD_DIV(MAXX_OUT,MAXY_OUT)
      REAL DDX_U(NX,NY),SD_DDX_U(NX,NY),DDY_V(NX,NY),SD_DDY_V(NX,NY)

C  Calculate the derivative of filtered U with respect to the first
C  dimension.
      CALL DDX_FILTER_2D(U,SD_U,MAXX_IN,MAXY_IN,NX,NY,DELX,BADDATA,WT_X,
     $MAX_REACH,REACH_X,WT_Y,REACH_Y,NX,NY,DDX_U,SD_DDX_U)

C  Calculate the derivative of filtered V with respect to the second
C  dimension.
      CALL DDY_FILTER_2D(V,SD_V,MAXX_IN,MAXY_IN,NX,NY,DELY,BADDATA,WT_X,
     $MAX_REACH,REACH_X,WT_Y,REACH_Y,NX,NY,DDY_V,SD_DDY_V)

C  Loop through the points on the plane and calculate the divergence and
C  its standard error.
      DO 1 IY=1,NY
         DO 2 IX=1,NX
            IF(DDX_U(IX,IY).NE.BADDATA.AND.
     $      SD_DDX_U(IX,IY).NE.BADDATA.AND.
     $      DDY_V(IX,IY).NE.BADDATA.AND.
     $      SD_DDY_V(IX,IY).NE.BADDATA)THEN
               DIV(IX,IY)=DDX_U(IX,IY)+DDY_V(IX,IY)
               SD_DIV(IX,IY)=SQRT(SD_DDX_U(IX,IY)**2+SD_DDY_V(IX,IY)**2)
            ELSE
               DIV(IX,IY)=BADDATA
               SD_DIV(IX,IY)=BADDATA
            ENDIF
2        CONTINUE
1     CONTINUE

C  Done.
      RETURN
      END
