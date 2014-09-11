      SUBROUTINE DIVERGENCE_FILTER_2D_DRIVER3(U,SD_U,V,SD_V,MAXX_IN,
     $MAXY_IN,MAXZ_IN,NX,NY,NZ,DELX,DELY,BADDATA,WT_X,MAX_REACH,REACH_X,
     $WT_Y,REACH_Y,MAXX_OUT,MAXY_OUT,MAXZ_OUT,DIV,SD_DIV)

C  Thomas Matejka NOAA/NSSL 8 March 1994

C  This subroutine is a third-dimensional driver for the subroutine
C  DIVERGENCE_FILTER_2D, which is applied in the first two dimensions.

C  Input:

C  U is a three-dimensional real array.  U(I,J,K) specifies the value,
C  at the Ith grid point in the first dimension, the Jth grid point in
C  the second dimension, and the Kth grid point in the third dimension,
C  of the velocity component field in the first dimension.  If it is
C  missing, it should equal BADDATA.  Errors in the values in U are
C  assumed to be uncorrelated with each other and with values in V.

C  SD_U is a three-dimensional real array.  SD_U(I,J,K) specifies the
C  standard deviation of U(I,J,K).  If it is missing, it should equal
C  BADDATA.  SD_U(I,J,K) must not be missing for U(I,J,K) to be usable.

C  V is a three-dimensional real array.  V(I,J,K) specifies the value,
C  at the Ith grid point in the first dimension, the Jth grid point in
C  the second dimension, and the Kth grid point in the third dimension,
C  of the velocity component field in the second dimension.  If it is
C  missing, it should equal BADDATA.  Errors in the values in V are
C  assumed to be uncorrelated with each other and with values in U.

C  SD_V is a three-dimensional real array.  SD_V(I,J,K) specifies the
C  standard deviation of V(I,J,K).  If it is missing, it should equal
C  BADDATA.  SD_V(I,J,K) must not be missing for V(I,J,K) to be usable.

C  MAXX_IN is an integer variable that specifies the first dimension of
C  U, SD_U, V, and SD_V in the calling program.

C  MAXY_IN is an integer variable that specifies the second dimension of
C  U, SD_U, V, and SD_V in the calling program.

C  MAXZ_IN is an integer variable that specifies the third dimension of
C  U, SD_U, V, and SD_V in the calling program.

C  NX is an integer variable that specifies the number of grid points in
C  the first dimension for U, SD_U, V, SD_V, DIV, and SD_DIV.

C  NY is an integer variable that specifies the number of grid points in
C  the second dimension for U, SD_U, V, SD_V, DIV, and SD_DIV.

C  NZ is an integer variable that specifies the number of grid points in
C  the third dimension for U, SD_U, V, SD_V, DIV, and SD_DIV.

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

C  MAXZ_OUT is an integer variable that specifies the third dimension of
C  DIV and SD_DIV in the calling program.

C  Output:

C  DIV is a three-dimensional real array.  DIV(I,J,K) returns the value,
C  at the Ith grid point in the first dimension, the Jth grid point in
C  the second dimension, and the Kth grid point in the third dimension,
C  of the divergence in the first two dimensions.  If it is missing, it
C  is returned as BADDATA.

C  SD_DIV is a three-dimensional real array.  SD_DIV(I,J,K) returns the
C  standard deviation of DIV(I,J,K).  If it is missing, it is returned
C  as BADDATA.

      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      INTEGER MAXX_IN,MAXY_IN,MAXZ_IN,MAXX_OUT,MAXY_OUT,MAXZ_OUT,NX,NY,
     $NZ,REACH_X,REACH_Y,MAX_REACH,IZ
      REAL DELX,DELY,BADDATA
      REAL WT_X(-MAX_REACH:MAX_REACH),WT_Y(-MAX_REACH:MAX_REACH)
      REAL U(MAXX_IN,MAXY_IN,MAXZ_IN),SD_U(MAXX_IN,MAXY_IN,MAXZ_IN),
     $V(MAXX_IN,MAXY_IN,MAXZ_IN),SD_V(MAXX_IN,MAXY_IN,MAXZ_IN)
      REAL DIV(MAXX_OUT,MAXY_OUT,MAXZ_OUT),
     $SD_DIV(MAXX_OUT,MAXY_OUT,MAXZ_OUT)

C  Loop through the planes and calculate the filtered divergence.
      DO 1 IZ=1,NZ
         CALL DIVERGENCE_FILTER_2D(U(1,1,IZ),SD_U(1,1,IZ),V(1,1,IZ),
     $   SD_V(1,1,IZ),MAXX_IN,MAXY_IN,NX,NY,DELX,DELY,BADDATA,WT_X,
     $   MAX_REACH,REACH_X,WT_Y,REACH_Y,MAXX_OUT,MAXY_OUT,DIV(1,1,IZ),
     $   SD_DIV(1,1,IZ))
1     CONTINUE

C  Done.
      RETURN
      END
