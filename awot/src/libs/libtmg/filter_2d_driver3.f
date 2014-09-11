      SUBROUTINE FILTER_2D_DRIVER3(A_IN,SD_A_IN,MAXX_IN,MAXY_IN,MAXZ_IN,
     $NX,NY,NZ,BADDATA,UNCORRELATED,WT_X,MAX_REACH,REACH_X,WT_Y,REACH_Y,
     $MAXX_OUT,MAXY_OUT,MAXZ_OUT,A_OUT,SD_A_OUT)

C  Thomas Matejka NOAA/NSSL 28 May 1999

C  This subroutine is a third-dimensional driver for subroutine
C  FILTER_2D, which is applied in the first two dimensions.

C  Input:

C  A_IN is a three-dimensional real array.  A_IN(I,J,K) specifies the
C  value, at the Ith grid point in the first dimension, the Jth grid
C  point in the second dimension, and the Kth grid point in the third
C  dimension, of the data field to be filtered.  If it is missing, it
C  should equal BADDATA.  Errors in the values in A_IN are assumed to be
C  uncorrelated.

C  SD_A_IN is a three-dimensional real array.  SD_A_IN(I,J,K) specifies
C  the standard deviation of A_IN(I,J,K).  If it is missing, it should
C  equal BADDATA.  SD_A_IN(I,J,K) must not be missing for A_IN(I,J,K) to
C  be usable.

C  MAXX_IN is an integer variable that specifies the first dimension of
C  A_IN and SD_A_IN in the calling program.

C  MAXY_IN is an integer variable that specifies the second dimension of
C  A_IN and SD_A_IN in the calling program.

C  MAXZ_IN is an integer variable that specifies the third dimension of
C  A_IN and SD_A_IN in the calling program.

C  NX is an integer variable that specifies the number of grid points in
C  the first dimension for A_IN, SD_A_IN, A_OUT, and SD_A_OUT.

C  NY is an integer variable that specifies the number of grid points in
C  the second dimension for A_IN, SD_A_IN, A_OUT, and SD_A_OUT.

C  NZ is an integer variable that specifies the number of grid points in
C  the third dimension for A_IN, SD_A_IN, A_OUT, and SD_A_OUT.

C  BADDATA is a real variable that indicates a missing value as
C  described.

C  UNCORRELATED is a logical variable.  If UNCORRELATED is .TRUE., then
C  the original data are assumed to be uncorrelated.  If UNCORRELATED is
C  .FALSE., then the original data are assumed to be fully correlated.

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
C  A_OUT and SD_A_OUT in the calling program.

C  MAXY_OUT is an integer variable that specifies the second dimension
C  of A_OUT and SD_A_OUT in the calling program.

C  MAXZ_OUT is an integer variable that specifies the third dimension of
C  A_OUT and SD_A_OUT in the calling program.

C  Output:

C  A_OUT is a three-dimensional real array.  A_OUT(I,J,K) returns the
C  value, at the Ith grid point in the first dimension, the Jth grid
C  point in the second dimension, and the third grid point in the third
C  dimension, of the data field that has been filtered.  If it is
C  missing, it is returned as BADDATA.  A_OUT may be identical to A_IN,
C  in which case A_IN will be overwritten.

C  SD_A_OUT is a three-dimensional real array.  SD_A_OUT(I,J,K) returns
C  the standard deviation of A_OUT(I,J,K).  If it is missing, it is
C  returned as BADDATA.  SD_A_OUT may be identical to SD_A_IN, in which
C  case SD_A_IN will be overwritten.

      IMPLICIT NONE
      LOGICAL UNCORRELATED
      INTEGER MAXX_IN,MAXY_IN,MAXZ_IN,MAXX_OUT,MAXY_OUT,MAXZ_OUT,NX,NY,
     $NZ,REACH_X,REACH_Y,IZ,MAX_REACH
      REAL BADDATA
      REAL WT_X(-MAX_REACH:MAX_REACH),WT_Y(-MAX_REACH:MAX_REACH)
      REAL A_IN(MAXX_IN,MAXY_IN,MAXZ_IN),
     $SD_A_IN(MAXX_IN,MAXY_IN,MAXZ_IN)
      REAL A_OUT(MAXX_OUT,MAXY_OUT,MAXZ_OUT),
     $SD_A_OUT(MAXX_OUT,MAXY_OUT,MAXZ_OUT)

C  Loop through the planes and calculate the filtered field values.
      DO IZ=1,NZ
         CALL FILTER_2D(A_IN(1,1,IZ),SD_A_IN(1,1,IZ),MAXX_IN,MAXY_IN,NX,
     $   NY,BADDATA,UNCORRELATED,WT_X,MAX_REACH,REACH_X,WT_Y,REACH_Y,
     $   MAXX_OUT,MAXY_OUT,A_OUT(1,1,IZ),SD_A_OUT(1,1,IZ))
      ENDDO

      END SUBROUTINE FILTER_2D_DRIVER3
