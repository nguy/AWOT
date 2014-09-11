      SUBROUTINE HOLE_FILL_1D_DRIVER12(A_IN,SD_A_IN,CANDIDATE,MAXX_IN,
     $MAXY_IN,MAXZ_IN,NX,NY,NZ,BADDATA,MAX_UNOCCUPIED_SIDES,
     $MAX_SEARCH_RADIUS,MIN_VALUES_PER_SIDE,MAX_VALUES_PER_SIDE,
     $MAXX_OUT,MAXY_OUT,MAXZ_OUT,A_OUT,SD_A_OUT,POINTS_MISSING,
     $POINTS_FILLED)

C  Thomas Matejka NOAA/NSSL 23 April 1997

C  This subroutine is a first and second-dimensional driver for
C  subroutine HOLE_FILL_1D, which is applied in the third dimension.

C  Input:

C  A_IN is a three-dimensional real array.  A_IN(I,J,K) specifies the
C  value, at the Ith grid point in the first dimension, the Jth grid
C  point in the second dimension, and the Kth grid point in the third
C  dimension, of the data field whose holes are to be filled.  If it is
C  missing, it should equal BADDATA.

C  SD_A_IN is a three-dimensional real array.  SD_A_IN(I,J,K) specifies
C  the standard deviation of A_IN(I,J,K).  If it is missing, it should
C  equal BADDATA.  SD_A_IN(I,J,K) must not be missing for A_IN(I,J,K) to
C  be usable.

C  CANDIDATE is a three-dimensional logical array.  If A_IN(I,J,K) is
C  missing, an attempt to fill it will be made if and only if
C  CANDIDATE(I,J,K) is .TRUE..

C  MAXX_IN is an integer variable that specifies the first dimension of
C  A_IN, SD_A_IN, and CANDIDATE in the calling program.

C  MAXY_IN is an integer variable that specifies the second dimension of
C  A_IN, SD_A_IN, and CANDIDATE in the calling program.

C  MAXZ_IN is an integer variable that specifies the third dimension of
C  A_IN, SD_A_IN, and CANDIDATE in the calling program.

C  NX is an integer variable that specifies the number of grid points in
C  the first dimension for A_IN, SD_A_IN, CANDIDATE, A_OUT, and
C  SD_A_OUT.

C  NY is an integer variable that specifies the number of grid points in
C  the second dimension for A_IN, SD_A_IN, CANDIDATE, A_OUT, and
C  SD_A_OUT.

C  NZ is an integer variable that specifies the number of grid points in
C  the third dimension for A_IN, SD_A_IN, CANDIDATE, A_OUT, and
C  SD_A_OUT.

C  BADDATA is a real variable that indicates a missing value as
C  described.

C  MAX_UNOCCUPIED_SIDES is an integer variable that specifies the
C  maximum number of sides without sufficient, close data surrounding a
C  point for filling of that point to occur.  MAX_UNOCCUPIED_SIDES must
C  be 0 or 1.  A value of 0 ensures that no extrapolation will ever
C  occur.

C  MAX_SEARCH_RADIUS is an integer variable that specifies how many
C  points away from the point to be filled the search for surrounding
C  data is conducted.  MAX_SEARCH_RADIUS must be greater than or equal
C  to 0.  MAX_SEARCH_RADIUS controls the size of holes to be filled.  To
C  fill holes up to N grid points long when MAX_UNOCCUPIED_SIDES is 0,
C  set MAX_SEARCH_RADIUS to N/2 + 1.  If MAX_SEARCH_RADIUS is 0, no
C  filling is performed.

C  MIN_VALUES_PER_SIDE is an integer variable that specifies the minimum
C  number of data points per side that are required for the side to be
C  considered occupied.  If fewer than MIN_VALUES_PER_SIDE points are
C  found on an side, they are still used to fill the hole if enough
C  occupied sides are found.  MIN_VALUES_PER_SIDE must be greater than
C  or equal to 1.

C  MAX_VALUES_PER_SIDE is an integer variable that specifies the maximum
C  number of surrounding data values allowed in each side.  Once
C  MAX_VALUES_PER_SIDE data values are found on either side, searches in
C  that side are stopped even if the search has not yet reached the
C  maximum search radius.  MAX_VALUES_PER_SIDE must be greater than or
C  equal to MIN_VALUES_PER_SIDE.

C  MAXX_OUT is an integer variable that specifies the first dimension of
C  A_OUT and SD_A_OUT in the calling program.

C  MAXY_OUT is an integer variable that specifies the second dimension
C  of A_OUT and SD_A_OUT in the calling program.

C  MAXZ_OUT is an integer variable that specifies the third dimension of
C  A_OUT and SD_A_OUT in the calling program.

C  Output:

C  A_OUT is a three-dimensional real array.  A_OUT(I,J,K) returns the
C  value, at the Ith grid point in the first dimension, the Jth grid
C  point in the second dimension, and the Kth grid point in the third
C  dimension, of the data field whose holes have been filled.  If it is
C  still missing, it is returned as BADDATA.  A_OUT may be identical to
C  A_IN, in which case A_IN will be overwritten.

C  SD_A_OUT is a three-dimensional real array.  SD_A_OUT(I,J,K) returns
C  the standard deviation of A_OUT(I,J,K).  If it is missing, it is
C  returned as BADDATA.  SD_A_OUT may be identical to SD_A_IN, in which
C  case SD_A_IN will be overwritten.

C  POINTS_MISSING is an integer variable that returns the original
C  number of missing data points.

C  POINTS_FILLED is an integer variable that returns the number of
C  points that were filled.

      IMPLICIT NONE
      INTEGER MAXX_IN,MAXY_IN,MAXZ_IN
      LOGICAL CANDIDATE(MAXX_IN,MAXY_IN,MAXZ_IN)
      INTEGER MAXX_OUT,MAXY_OUT,MAXZ_OUT,MAX_UNOCCUPIED_SIDES,
     $POINTS_FILLED,POINTS_MISSING,IX,IY,IZ,MAX_SEARCH_RADIUS,
     $MIN_VALUES_PER_SIDE,MAX_VALUES_PER_SIDE,NX,NY,NZ,
     $POINTS_MISSING_COLUMN,POINTS_FILLED_COLUMN
      REAL BADDATA
      LOGICAL CANDIDATE_1D(NZ)
      REAL A_1D(NZ),SD_A_1D(NZ)
      REAL A_IN(MAXX_IN,MAXY_IN,MAXZ_IN),
     $SD_A_IN(MAXX_IN,MAXY_IN,MAXZ_IN)
      REAL A_OUT(MAXX_OUT,MAXY_OUT,MAXZ_OUT),
     $SD_A_OUT(MAXX_OUT,MAXY_OUT,MAXZ_OUT)

C  Loop through the columns and fill data holes.
      POINTS_MISSING=0
      POINTS_FILLED=0
      DO 1 IY=1,NY
         DO 2 IX=1,NX
            DO 3 IZ=1,NZ
               A_1D(IZ)=A_IN(IX,IY,IZ)
               SD_A_1D(IZ)=SD_A_IN(IX,IY,IZ)
               CANDIDATE_1D(IZ)=CANDIDATE(IX,IY,IZ)
3           CONTINUE
            CALL HOLE_FILL_1D(A_1D,SD_A_1D,CANDIDATE_1D,NZ,BADDATA,
     $      MAX_UNOCCUPIED_SIDES,MAX_SEARCH_RADIUS,MIN_VALUES_PER_SIDE,
     $      MAX_VALUES_PER_SIDE,A_1D,SD_A_1D,POINTS_MISSING_COLUMN,
     $      POINTS_FILLED_COLUMN)
            DO 4 IZ=1,NZ
               A_OUT(IX,IY,IZ)=A_1D(IZ)
               SD_A_OUT(IX,IY,IZ)=SD_A_1D(IZ)
4           CONTINUE
            POINTS_MISSING=POINTS_MISSING+POINTS_MISSING_COLUMN
            POINTS_FILLED=POINTS_FILLED+POINTS_FILLED_COLUMN
2        CONTINUE
1     CONTINUE

C  Done.
      RETURN
      END
