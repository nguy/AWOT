      SUBROUTINE HOLE_FILL_EXTEND_2D_DRIVER3(A_IN,SD_A_IN,CANDIDATE,
     $MAXX_IN,MAXY_IN,MAXZ_IN,DELX,DELY,NX,NY,NZ,BADDATA,
     $MAX_CONSECUTIVE_UNOCCUPIED_OCTANTS,MAX_SEARCH_RADIUS,
     $MIN_VALUES_PER_OCTANT,MAX_VALUES_PER_OCTANT,MAXX_OUT,MAXY_OUT,
     $MAXZ_OUT,A_OUT,SD_A_OUT,POINTS_MISSING,POINTS_FILLED)

C  Thomas Matejka NOAA/NSSL 16 February 1995

C  This subroutine is a third-dimensional driver for subroutine
C  HOLE_FILL_EXTEND_2D, which is applied in the first two dimensions.

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

C  DELX is a real variable that specifies the grid point spacing in the
C  first dimension for A_IN, SD_A_IN, CANDIDATE, A_OUT, and SD_A_OUT.

C  DELY is a real variable that specifies the grid point spacing in the
C  second dimension for A_IN, SD_A_IN, CANDIDATE, A_OUT, and SD_A_OUT.

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

C  MAX_CONSECUTIVE_UNOCCUPIED_OCTANTS is an integer variable that
C  specifies the maximum number of consecutive octants without
C  sufficient, close data surrounding a point for filling of that point
C  to occur.  MAX_CONSECUTIVE_UNOCCUPIED_OCTANTS must be from 0 to 7.  A
C  value of 2 ensures that no extrapolation will ever occur.  A value of
C  3 usually prevents extrapolation from occurring.

C  MAX_SEARCH_RADIUS is an integer variable that specifies how many
C  points away from the point to be filled the search for surrounding
C  data is conducted.  MAX_SEARCH_RADIUS must be greater than or equal
C  to 0.  MAX_SEARCH_RADIUS controls the size of holes to be filled.  To
C  fill square holes up to N grid points in diameter and long channels
C  up to N grid points in width when MAX_CONSECUTIVE_UNOCCUPIED_OCTANTS
C  is 2 or 3, set MAX_SEARCH_RADIUS to N/2 + 1.  If MAX_SEARCH_RADIUS is
C  0, no filling is performed.

C  MIN_VALUES_PER_OCTANT is an integer variable that specifies the
C  minimum number of data points per octant that are required for the
C  octant to be considered occupied.  If fewer than
C  MIN_VALUES_PER_OCTANT points are found in an octant, they are still
C  used to fill the hole if enough occupied octants are found.
C  MIN_VALUES_PER_OCTANT must be greater than or equal to 1.

C  MAX_VALUES_PER_OCTANT is an integer variable that specifies the
C  maximum number of surrounding data values allowed in each octant.
C  Once MAX_VALUES_PER_OCTANT data values are found in any octant,
C  searches in that octant are stopped even if the search has not yet
C  reached the maximum search radius.  MAX_VALUES_PER_OCTANT must be
C  greater than or equal to MIN_VALUES_PER_OCTANT.

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
      INTEGER MAXX_OUT,MAXY_OUT,MAXZ_OUT,
     $MAX_CONSECUTIVE_UNOCCUPIED_OCTANTS,POINTS_FILLED,POINTS_MISSING,
     $IZ,MAX_SEARCH_RADIUS,MIN_VALUES_PER_OCTANT,MAX_VALUES_PER_OCTANT,
     $NX,NY,NZ,POINTS_MISSING_PLANE,POINTS_FILLED_PLANE
      REAL DELX,DELY,BADDATA
      REAL A_IN(MAXX_IN,MAXY_IN,MAXZ_IN),
     $SD_A_IN(MAXX_IN,MAXY_IN,MAXZ_IN)
      REAL A_OUT(MAXX_OUT,MAXY_OUT,MAXZ_OUT),
     $SD_A_OUT(MAXX_OUT,MAXY_OUT,MAXZ_OUT)

C  Loop through the planes and fill data holes.
      POINTS_MISSING=0
      POINTS_FILLED=0
      DO 1 IZ=1,NZ
         CALL HOLE_FILL_EXTEND_2D(A_IN(1,1,IZ),SD_A_IN(1,1,IZ),
     $   CANDIDATE(1,1,IZ),MAXX_IN,MAXY_IN,DELX,DELY,NX,NY,BADDATA,
     $   MAX_CONSECUTIVE_UNOCCUPIED_OCTANTS,MAX_SEARCH_RADIUS,
     $   MIN_VALUES_PER_OCTANT,MAX_VALUES_PER_OCTANT,MAXX_OUT,MAXY_OUT,
     $   A_OUT(1,1,IZ),SD_A_OUT(1,1,IZ),POINTS_MISSING_PLANE,
     $   POINTS_FILLED_PLANE)
         POINTS_MISSING=POINTS_MISSING+POINTS_MISSING_PLANE
         POINTS_FILLED=POINTS_FILLED+POINTS_FILLED_PLANE
1     CONTINUE

C  Done.
      RETURN
      END
