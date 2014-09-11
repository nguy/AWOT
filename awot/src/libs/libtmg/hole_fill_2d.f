      SUBROUTINE HOLE_FILL_2D(A_IN,SD_A_IN,CANDIDATE,MAXX_IN,MAXY_IN,
     $DELX,DELY,NX,NY,BADDATA,MAX_CONSECUTIVE_UNOCCUPIED_OCTANTS,
     $MAX_SEARCH_RADIUS,MIN_VALUES_PER_OCTANT,MAX_VALUES_PER_OCTANT,
     $MAXX_OUT,MAXY_OUT,A_OUT,SD_A_OUT,POINTS_MISSING,POINTS_FILLED)

C  Thomas Matejka NOAA/NSSL 18 August 1994

C  This subroutine fills holes in a two-dimensional data field by
C  two-dimensional, first-order, weighted linear regression.  Holes in
C  the standard deviation of the data field are filled with the square
C  root of the weighted mean of the variance of the data used in the
C  regression.  This subroutine performs two passes of hole filling.
C  Performing two passes eliminates artifacts such as the filling of the
C  centers of holes too large to be completely filled.  Note that
C  MAX_SEARCH_RADIUS should be smaller than that used in one pass hole
C  filling to get similar behavior.

C  Input:

C  A_IN is a two-dimensional real array.  A_IN(I,J) specifies the value,
C  at the Ith grid point in the first dimension and the Jth grid point
C  in the second dimension, of the data field whose holes are to be
C  filled.  If it is missing, it should equal BADDATA.

C  SD_A_IN is a two-dimensional real array.  SD_A_IN(I,J) specifies the
C  standard deviation of A_IN(I,J).  If it is missing, it should equal
C  BADDATA.  SD_A_IN(I,J) must not be missing for A_IN(I,J) to be
C  usable.

C  CANDIDATE is a two-dimensional logical array.  If A_IN(I,J) is
C  missing, an attempt to fill it will be made if and only if
C  CANDIDATE(I,J) is .TRUE..

C  MAXX_IN is an integer variable that specifies the first dimension of
C  A_IN, SD_A_IN, and CANDIDATE in the calling program.

C  MAXY_IN is an integer variable that specifies the second dimension of
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

C  Output:

C  A_OUT is a two-dimensional real array.  A_OUT(I,J) returns the value,
C  at the Ith grid point in the first dimension and the Jth grid point
C  in the second dimension, of the data field whose holes have been
C  filled.  If it is still missing, it is returned as BADDATA.  A_OUT
C  may be identical to A_IN, in which case A_IN will be overwritten.

C  SD_A_OUT is a two-dimensional real array.  SD_A_OUT(I,J) returns the
C  standard deviation of A_OUT(I,J).  If it is missing, it is returned
C  as BADDATA.  SD_A_OUT may be identical to SD_A_IN, in which case
C  SD_A_IN will be overwritten.

C  POINTS_MISSING is an integer variable that returns the original
C  number of missing data points.

C  POINTS_FILLED is an integer variable that returns the number of
C  points that were filled.

      IMPLICIT NONE
      INTEGER MAXX_IN,MAXY_IN
      LOGICAL CANDIDATE(MAXX_IN,MAXY_IN)
      INTEGER MAXX_OUT,MAXY_OUT,MAX_CONSECUTIVE_UNOCCUPIED_OCTANTS,
     $POINTS_FILLED,POINTS_MISSING,MAX_SEARCH_RADIUS,
     $MIN_VALUES_PER_OCTANT,MAX_VALUES_PER_OCTANT,NX,NY,
     $POINTS_FILLED_PASS,POINTS_MISSING_PASS
      REAL DELX,DELY,BADDATA
      REAL A_IN(MAXX_IN,MAXY_IN),SD_A_IN(MAXX_IN,MAXY_IN)
      REAL A_OUT(MAXX_OUT,MAXY_OUT),SD_A_OUT(MAXX_OUT,MAXY_OUT)

C  Perform two passes of hole filling.
      CALL HOLE_FILL_2D_ONE_PASS(A_IN,SD_A_IN,CANDIDATE,MAXX_IN,
     $MAXY_IN,DELX,DELY,NX,NY,BADDATA,
     $MAX_CONSECUTIVE_UNOCCUPIED_OCTANTS,MAX_SEARCH_RADIUS,
     $MIN_VALUES_PER_OCTANT,MAX_VALUES_PER_OCTANT,MAXX_OUT,MAXY_OUT,
     $A_OUT,SD_A_OUT,POINTS_MISSING_PASS,POINTS_FILLED_PASS)
      POINTS_MISSING=POINTS_MISSING_PASS
      POINTS_FILLED=POINTS_FILLED_PASS
      CALL HOLE_FILL_2D_ONE_PASS(A_OUT,SD_A_OUT,CANDIDATE,MAXX_OUT,
     $MAXY_OUT,DELX,DELY,NX,NY,BADDATA,
     $MAX_CONSECUTIVE_UNOCCUPIED_OCTANTS,MAX_SEARCH_RADIUS,
     $MIN_VALUES_PER_OCTANT,MAX_VALUES_PER_OCTANT,MAXX_OUT,MAXY_OUT,
     $A_OUT,SD_A_OUT,POINTS_MISSING_PASS,POINTS_FILLED_PASS)
      POINTS_FILLED=POINTS_FILLED+POINTS_FILLED_PASS

C  Done.
      RETURN
      END
