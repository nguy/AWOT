      SUBROUTINE HOLE_FILL_EXTEND_1D(A_IN,SD_A_IN,CANDIDATE,NX,BADDATA,
     $MAX_UNOCCUPIED_SIDES,MAX_SEARCH_RADIUS,MIN_VALUES_PER_SIDE,
     $MAX_VALUES_PER_SIDE,A_OUT,SD_A_OUT,POINTS_MISSING,POINTS_FILLED)

C  Thomas Matejka NOAA/NSSL 16 February 1995

C  This subroutine fills holes in a one-dimensional data field by
C  one-dimensional, zeroth-order, weighted linear regression.  Holes in
C  the standard deviation of the data field are filled with the square
C  root of the weighted mean of the variance of the data used in the
C  regression.  This subroutine performs two passes of hole filling.
C  Performing two passes eliminates artifacts such as the filling of the
C  centers of holes too large to be completely filled.  Note that
C  MAX_SEARCH_RADIUS should be smaller than that used in one pass hole
C  filling to get similar behavior.

C  Input:

C  A_IN is a one-dimensional real array.  A_IN(I) specifies the value,
C  at the Ith grid point, of the data field whose holes are to be
C  filled.  If it is missing, it should equal BADDATA.

C  SD_A_IN is a one-dimensional real array.  SD_A_IN(I) specifies the
C  standard deviation of A_IN(I).  If it is missing, it should equal
C  BADDATA.  SD_A_IN(I) must not be missing for A_IN(I) to be usable.

C  CANDIDATE is a one-dimensional logical array.  If A_IN(I) is missing,
C  an attempt to fill it will be made if and only if CANDIDATE(I) is
C  .TRUE..

C  NX is an integer variable that specifies the number of grid points
C  for A_IN, SD_A_IN, CANDIDATE, A_OUT, and SD_A_OUT.

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

C  Output:

C  A_OUT is a one-dimensional real array.  A_OUT(I) returns the value,
C  at the Ith grid point, of the data field whose holes have been
C  filled.  If it is still missing, it is returned as BADDATA.  A_OUT
C  may be identical to A_IN, in which case A_IN will be overwritten.

C  SD_A_OUT is a one-dimensional real array.  SD_A_OUT(I) returns the
C  standard deviation of A_OUT(I).  If it is missing, it is returned as
C  BADDATA.  SD_A_OUT may be identical to SD_A_IN, in which case SD_A_IN
C  will be overwritten.

C  POINTS_MISSING is an integer variable that returns the original
C  number of missing data points.

C  POINTS_FILLED is an integer variable that returns the number of
C  points that were filled.

      IMPLICIT NONE
      INTEGER NX
      LOGICAL CANDIDATE(NX)
      INTEGER MAX_UNOCCUPIED_SIDES,POINTS_FILLED,POINTS_MISSING,
     $MAX_SEARCH_RADIUS,MIN_VALUES_PER_SIDE,MAX_VALUES_PER_SIDE,
     $POINTS_FILLED_PASS,POINTS_MISSING_PASS
      REAL BADDATA
      REAL A_IN(NX),SD_A_IN(NX),A_OUT(NX),SD_A_OUT(NX)

C  Perform two passes of hole filling.
      CALL HOLE_FILL_EXTEND_1D_ONE_PASS(A_IN,SD_A_IN,CANDIDATE,NX,
     $BADDATA,MAX_UNOCCUPIED_SIDES,MAX_SEARCH_RADIUS,
     $MIN_VALUES_PER_SIDE,MAX_VALUES_PER_SIDE,A_OUT,SD_A_OUT,
     $POINTS_MISSING_PASS,POINTS_FILLED_PASS)
      POINTS_MISSING=POINTS_MISSING_PASS
      POINTS_FILLED=POINTS_FILLED_PASS
      CALL HOLE_FILL_EXTEND_1D_ONE_PASS(A_OUT,SD_A_OUT,CANDIDATE,NX,
     $BADDATA,MAX_UNOCCUPIED_SIDES,MAX_SEARCH_RADIUS,
     $MIN_VALUES_PER_SIDE,MAX_VALUES_PER_SIDE,A_OUT,SD_A_OUT,
     $POINTS_MISSING_PASS,POINTS_FILLED_PASS)
      POINTS_FILLED=POINTS_FILLED+POINTS_FILLED_PASS

C  Done.
      RETURN
      END
