      SUBROUTINE HOLE_FILL_1D_ONE_PASS(A_IN,SD_A_IN,CANDIDATE,NX,
     $BADDATA,MAX_UNOCCUPIED_SIDES,MAX_SEARCH_RADIUS,
     $MIN_VALUES_PER_SIDE,MAX_VALUES_PER_SIDE,A_OUT,SD_A_OUT,
     $POINTS_MISSING,POINTS_FILLED)

C  Thomas Matejka NOAA/NSSL 13 September 1995

C  This subroutine fills holes in a one-dimensional data field by
C  one-dimensional, first-order, weighted linear regression.  Holes in
C  the standard deviation of the data field are filled with the square
C  root of the weighted mean of the variance of the data used in the
C  regression.

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
C  set MAX_SEARCH_RADIUS to N.  If MAX_SEARCH_RADIUS is 0, no filling is
C  performed.

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
      INCLUDE 'tmmlib.inc'
      LOGICAL SIDE_OCCUPIED(2)
      LOGICAL SUCCESS
      INTEGER NX
      LOGICAL CANDIDATE(NX)
      INTEGER MAX_UNOCCUPIED_SIDES,POINTS_FILLED,POINTS_MISSING,
     $IX,TOTAL_POINTS,UNOCCUPIED_SIDES,MAX_SEARCH_RADIUS,
     $MIN_VALUES_PER_SIDE,MAX_VALUES_PER_SIDE,POINTS_FOUND,IPOINT,
     $EFFECTIVE_TOTAL_POINTS
      INTEGER IX_FOUND(NX)
      REAL X,SSQ,BADDATA,SUM,dum
      REAL COEFF(2),VAR_COEFF(2)
      REAL RESPONSE_VECTOR(NX),VAR_RESPONSE_VECTOR(NX),A_TEMP(NX),
     $SD_A_TEMP(NX)
      REAL REG_A_INV(2,2)
      REAL A_IN(NX),SD_A_IN(NX),A_OUT(NX),SD_A_OUT(NX)
      REAL PREDICTOR_ARRAY(NX,2)

C  Copy the input fields in case they are to be overwritten.
      DO 1 IX=1,NX
         A_TEMP(IX)=A_IN(IX)
         SD_A_TEMP(IX)=SD_A_IN(IX)
1     CONTINUE

C  Initialize counters.
      POINTS_MISSING=0
      POINTS_FILLED=0

C  Loop through the points.
      DO 2 IX=1,NX

C  Check whether to fill the field value.
         IF(A_TEMP(IX).NE.BADDATA.AND.
     $   SD_A_TEMP(IX).NE.BADDATA)THEN
            A_OUT(IX)=A_TEMP(IX)
            SD_A_OUT(IX)=SD_A_TEMP(IX)
         ELSE
            POINTS_MISSING=POINTS_MISSING+1
            IF(CANDIDATE(IX))THEN
               IF(MAX_SEARCH_RADIUS.GE.1)THEN
                  UNOCCUPIED_SIDES=0
                  TOTAL_POINTS=0

C  Search side 1.
                  CALL SIDE1(A_TEMP,SD_A_TEMP,NX,BADDATA,
     $            MAX_SEARCH_RADIUS,MAX_VALUES_PER_SIDE,IX,POINTS_FOUND,
     $            IX_FOUND(TOTAL_POINTS+1))
                  IF(POINTS_FOUND.LT.MIN_VALUES_PER_SIDE)THEN
                     UNOCCUPIED_SIDES=UNOCCUPIED_SIDES+1
                     IF(UNOCCUPIED_SIDES.GT.MAX_UNOCCUPIED_SIDES)THEN
                        A_OUT(IX)=BADDATA
                        SD_A_OUT(IX)=BADDATA
                        GOTO 3
                     ENDIF
                     SIDE_OCCUPIED(1)=.FALSE.
                  ELSE
                     SIDE_OCCUPIED(1)=.TRUE.
                  ENDIF
                  TOTAL_POINTS=TOTAL_POINTS+POINTS_FOUND

C  Search side 2.
                  CALL SIDE2(A_TEMP,SD_A_TEMP,NX,BADDATA,
     $            MAX_SEARCH_RADIUS,MAX_VALUES_PER_SIDE,IX,POINTS_FOUND,
     $            IX_FOUND(TOTAL_POINTS+1))
                  IF(POINTS_FOUND.LT.MIN_VALUES_PER_SIDE)THEN
                     UNOCCUPIED_SIDES=UNOCCUPIED_SIDES+1
                     IF(UNOCCUPIED_SIDES.GT.
     $               MAX_UNOCCUPIED_SIDES)THEN
                        A_OUT(IX)=BADDATA
                        SD_A_OUT(IX)=BADDATA
                        GOTO 3
                     ENDIF
                     SIDE_OCCUPIED(2)=.FALSE.
                  ELSE
                     SIDE_OCCUPIED(2)=.TRUE.
                  ENDIF
                  TOTAL_POINTS=TOTAL_POINTS+POINTS_FOUND

C  Check whether there are too many unoccupied sides.
                  IF(UNOCCUPIED_SIDES.GT.MAX_UNOCCUPIED_SIDES)THEN
                     A_OUT(IX)=BADDATA
                     SD_A_OUT(IX)=BADDATA
                     GOTO 3
                  ENDIF

C  Load the data into arrays for the regression.
                  DO 4 IPOINT=1,TOTAL_POINTS
                     X=FLOAT(IX_FOUND(IPOINT)-IX)
                     PREDICTOR_ARRAY(IPOINT,1)=1.
                     PREDICTOR_ARRAY(IPOINT,2)=X
                     RESPONSE_VECTOR(IPOINT)=A_TEMP(IX_FOUND(IPOINT))
                     VAR_RESPONSE_VECTOR(IPOINT)=
     $               SD_A_TEMP(IX_FOUND(IPOINT))**2
4                 CONTINUE

C  Perform the regression.
                  CALL LLS_VAR(2,TOTAL_POINTS,PREDICTOR_ARRAY,NX,
     $            RESPONSE_VECTOR,VAR_RESPONSE_VECTOR,.FALSE.,dum,.FALSE.,
     $            1.,.TRUE.,HOLE_FILL_SINGULAR_THRESHOLD,0,2,
     $            EFFECTIVE_TOTAL_POINTS,COEFF,VAR_COEFF,SSQ,REG_A_INV,
     $            SUCCESS)
                  IF(SUCCESS)THEN

C  Evaluate the regression at the point to be filled, and estimate the
C  standard deviation of the filled point as the square root of the
C  weighted mean of the variance of the data used in the regression.
                     A_OUT(IX)=COEFF(1)
                     SUM=0.
                     DO 5 IPOINT=1,TOTAL_POINTS
                        SUM=SUM+1./VAR_RESPONSE_VECTOR(IPOINT)
5                    CONTINUE
                     SD_A_OUT(IX)=SQRT(FLOAT(TOTAL_POINTS)/SUM)
                     POINTS_FILLED=POINTS_FILLED+1
                  ELSE
                     A_OUT(IX)=BADDATA
                     SD_A_OUT(IX)=BADDATA
                  ENDIF
               ELSE
                  A_OUT(IX)=BADDATA
                  SD_A_OUT(IX)=BADDATA
               ENDIF
            ELSE
               A_OUT(IX)=BADDATA
               SD_A_OUT(IX)=BADDATA
            ENDIF
         ENDIF
3        CONTINUE
2     CONTINUE

C  Done.
      RETURN
      END
