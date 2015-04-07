      SUBROUTINE HOLE_FILL_2D_ONE_PASS(A_IN,SD_A_IN,CANDIDATE,MAXX_IN,
     $MAXY_IN,DELX,DELY,NX,NY,BADDATA,
     $MAX_CONSECUTIVE_UNOCCUPIED_OCTANTS,MAX_SEARCH_RADIUS,
     $MIN_VALUES_PER_OCTANT,MAX_VALUES_PER_OCTANT,MAXX_OUT,MAXY_OUT,
     $A_OUT,SD_A_OUT,POINTS_MISSING,POINTS_FILLED)

C  Thomas Matejka NOAA/NSSL 9 July 1997

C  This subroutine fills holes in a two-dimensional data field by
C  two-dimensional, first-order, weighted linear regression.  Holes in
C  the standard deviation of the data field are filled with the square
C  root of the weighted mean of the variance of the data used in the
C  regression.

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
C  is 2 or 3, set MAX_SEARCH_RADIUS to N.  If MAX_SEARCH_RADIUS is 0, no
C  filling is performed.

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
      INCLUDE 'tmmlib.inc'
      LOGICAL OCTANT_OCCUPIED(8)
      LOGICAL SUCCESS
      INTEGER MAXX_IN,MAXY_IN
      LOGICAL CANDIDATE(MAXX_IN,MAXY_IN)
      INTEGER MAXX_OUT,MAXY_OUT,MAX_CONSECUTIVE_UNOCCUPIED_OCTANTS,
     $POINTS_FILLED,POINTS_MISSING,IX,IY,TOTAL_POINTS,
     $CONSECUTIVE_UNOCCUPIED_OCTANTS,MAX_SEARCH_RADIUS,
     $MIN_VALUES_PER_OCTANT,MAX_VALUES_PER_OCTANT,POINTS_FOUND,IOCTANT,
     $IPOINT,EFFECTIVE_TOTAL_POINTS,NX,NY
      INTEGER IX_FOUND(NX*NY),IY_FOUND(NX*NY)
      REAL DELX,DELY,X,Y,SSQ,BADDATA,SUM,dum
      REAL COEFF(3),VAR_COEFF(3)
      REAL RESPONSE_VECTOR(NX*NY),VAR_RESPONSE_VECTOR(NX*NY)
      REAL REG_A_INV(3,3)
      REAL A_IN(MAXX_IN,MAXY_IN),SD_A_IN(MAXX_IN,MAXY_IN)
      REAL A_OUT(MAXX_OUT,MAXY_OUT),SD_A_OUT(MAXX_OUT,MAXY_OUT)
      REAL A_TEMP(NX,NY),SD_A_TEMP(NX,NY)
      REAL PREDICTOR_ARRAY(NX*NY,3)

C  Copy the input fields in case they are to be overwritten.
      DO 1 IY=1,NY
         DO 2 IX=1,NX
            A_TEMP(IX,IY)=A_IN(IX,IY)
            SD_A_TEMP(IX,IY)=SD_A_IN(IX,IY)
2        CONTINUE
1     CONTINUE

C  Initialize counters.
      POINTS_MISSING=0
      POINTS_FILLED=0

C  Loop through the points on the plane.
      DO 3 IY=1,NY
         DO 4 IX=1,NX

C  Check whether to fill the field value.
            IF(A_TEMP(IX,IY).NE.BADDATA.AND.
     $      SD_A_TEMP(IX,IY).NE.BADDATA)THEN
               A_OUT(IX,IY)=A_TEMP(IX,IY)
               SD_A_OUT(IX,IY)=SD_A_TEMP(IX,IY)
            ELSE
               POINTS_MISSING=POINTS_MISSING+1
               IF(CANDIDATE(IX,IY))THEN
                  IF(MAX_SEARCH_RADIUS.GE.1)THEN
                     CONSECUTIVE_UNOCCUPIED_OCTANTS=0
                     TOTAL_POINTS=0

C  Search octant 1.
                     CALL OCTANT1(A_TEMP,SD_A_TEMP,NX,NY,NX,NY,BADDATA,
     $               MAX_SEARCH_RADIUS,MAX_VALUES_PER_OCTANT,IX,IY,
     $               POINTS_FOUND,IX_FOUND(TOTAL_POINTS+1),
     $               IY_FOUND(TOTAL_POINTS+1))
                     IF(POINTS_FOUND.LT.MIN_VALUES_PER_OCTANT)THEN
                        CONSECUTIVE_UNOCCUPIED_OCTANTS=
     $                  CONSECUTIVE_UNOCCUPIED_OCTANTS+1
                        IF(CONSECUTIVE_UNOCCUPIED_OCTANTS.GT.
     $                  MAX_CONSECUTIVE_UNOCCUPIED_OCTANTS)THEN
                           A_OUT(IX,IY)=BADDATA
                           SD_A_OUT(IX,IY)=BADDATA
                           GOTO 5
                        ENDIF
                        OCTANT_OCCUPIED(1)=.FALSE.
                     ELSE
                        CONSECUTIVE_UNOCCUPIED_OCTANTS=0
                        OCTANT_OCCUPIED(1)=.TRUE.
                     ENDIF
                     TOTAL_POINTS=TOTAL_POINTS+POINTS_FOUND

C  Search octant 2.
                     CALL OCTANT2(A_TEMP,SD_A_TEMP,NX,NY,NX,NY,BADDATA,
     $               MAX_SEARCH_RADIUS,MAX_VALUES_PER_OCTANT,IX,IY,
     $               POINTS_FOUND,IX_FOUND(TOTAL_POINTS+1),
     $               IY_FOUND(TOTAL_POINTS+1))
                     IF(POINTS_FOUND.LT.MIN_VALUES_PER_OCTANT)THEN
                        CONSECUTIVE_UNOCCUPIED_OCTANTS=
     $                  CONSECUTIVE_UNOCCUPIED_OCTANTS+1
                        IF(CONSECUTIVE_UNOCCUPIED_OCTANTS.GT.
     $                  MAX_CONSECUTIVE_UNOCCUPIED_OCTANTS)THEN
                           A_OUT(IX,IY)=BADDATA
                           SD_A_OUT(IX,IY)=BADDATA
                           GOTO 5
                        ENDIF
                        OCTANT_OCCUPIED(2)=.FALSE.
                     ELSE
                        CONSECUTIVE_UNOCCUPIED_OCTANTS=0
                        OCTANT_OCCUPIED(2)=.TRUE.
                     ENDIF
                     TOTAL_POINTS=TOTAL_POINTS+POINTS_FOUND

C  Search octant 3.
                     CALL OCTANT3(A_TEMP,SD_A_TEMP,NX,NY,NX,NY,BADDATA,
     $               MAX_SEARCH_RADIUS,MAX_VALUES_PER_OCTANT,IX,IY,
     $               POINTS_FOUND,IX_FOUND(TOTAL_POINTS+1),
     $               IY_FOUND(TOTAL_POINTS+1))
                     IF(POINTS_FOUND.LT.MIN_VALUES_PER_OCTANT)THEN
                        CONSECUTIVE_UNOCCUPIED_OCTANTS=
     $                  CONSECUTIVE_UNOCCUPIED_OCTANTS+1
                        IF(CONSECUTIVE_UNOCCUPIED_OCTANTS.GT.
     $                  MAX_CONSECUTIVE_UNOCCUPIED_OCTANTS)THEN
                           A_OUT(IX,IY)=BADDATA
                           SD_A_OUT(IX,IY)=BADDATA
                           GOTO 5
                        ENDIF
                        OCTANT_OCCUPIED(3)=.FALSE.
                     ELSE
                        CONSECUTIVE_UNOCCUPIED_OCTANTS=0
                        OCTANT_OCCUPIED(3)=.TRUE.
                     ENDIF
                     TOTAL_POINTS=TOTAL_POINTS+POINTS_FOUND

C  Search octant 4.
                     CALL OCTANT4(A_TEMP,SD_A_TEMP,NX,NY,NX,NY,BADDATA,
     $               MAX_SEARCH_RADIUS,MAX_VALUES_PER_OCTANT,IX,IY,
     $               POINTS_FOUND,IX_FOUND(TOTAL_POINTS+1),
     $               IY_FOUND(TOTAL_POINTS+1))
                     IF(POINTS_FOUND.LT.MIN_VALUES_PER_OCTANT)THEN
                        CONSECUTIVE_UNOCCUPIED_OCTANTS=
     $                  CONSECUTIVE_UNOCCUPIED_OCTANTS+1
                        IF(CONSECUTIVE_UNOCCUPIED_OCTANTS.GT.
     $                  MAX_CONSECUTIVE_UNOCCUPIED_OCTANTS)THEN
                           A_OUT(IX,IY)=BADDATA
                           SD_A_OUT(IX,IY)=BADDATA
                           GOTO 5
                        ENDIF
                        OCTANT_OCCUPIED(4)=.FALSE.
                     ELSE
                        CONSECUTIVE_UNOCCUPIED_OCTANTS=0
                        OCTANT_OCCUPIED(4)=.TRUE.
                     ENDIF
                     TOTAL_POINTS=TOTAL_POINTS+POINTS_FOUND

C  Search octant 5.
                     CALL OCTANT5(A_TEMP,SD_A_TEMP,NX,NY,NX,NY,BADDATA,
     $               MAX_SEARCH_RADIUS,MAX_VALUES_PER_OCTANT,IX,IY,
     $               POINTS_FOUND,IX_FOUND(TOTAL_POINTS+1),
     $               IY_FOUND(TOTAL_POINTS+1))
                     IF(POINTS_FOUND.LT.MIN_VALUES_PER_OCTANT)THEN
                        CONSECUTIVE_UNOCCUPIED_OCTANTS=
     $                  CONSECUTIVE_UNOCCUPIED_OCTANTS+1
                        IF(CONSECUTIVE_UNOCCUPIED_OCTANTS.GT.
     $                  MAX_CONSECUTIVE_UNOCCUPIED_OCTANTS)THEN
                           A_OUT(IX,IY)=BADDATA
                           SD_A_OUT(IX,IY)=BADDATA
                           GOTO 5
                        ENDIF
                        OCTANT_OCCUPIED(5)=.FALSE.
                     ELSE
                        CONSECUTIVE_UNOCCUPIED_OCTANTS=0
                        OCTANT_OCCUPIED(5)=.TRUE.
                     ENDIF
                     TOTAL_POINTS=TOTAL_POINTS+POINTS_FOUND

C  Search octant 6.
                     CALL OCTANT6(A_TEMP,SD_A_TEMP,NX,NY,NX,NY,BADDATA,
     $               MAX_SEARCH_RADIUS,MAX_VALUES_PER_OCTANT,IX,IY,
     $               POINTS_FOUND,IX_FOUND(TOTAL_POINTS+1),
     $               IY_FOUND(TOTAL_POINTS+1))
                     IF(POINTS_FOUND.LT.MIN_VALUES_PER_OCTANT)THEN
                        CONSECUTIVE_UNOCCUPIED_OCTANTS=
     $                  CONSECUTIVE_UNOCCUPIED_OCTANTS+1
                        IF(CONSECUTIVE_UNOCCUPIED_OCTANTS.GT.
     $                  MAX_CONSECUTIVE_UNOCCUPIED_OCTANTS)THEN
                           A_OUT(IX,IY)=BADDATA
                           SD_A_OUT(IX,IY)=BADDATA
                           GOTO 5
                        ENDIF
                        OCTANT_OCCUPIED(6)=.FALSE.
                     ELSE
                        CONSECUTIVE_UNOCCUPIED_OCTANTS=0
                        OCTANT_OCCUPIED(6)=.TRUE.
                     ENDIF
                     TOTAL_POINTS=TOTAL_POINTS+POINTS_FOUND

C  Search octant 7.
                     CALL OCTANT7(A_TEMP,SD_A_TEMP,NX,NY,NX,NY,BADDATA,
     $               MAX_SEARCH_RADIUS,MAX_VALUES_PER_OCTANT,IX,IY,
     $               POINTS_FOUND,IX_FOUND(TOTAL_POINTS+1),
     $               IY_FOUND(TOTAL_POINTS+1))
                     IF(POINTS_FOUND.LT.MIN_VALUES_PER_OCTANT)THEN
                        CONSECUTIVE_UNOCCUPIED_OCTANTS=
     $                  CONSECUTIVE_UNOCCUPIED_OCTANTS+1
                        IF(CONSECUTIVE_UNOCCUPIED_OCTANTS.GT.
     $                  MAX_CONSECUTIVE_UNOCCUPIED_OCTANTS)THEN
                           A_OUT(IX,IY)=BADDATA
                           SD_A_OUT(IX,IY)=BADDATA
                           GOTO 5
                        ENDIF
                        OCTANT_OCCUPIED(7)=.FALSE.
                     ELSE
                        CONSECUTIVE_UNOCCUPIED_OCTANTS=0
                        OCTANT_OCCUPIED(7)=.TRUE.
                     ENDIF
                     TOTAL_POINTS=TOTAL_POINTS+POINTS_FOUND

C  Search octant 8.
                     CALL OCTANT8(A_TEMP,SD_A_TEMP,NX,NY,NX,NY,BADDATA,
     $               MAX_SEARCH_RADIUS,MAX_VALUES_PER_OCTANT,IX,IY,
     $               POINTS_FOUND,IX_FOUND(TOTAL_POINTS+1),
     $               IY_FOUND(TOTAL_POINTS+1))
                     IF(POINTS_FOUND.LT.MIN_VALUES_PER_OCTANT)THEN
                        CONSECUTIVE_UNOCCUPIED_OCTANTS=
     $                  CONSECUTIVE_UNOCCUPIED_OCTANTS+1
                        IF(CONSECUTIVE_UNOCCUPIED_OCTANTS.GT.
     $                  MAX_CONSECUTIVE_UNOCCUPIED_OCTANTS)THEN
                           A_OUT(IX,IY)=BADDATA
                           SD_A_OUT(IX,IY)=BADDATA
                           GOTO 5
                        ENDIF
                        OCTANT_OCCUPIED(8)=.FALSE.
                     ELSE
                        CONSECUTIVE_UNOCCUPIED_OCTANTS=0
                        OCTANT_OCCUPIED(8)=.TRUE.
                     ENDIF
                     TOTAL_POINTS=TOTAL_POINTS+POINTS_FOUND

C  Check whether there are too many consecutive unoccupied octants
C  including octants 8 and 1.
                     IF(.NOT.OCTANT_OCCUPIED(8))THEN
                        DO 6 IOCTANT=1,6
                           IF(.NOT.OCTANT_OCCUPIED(IOCTANT))THEN
                              CONSECUTIVE_UNOCCUPIED_OCTANTS=
     $                        CONSECUTIVE_UNOCCUPIED_OCTANTS+1
                              IF(CONSECUTIVE_UNOCCUPIED_OCTANTS.GT.
     $                        MAX_CONSECUTIVE_UNOCCUPIED_OCTANTS)THEN
                                 A_OUT(IX,IY)=BADDATA
                                 SD_A_OUT(IX,IY)=BADDATA
                                 GOTO 5
                              ENDIF
                           ELSE
                              GOTO 7
                           ENDIF
6                       CONTINUE
                     ENDIF
7                    CONTINUE

C  Load the data into arrays for the regression.
                     DO 8 IPOINT=1,TOTAL_POINTS
                        X=FLOAT(IX_FOUND(IPOINT)-IX)*DELX
                        Y=FLOAT(IY_FOUND(IPOINT)-IY)*DELY
                        PREDICTOR_ARRAY(IPOINT,1)=1.
                        PREDICTOR_ARRAY(IPOINT,2)=X
                        PREDICTOR_ARRAY(IPOINT,3)=Y
                        RESPONSE_VECTOR(IPOINT)=
     $                  A_TEMP(IX_FOUND(IPOINT),IY_FOUND(IPOINT))
                        VAR_RESPONSE_VECTOR(IPOINT)=
     $                  SD_A_TEMP(IX_FOUND(IPOINT),IY_FOUND(IPOINT))**2
8                    CONTINUE

C  Perform the regression.
                     CALL LLS_VAR(3,TOTAL_POINTS,PREDICTOR_ARRAY,
     $               NX*NY,RESPONSE_VECTOR,
     $               VAR_RESPONSE_VECTOR,.FALSE.,dum,.FALSE.,1.,.TRUE.,
     $               HOLE_FILL_SINGULAR_THRESHOLD,0,3,
     $               EFFECTIVE_TOTAL_POINTS,COEFF,VAR_COEFF,SSQ,
     $               REG_A_INV,SUCCESS)
                     IF(SUCCESS)THEN

C  Evaluate the regression at the point to be filled, and estimate the
C  standard deviation of the filled point as the square root of the
C  weighted mean of the variance of the data used in the regression.
                        A_OUT(IX,IY)=COEFF(1)
                        SUM=0.
                        DO 9 IPOINT=1,TOTAL_POINTS
                           SUM=SUM+1./VAR_RESPONSE_VECTOR(IPOINT)
9                       CONTINUE
                        SD_A_OUT(IX,IY)=SQRT(FLOAT(TOTAL_POINTS)/SUM)
                        POINTS_FILLED=POINTS_FILLED+1
                     ELSE
                        A_OUT(IX,IY)=BADDATA
                        SD_A_OUT(IX,IY)=BADDATA
                     ENDIF
                  ELSE
                     A_OUT(IX,IY)=BADDATA
                     SD_A_OUT(IX,IY)=BADDATA
                  ENDIF
               ELSE
                  A_OUT(IX,IY)=BADDATA
                  SD_A_OUT(IX,IY)=BADDATA
               ENDIF
            ENDIF
5           CONTINUE
4        CONTINUE
3     CONTINUE

C  Done.
      RETURN
      END
