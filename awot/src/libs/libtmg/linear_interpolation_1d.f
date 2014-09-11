      REAL FUNCTION LINEAR_INTERPOLATION_1D(A,XMIN,DELX,NX,BADDATA,X,DB,
     $STRICT_CONTROL,EXTEND_X)

C  Thomas Matejka NOAA/NSSL 11 August 1995

C  This function returns an interpolated value from the one-dimensional
C  data field A for the point X.

C  The data field A is defined at NX grid points.  XMIN is the minimum
C  grid coordinate.  The grid points are spaced DELX.

C  When DB is .FALSE., the interpolation is performed directly on the
C  data in A.  When DB is .TRUE., the data in A are assumed to be in
C  decibels, and the interpolation is performed on the linear, not the
C  decibel, values and then converted back to decibels.

C  BADDATA is the value used to indicate missing data in A.  When an
C  interpolated value cannot be obtained, BADDATA is returned.

C  STRICT_CONTROL should be an integer from 1 to 2 that specifies how
C  many of the two surrounding values must be present for the
C  interpolation to proceed.

C  EXTEND_X controls how the interpolation is handled when X lies
C  outside of the domain of A.  If EXTEND_X is .TRUE., then the value of
C  A at the edge of the domain is used.  If EXTEND_X is .FALSE., then
C  BADDATA is returned.

      IMPLICIT NONE
      REAL LINEAR_INTERPOLATION
      LOGICAL DB,EXTEND_X
      INTEGER NX,I,COUNT,STRICT_CONTROL
      INTEGER IX_FACE(2)
      REAL XMIN,DELX,X,BADDATA
      REAL CORNER_DATA_1D(2),X_FACE(2)
      REAL A(NX)

C  Find the grid point numbers and coordinates that bracket X.
      CALL FIND_GRIDBOX_1D(X,XMIN,DELX,NX,IX_FACE,X_FACE)

C  Find the field values at the two surrounding points.
      DO I=1,2
         CORNER_DATA_1D(I)=A(IX_FACE(I))
      ENDDO

C  Count the number of good values.
      COUNT=0
      DO I=1,2
         IF(CORNER_DATA_1D(I).NE.BADDATA)THEN
            COUNT=COUNT+1
         ENDIF
      ENDDO

C  If there are enough good values, proceed with the interpolation.
      IF(COUNT.GE.STRICT_CONTROL)THEN

C  Interpolate.
         LINEAR_INTERPOLATION_1D=LINEAR_INTERPOLATION(X_FACE(1),
     $   CORNER_DATA_1D(1),X_FACE(2),CORNER_DATA_1D(2),X,BADDATA,DB,
     $   EXTEND_X)
      ELSE

C  There are not enough good values.
         LINEAR_INTERPOLATION_1D=BADDATA
      ENDIF

C  Done.
      RETURN
      END
