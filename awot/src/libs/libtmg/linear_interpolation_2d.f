      REAL FUNCTION LINEAR_INTERPOLATION_2D(A,MAXX,MAXY,XMIN,YMIN,DELX,
     $DELY,NX,NY,BADDATA,X,Y,DB,STRICT_CONTROL,EXTEND_X,EXTEND_Y)

C  Thomas Matejka NOAA/NSSL 11 August 1995

C  This function returns an interpolate value from the two-dimensional
C  data field A for the point (X,Y).

C  MAXX and MAXY are the first and second dimensions of A in the calling
C  program.

C  The data field A is defined at NX and NY grid points in the first and
C  second dimensions.  XMIN and YMIN are the minimum grid coordinates in
C  the first and second dimensions.  The grid points are spaced DELX and
C  DELY in the first and second dimensions.

C  When DB is .FALSE., the interpolation is performed directly on the
C  data in A.  When DB is .TRUE., the data in A are assumed to be in
C  decibels, and the interpolation is performed on the linear, not the
C  decibel, values and then converted back to decibels.

C  BADDATA is the value used to indicate missing data in A.  When an
C  interpolated value cannot be obtained, BADDATA is returned.

C  STRICT_CONTROL should be an integer from 1 to 4 that specifies how
C  many of the four surrounding values must be present for the
C  interpolation to proceed.

C  EXTEND_X controls how the interpolation is handled when X lies
C  outside of the domain of A.  If EXTEND_X is .TRUE., then the value of
C  A at the edge of the domain is used.  If EXTEND_X is .FALSE., then
C  BADDATA is returned

C  EXTEND_Y controls how the interpolation is handled when Y lies
C  outside of the domain of A.  If EXTEND_Y is .TRUE., then the value of
C  A at the edge of the domain is used.  If EXTEND_Y is .FALSE., then
C  BADDATA is returned

      IMPLICIT NONE
      REAL LINEAR_INTERPOLATION
      LOGICAL DB,EXTEND_X,EXTEND_Y
      INTEGER MAXX,MAXY,NX,NY,I,J,COUNT,STRICT_CONTROL
      INTEGER IX_FACE(2),IY_FACE(2)
      REAL XMIN,YMIN,DELX,DELY,X,Y,BADDATA
      REAL CORNER_DATA_1D(2),X_FACE(2),Y_FACE(2)
      REAL CORNER_DATA_2D(2,2)
      REAL A(MAXX,MAXY)

C  Find the grid point numbers and coordinates that bracket (X,Y).
      CALL FIND_GRIDBOX_1D(X,XMIN,DELX,NX,IX_FACE,X_FACE)
      CALL FIND_GRIDBOX_1D(Y,YMIN,DELY,NY,IY_FACE,Y_FACE)

C  Find the field values at the four surrounding points.
      DO J=1,2
         DO I=1,2
            CORNER_DATA_2D(I,J)=A(IX_FACE(I),IY_FACE(J))
         ENDDO
      ENDDO

C  Count the number of good values.
      COUNT=0
      DO J=1,2
         DO I=1,2
            IF(CORNER_DATA_2D(I,J).NE.BADDATA)THEN
               COUNT=COUNT+1
            ENDIF
         ENDDO
      ENDDO

C  If there are enough good values, proceed with the interpolation.
      IF(COUNT.GE.STRICT_CONTROL)THEN

C  Interpolate in the first dimension along the two edges with constant
C  coordinates in the second dimension.
         DO J=1,2
            CORNER_DATA_1D(J)=LINEAR_INTERPOLATION(X_FACE(1),
     $      CORNER_DATA_2D(1,J),X_FACE(2),CORNER_DATA_2D(2,J),X,BADDATA,
     $      DB,EXTEND_X)
         ENDDO

C  Interpolate in the second dimension along the one resulting edge.
         LINEAR_INTERPOLATION_2D=LINEAR_INTERPOLATION(Y_FACE(1),
     $   CORNER_DATA_1D(1),Y_FACE(2),CORNER_DATA_1D(2),Y,BADDATA,DB,
     $   EXTEND_Y)
      ELSE

C  There are not enough good values.
         LINEAR_INTERPOLATION_2D=BADDATA
      ENDIF

C  Done.
      RETURN
      END
