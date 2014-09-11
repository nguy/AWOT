      REAL FUNCTION INTERP_2D(A,
     $MAXX,MAXY,
     $XMIN,DELX,NX,
     $YMIN,DELY,NY,
     $BADDATA,
     $X,Y,
     $EXTRAP_MODE_X,EXTRAP_MODE_Y,
     $INTERP_ONE_GOOD_X,INTERP_ONE_GOOD_Y)

C  Thomas Matejka NOAA/NSSL 23 June 1998
C=======================================================================
C  This real function performs linear interpolation on a two-dimensional
C  data field and returns the interpolated value.  If an interpolated
C  value cannot be calculated, it returns BADDATA.

C  Input:

C  A (2d real array 1:MAXX,1:MAXY) specifies the data field to be
C  interpolated.

C  MAXX (integer) specifies the first dimension of A in the calling
C  program.

C  MAXY (integer) specifies the second dimension of A in the calling
C  program.

C  XMIN (real).  XMIN specifies the minimum grid coordinate in the first
C  dimension.

C  YMIN (real).  YMIN specifies the minimum grid coordinate in the
C  second dimension.

C  DELX (real).  DELX specifies the grid point spacing in the first
C  dimension.

C  DELY (real).  DELY specifies the grid point spacing in the second
C  dimension.

C  NX (integer).  NX specifies the number of grid points in the first
C  dimension.

C  NY (integer).  NY specifies the number of grid points in the second
C  dimension.

C  BADDATA (real) indicates missing values as described.

C  X (real) specifies the coordinate in the first dimension to
C  interpolate to.

C  Y (real) specifies the coordinate in the second dimension to
C  interpolate to.

C  EXTRAP_MODE_X (character string) controls what the function returns
C  when X lies outside of the grid domain in the first dimension.  If it
C  is LINEAR, extrapolation is performed with the linear trend used for
C  interpolation.  If it is EXTEND, the nearest value on the edge of the
C  grid is returned.  If it is NOTHING, BADDATA is returned.

C  EXTRAP_MODE_Y (character string) controls what the function returns
C  when Y lies outside of the grid domain in the second dimension.  If
C  it is LINEAR, extrapolation is performed with the linear trend used
C  for interpolation.  If it is EXTEND, the nearest value on the edge of
C  the grid is returned.  If it is NOTHING, BADDATA is returned.

C  INTERP_ONE_GOOD_X (logical) controls what the function returns when X
C  lies inside of the grid domain in the first dimension but when only
C  one of the enclosing value exists.  If it is .TRUE., the one good
C  value is returned.  If it is .FALSE., BADDATA is returned.

C  INTERP_ONE_GOOD_Y (logical) controls what the function returns when Y
C  lies inside of the grid domain in the second dimension but when only
C  one of the enclosing value exists.  If it is .TRUE., the one good
C  value is returned.  If it is .FALSE., BADDATA is returned.
C=======================================================================
      IMPLICIT NONE
      REAL,EXTERNAL::LINEAR_INTERP
      CHARACTER(LEN=7)::EXTRAP_MODE_X,EXTRAP_MODE_Y
      LOGICAL::INTERP_ONE_GOOD_X,INTERP_ONE_GOOD_Y
      INTEGER::MAXX,MAXY,
     $NX,NY,I,J
      INTEGER,DIMENSION(2)::IX_FACE,IY_FACE
      REAL::XMIN,YMIN,
     $DELX,DELY,
     $X,Y,
     $BADDATA
      REAL,DIMENSION(2)::CORNER_DATA_1D,
     $X_FACE,Y_FACE
      REAL,DIMENSION(2,2)::CORNER_DATA_2D
      REAL,DIMENSION(MAXX,MAXY)::A

C  Find the grid point numbers and coordinates that enclose (X,Y).
      CALL GRIDBOX_1D(XMIN,DELX,NX,X,IX_FACE,X_FACE)
      CALL GRIDBOX_1D(YMIN,DELY,NY,Y,IY_FACE,Y_FACE)

C  Find the field values at the four surrounding points.
      DO J=1,2
         DO I=1,2
            CORNER_DATA_2D(I,J)=A(IX_FACE(I),IY_FACE(J))
         ENDDO
      ENDDO

C  Interpolate in the first dimension along the two edges with constant
C  coordinates in the second dimension.
      DO J=1,2
         CORNER_DATA_1D(J)=LINEAR_INTERP(X_FACE(1),CORNER_DATA_2D(1,J),
     $   X_FACE(2),CORNER_DATA_2D(2,J),
     $   X,
     $   BADDATA,
     $   EXTRAP_MODE_X,INTERP_ONE_GOOD_X)
      ENDDO

C  Interpolate in the second dimension.
      INTERP_2D=LINEAR_INTERP(Y_FACE(1),CORNER_DATA_1D(1),
     $Y_FACE(2),CORNER_DATA_1D(2),
     $Y,
     $BADDATA,
     $EXTRAP_MODE_Y,INTERP_ONE_GOOD_Y)

      END FUNCTION INTERP_2D
