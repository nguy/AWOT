      REAL FUNCTION INTERP_1D(A,
     $XMIN,DELX,NX,
     $BADDATA,
     $X,
     $EXTRAP_MODE_X,
     $INTERP_ONE_GOOD_X)

C  Thomas Matejka NOAA/NSSL 23 June 1998
C=======================================================================
C  This real function performs linear interpolation on a one-dimensional
C  data field and returns the interpolated value.  If an interpolated
C  value cannot be calculated, it returns BADDATA.

C  Input:

C  A (1d real array 1:NX) specifies the data field to be interpolated.

C  XMIN (real).  XMIN specifies the minimum grid coordinate.

C  DELX (real).  DELX specifies the grid point spacing.

C  NX (integer).  NX specifies the number of grid points.

C  BADDATA (real) indicates missing values as described.

C  X (real) specifies the coordinate to interpolate to.

C  EXTRAP_MODE_X (character string) controls what the function returns
C  when X lies outside of the grid domain.  If it is LINEAR,
C  extrapolation is performed with the linear trend used for
C  interpolation.  If it is EXTEND, the nearest value on the edge of the
C  grid is returned.  If it is NOTHING, BADDATA is returned.

C  INTERP_ONE_GOOD_X (logical) controls what the function returns when X
C  lies inside of the grid domain but when only one of the enclosing
C  value exists.  If it is .TRUE., the one good value is returned.  If
C  it is .FALSE., BADDATA is returned.
C=======================================================================
      IMPLICIT NONE
      REAL,EXTERNAL::LINEAR_INTERP
      CHARACTER(LEN=7)::EXTRAP_MODE_X
      LOGICAL::INTERP_ONE_GOOD_X
      INTEGER::NX,I
      INTEGER,DIMENSION(2)::IX_FACE
      REAL::XMIN,
     $DELX,
     $X,
     $BADDATA
      REAL,DIMENSION(2)::CORNER_DATA_1D,
     $X_FACE
      REAL,DIMENSION(NX)::A

C  Find the grid point numbers and coordinates that enclose X.
      CALL GRIDBOX_1D(XMIN,DELX,NX,X,IX_FACE,X_FACE)

C  Find the field values at the two surrounding points.
      DO I=1,2
         CORNER_DATA_1D(I)=A(IX_FACE(I))
      ENDDO

C  Interpolate.
      INTERP_1D=LINEAR_INTERP(X_FACE(1),CORNER_DATA_1D(1),
     $X_FACE(2),CORNER_DATA_1D(2),
     $X,
     $BADDATA,
     $EXTRAP_MODE_X,INTERP_ONE_GOOD_X)

      END FUNCTION INTERP_1D
