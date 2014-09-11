      REAL FUNCTION INTERP_3D(A,
     $MAXX,MAXY,MAXZ,
     $XMIN,DELX,NX,
     $YMIN,DELY,NY,
     $ZMIN,DELZ,NZ,
     $BADDATA,
     $X,Y,Z,
     $EXTRAP_MODE_X,EXTRAP_MODE_Y,EXTRAP_MODE_Z,
     $INTERP_ONE_GOOD_X,INTERP_ONE_GOOD_Y,INTERP_ONE_GOOD_Z)

C  Thomas Matejka NOAA/NSSL 23 June 1998
C=======================================================================
C  This real function performs linear interpolation on a
C  three-dimensional data field and returns the interpolated value.  If
C  an interpolated value cannot be calculated, it returns BADDATA.

C  Input:

C  A (3d real array 1:MAXX,1:MAXY,1:MAXZ) specifies the data field to be
C  interpolated.

C  MAXX (integer) specifies the first dimension of A in the calling
C  program.

C  MAXY (integer) specifies the second dimension of A in the calling
C  program.

C  MAXZ (integer) specifies the third dimension of A in the calling
C  program.

C  XMIN (real).  XMIN specifies the minimum grid coordinate in the first
C  dimension.

C  YMIN (real).  YMIN specifies the minimum grid coordinate in the
C  second dimension.

C  ZMIN (real).  ZMIN specifies the minimum grid coordinate in the third
C  dimension.

C  DELX (real).  DELX specifies the grid point spacing in the first
C  dimension.

C  DELY (real).  DELY specifies the grid point spacing in the second
C  dimension.

C  DELZ (real).  DELZ specifies the grid point spacing in the third
C  dimension.

C  NX (integer).  NX specifies the number of grid points in the first
C  dimension.

C  NY (integer).  NY specifies the number of grid points in the second
C  dimension.

C  NZ (integer).  NZ specifies the number of grid points in the third
C  dimension.

C  BADDATA (real) indicates missing values as described.

C  X (real) specifies the coordinate in the first dimension to
C  interpolate to.

C  Y (real) specifies the coordinate in the second dimension to
C  interpolate to.

C  Z (real) specifies the coordinate in the third dimension to
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

C  EXTRAP_MODE_Z (character string) controls what the function returns
C  when Z lies outside of the grid domain in the third dimension.  If it
C  is LINEAR, extrapolation is performed with the linear trend used for
C  interpolation.  If it is EXTEND, the nearest value on the edge of the
C  grid is returned.  If it is NOTHING, BADDATA is returned.

C  INTERP_ONE_GOOD_X (logical) controls what the function returns when X
C  lies inside of the grid domain in the first dimension but when only
C  one of the enclosing value exists.  If it is .TRUE., the one good
C  value is returned.  If it is .FALSE., BADDATA is returned.

C  INTERP_ONE_GOOD_Y (logical) controls what the function returns when Y
C  lies inside of the grid domain in the second dimension but when only
C  one of the enclosing value exists.  If it is .TRUE., the one good
C  value is returned.  If it is .FALSE., BADDATA is returned.

C  INTERP_ONE_GOOD_Z (logical) controls what the function returns when Z
C  lies inside of the grid domain in the third dimension but when only
C  one of the enclosing value exists.  If it is .TRUE., the one good
C  value is returned.  If it is .FALSE., BADDATA is returned.
C=======================================================================
      IMPLICIT NONE
      REAL,EXTERNAL::LINEAR_INTERP
      CHARACTER(LEN=7)::EXTRAP_MODE_X,EXTRAP_MODE_Y,EXTRAP_MODE_Z
      LOGICAL::INTERP_ONE_GOOD_X,INTERP_ONE_GOOD_Y,INTERP_ONE_GOOD_Z
      INTEGER::MAXX,MAXY,MAXZ,
     $NX,NY,NZ,I,J,K
      INTEGER,DIMENSION(2)::IX_FACE,IY_FACE,IZ_FACE
      REAL::XMIN,YMIN,ZMIN,
     $DELX,DELY,DELZ,
     $X,Y,Z,
     $BADDATA
      REAL,DIMENSION(2)::CORNER_DATA_1D,
     $X_FACE,Y_FACE,Z_FACE
      REAL,DIMENSION(2,2)::CORNER_DATA_2D
      REAL,DIMENSION(2,2,2)::CORNER_DATA_3D
      REAL,DIMENSION(MAXX,MAXY,MAXZ)::A

C  Find the grid point numbers and coordinates that enclose (X,Y,Z).
      CALL GRIDBOX_1D(XMIN,DELX,NX,X,IX_FACE,X_FACE)
      CALL GRIDBOX_1D(YMIN,DELY,NY,Y,IY_FACE,Y_FACE)
      CALL GRIDBOX_1D(ZMIN,DELZ,NZ,Z,IZ_FACE,Z_FACE)

C  Find the field values at the eight surrounding points.
      DO K=1,2
         DO J=1,2
            DO I=1,2
               CORNER_DATA_3D(I,J,K)=A(IX_FACE(I),IY_FACE(J),IZ_FACE(K))
            ENDDO
         ENDDO
      ENDDO

C  Interpolate in the first dimension along the four edges with constant
C  coordinates in the second and third dimensions.
      DO K=1,2
         DO J=1,2
            CORNER_DATA_2D(J,K)=LINEAR_INTERP(
     $      X_FACE(1),CORNER_DATA_3D(1,J,K),
     $      X_FACE(2),CORNER_DATA_3D(2,J,K),
     $      X,
     $      BADDATA,
     $      EXTRAP_MODE_X,INTERP_ONE_GOOD_X)
         ENDDO
      ENDDO

C  Interpolate in the second dimension along the two resulting edges
C  with constant coordinates in the third dimension.
      DO K=1,2
         CORNER_DATA_1D(K)=LINEAR_INTERP(Y_FACE(1),CORNER_DATA_2D(1,K),
     $   Y_FACE(2),CORNER_DATA_2D(2,K),
     $   Y,
     $   BADDATA,
     $   EXTRAP_MODE_Y,INTERP_ONE_GOOD_Y)
      ENDDO

C  Interpolate in the third dimension.
      INTERP_3D=LINEAR_INTERP(Z_FACE(1),CORNER_DATA_1D(1),
     $Z_FACE(2),CORNER_DATA_1D(2),
     $Z,
     $BADDATA,
     $EXTRAP_MODE_Z,INTERP_ONE_GOOD_Z)

      END FUNCTION INTERP_3D
