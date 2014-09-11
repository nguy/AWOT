      REAL FUNCTION INTERP_4D(A,
     $MAXX,MAXY,MAXZ,MAXT,
     $XMIN_ARRAY,DELX_ARRAY,NX_ARRAY,
     $YMIN_ARRAY,DELY_ARRAY,NY_ARRAY,
     $ZMIN_ARRAY,DELZ_ARRAY,NZ_ARRAY,
     $T_DATA,NT,
     $USYS,VSYS,WSYS,
     $BADDATA,
     $X,Y,Z,T,
     $EXTRAP_MODE_X,EXTRAP_MODE_Y,EXTRAP_MODE_Z,EXTRAP_MODE_T,
     $INTERP_ONE_GOOD_X,INTERP_ONE_GOOD_Y,INTERP_ONE_GOOD_Z,
     $INTERP_ONE_GOOD_T)

C  Thomas Matejka NOAA/NSSL 29 April 1998
C=======================================================================
C  This real function performs interpolation on a four-dimensional data
C  field and returns the interpolated value.  If an interpolated value
C  cannot be calculated, it returns BADDATA.  The first three dimensions
C  are assumed to be spatial, and the fourth, time.  The interpolation
C  is conducted in a moving frame of reference.  Linear interpolation is
C  used to interpolate in the spatial dimensions, and cubic spline
C  interpolation is used in the time dimension.

C  Input:

C  A (4d real array 1:MAXX,1:MAXY,1:MAXZ,1:MAXT) specifies the data
C  field to be interpolated.

C  MAXX (integer) specifies the first dimension of A in the calling
C  program.

C  MAXY (integer) specifies the second dimension of A in the calling
C  program.

C  MAXZ (integer) specifies the third dimension of A in the calling
C  program.

C  MAXT (integer) specifies the fourth dimension of A in the calling
C  program.

C  XMIN_ARRAY (1d real array 1:NT).  XMIN_ARRAY(I) specifies the minimum
C  grid coordinate in the first dimension at the Ith data time.

C  YMIN_ARRAY (1d real array 1:NT).  YMIN_ARRAY(I) specifies the minimum
C  grid coordinate in the second dimension at the Ith data time.

C  ZMIN_ARRAY (1d real array 1:NT).  ZMIN_ARRAY(I) specifies the minimum
C  grid coordinate in the third dimension at the Ith data time.

C  DELX_ARRAY (1d real array 1:NT).  DELX_ARRAY(I) specifies the grid
C  point spacing in the first dimension at the Ith data time.

C  DELY_ARRAY (1d real array 1:NT).  DELY_ARRAY(I) specifies the grid
C  point spacing in the second dimension at the Ith data time.

C  DELZ_ARRAY (1d real array 1:NT).  DELZ_ARRAY(I) specifies the grid
C  point spacing in the third dimension at the Ith data time.

C  NX_ARRAY (1d integer array 1:NT).  NX_ARRAY(I) specifies the number
C  of grid points in the first dimension at the Ith data time.

C  NY_ARRAY (1d integer array 1:NT).  NY_ARRAY(I) specifies the number
C  of grid points in the second dimension at the Ith data time.

C  NZ_ARRAY (1d integer array 1:NT).  NZ_ARRAY(I) specifies the number
C  of grid points in the third dimension at the Ith data time.

C  T_DATA (1d real array 1:NT).  T_DATA(I) specifies the Ith data time.

C  NT (integer) specifies the number of data times.

C  USYS (real) specifies the component of velocity in the first
C  dimension of the frame of reference in which to perform the
C  interpolation.

C  VSYS (real) specifies the component of velocity in the second
C  dimension of the frame of reference in which to perform the
C  interpolation.

C  WSYS (real) specifies the component of velocity in the third
C  dimension of the frame of reference in which to perform the
C  interpolation.

C  BADDATA (real) indicates missing values as described.

C  X (real) specifies the coordinate in the first dimension to
C  interpolate to.

C  Y (real) specifies the coordinate in the second dimension to
C  interpolate to.

C  Z (real) specifies the coordinate in the third dimension to
C  interpolate to.

C  T (real) specifies the coordinate in the fourth dimension to
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

C  EXTRAP_MODE_T (character string) controls what the function returns
C  when T lies outside of the grid domain in the fourth dimension.  If
C  it is EXTEND, the nearest value on the edge of the grid is returned.
C  If it is NOTHING, BADDATA is returned.

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

C  INTERP_ONE_GOOD_T (logical) controls what the function returns when T
C  lies inside of the grid domain in the fourth dimension but when only
C  one of the enclosing value exists.  If it is .TRUE., the one good
C  value is returned.  If it is .FALSE., BADDATA is returned.
C=======================================================================
      IMPLICIT NONE
      REAL,EXTERNAL::LINEAR_INTERP,CUBIC_SPLINE_INTERP
      CHARACTER(LEN=7)::EXTRAP_MODE_X,EXTRAP_MODE_Y,EXTRAP_MODE_Z,
     $EXTRAP_MODE_T
      LOGICAL::INTERP_ONE_GOOD_X,INTERP_ONE_GOOD_Y,INTERP_ONE_GOOD_Z,
     $INTERP_ONE_GOOD_T
      INTEGER::MAXX,MAXY,MAXZ,MAXT,NT,
     $I,J,K,L,L_MIN,L_MAX
      INTEGER,DIMENSION(2)::IT_FACE
      INTEGER,DIMENSION(NT)::NX_ARRAY,NY_ARRAY,NZ_ARRAY
      INTEGER,DIMENSION(2,NT)::IX_FACE,IY_FACE,IZ_FACE
      REAL::USYS,VSYS,WSYS,
     $X,Y,Z,T,
     $BADDATA
      REAL,DIMENSION(2)::T_FACE
      REAL,DIMENSION(NT)::CORNER_DATA_1D,
     $X_DATATIME,Y_DATATIME,Z_DATATIME,
     $T_DATA,
     $XMIN_ARRAY,YMIN_ARRAY,ZMIN_ARRAY,
     $DELX_ARRAY,DELY_ARRAY,DELZ_ARRAY
      REAL,DIMENSION(2,NT)::CORNER_DATA_2D,
     $X_FACE,Y_FACE,Z_FACE
      REAL,DIMENSION(2,2,NT)::CORNER_DATA_3D
      REAL,DIMENSION(2,2,2,NT)::CORNER_DATA_4D
      REAL,DIMENSION(MAXX,MAXY,MAXZ,MAXT)::A

C  Loop through the data times.
      DO L=1,NT

C  Calculate the coordinates that are obtained when the position
C  (X,Y,Z,T) is projected in the moving frame of reference onto the data
C  time.
         X_DATATIME(L)=X+USYS*(T_DATA(L)-T)
         Y_DATATIME(L)=Y+VSYS*(T_DATA(L)-T)
         Z_DATATIME(L)=Z+WSYS*(T_DATA(L)-T)

C  Find the grid point numbers and coordinates that enclose the
C  projected positions.
         CALL GRIDBOX_1D(XMIN_ARRAY(L),DELX_ARRAY(L),NX_ARRAY(L),
     $   X_DATATIME(L),IX_FACE(1,L),X_FACE(1,L))
         CALL GRIDBOX_1D(YMIN_ARRAY(L),DELY_ARRAY(L),NY_ARRAY(L),
     $   Y_DATATIME(L),IY_FACE(1,L),Y_FACE(1,L))
         CALL GRIDBOX_1D(ZMIN_ARRAY(L),DELZ_ARRAY(L),NZ_ARRAY(L),
     $   Z_DATATIME(L),IZ_FACE(1,L),Z_FACE(1,L))

C  Find the field values at the eight surrounding points.
         DO K=1,2
            DO J=1,2
               DO I=1,2
                  CORNER_DATA_4D(I,J,K,L)=A(IX_FACE(I,L),IY_FACE(J,L),
     $            IZ_FACE(K,L),L)
               ENDDO
            ENDDO
         ENDDO

C  Interpolate in the first dimension along the eight edges with
C  constant coordinates in the second, third, and fourth dimensions.
         DO K=1,2
            DO J=1,2
               CORNER_DATA_3D(J,K,L)=LINEAR_INTERP(
     $         X_FACE(1,L),CORNER_DATA_4D(1,J,K,L),
     $         X_FACE(2,L),CORNER_DATA_4D(2,J,K,L),
     $         X_DATATIME(L),
     $         BADDATA,
     $         EXTRAP_MODE_X,INTERP_ONE_GOOD_X)
            ENDDO
         ENDDO

C  Interpolate in the second dimension along the four resulting edges
C  with constant coordinates in the third and fourth dimensions.
         DO K=1,2
            CORNER_DATA_2D(K,L)=LINEAR_INTERP(
     $      Y_FACE(1,L),CORNER_DATA_3D(1,K,L),
     $      Y_FACE(2,L),CORNER_DATA_3D(2,K,L),
     $      Y_DATATIME(L),
     $      BADDATA,
     $      EXTRAP_MODE_Y,INTERP_ONE_GOOD_Y)
         ENDDO

C  Interpolate in the third dimension along the two resulting edges with
C  a contant coordinate in the fourth dimension.
         CORNER_DATA_1D(L)=LINEAR_INTERP(
     $   Z_FACE(1,L),CORNER_DATA_2D(1,L),
     $   Z_FACE(2,L),CORNER_DATA_2D(2,L),
     $   Z_DATATIME(L),
     $   BADDATA,
     $   EXTRAP_MODE_Z,INTERP_ONE_GOOD_Z)
      ENDDO

C  Find the data times that enclose T.
      CALL IRREGULAR_GRIDBOX_1D(NT,T_DATA,T,IT_FACE,T_FACE)

C  Proceed with the interpolation based on by how many good values T is
C  enclosed.
      IF(CORNER_DATA_1D(IT_FACE(1)).EQ.BADDATA.AND.
     $CORNER_DATA_1D(IT_FACE(2)).EQ.BADDATA)THEN
         INTERP_4D=BADDATA
      ELSEIF(CORNER_DATA_1D(IT_FACE(1)).NE.BADDATA.AND.
     $CORNER_DATA_1D(IT_FACE(2)).EQ.BADDATA)THEN
         IF(INTERP_ONE_GOOD_T)THEN
            INTERP_4D=CORNER_DATA_1D(IT_FACE(1))
         ELSE
            INTERP_4D=BADDATA
         ENDIF
      ELSEIF(CORNER_DATA_1D(IT_FACE(1)).EQ.BADDATA.AND.
     $CORNER_DATA_1D(IT_FACE(2)).NE.BADDATA)THEN
         IF(INTERP_ONE_GOOD_T)THEN
            INTERP_4D=CORNER_DATA_1D(IT_FACE(2))
         ELSE
            INTERP_4D=BADDATA
         ENDIF
      ELSE
         DO L=IT_FACE(1),1,-1
            IF(CORNER_DATA_1D(L).NE.BADDATA)THEN
               L_MIN=L
            ELSE
               EXIT
            ENDIF
         ENDDO
         DO L=IT_FACE(2),NT
            IF(CORNER_DATA_1D(L).NE.BADDATA)THEN
               L_MAX=L
            ELSE
               EXIT
            ENDIF
         ENDDO
         INTERP_4D=CUBIC_SPLINE_INTERP(L_MAX-L_MIN+1,T_DATA(L_MIN),
     $   CORNER_DATA_1D(L_MIN),T,
     $   BADDATA,
     $   EXTRAP_MODE_T)
      ENDIF

      END FUNCTION INTERP_4D
