      REAL FUNCTION LINEAR_INTERPOLATION_3D(A,MAXX,MAXY,MAXZ,XMIN,YMIN,
     $ZMIN,DELX,DELY,DELZ,NX,NY,NZ,BADDATA,X,Y,Z,DB,STRICT_CONTROL,
     $EXTEND_X,EXTEND_Y,EXTEND_Z)

C  Thomas Matejka NOAA/NSSL 11 August 1995

C  This function returns an interpolated value from the
C  three-dimensional data field A for the point (X,Y,Z).

C  MAXX, MAXY, and MAXZ are the first, second, and third dimensions of A
C  in the calling program.

C  The data field A is defined at NX, NY, NZ, and NT grid points in the
C  first, second, third, and fourth dimensions.  XMIN, YMIN, and ZMIN
C  are the minimum grid coordinates in the first, second, and third
C  dimensions.  The grid points are spaced DELX, DELY, and DELZ in the
C  first, second, and third dimensions.

C  When DB is .FALSE., the interpolation is performed directly on the
C  data in A.  When DB is .TRUE., the data in A are assumed to be in
C  decibels, and the interpolation is performed on the linear, not the
C  decibel, values and then converted back to decibels.

C  BADDATA is the value used to indicate missing data in A.  When an
C  interpolated value cannot be obtained, BADDATA is returned.

C  STRICT_CONTROL should be an integer from 1 to 8 that specifies how
C  many of the eight surrounding values must be present for the
C  interpolation to proceed.

C  EXTEND_X controls how the interpolation is handled when X lies
C  outside of the domain of A.  If EXTEND_X is .TRUE., then the value of
C  A at the edge of the domain is used.  If EXTEND_X is .FALSE., then
C  BADDATA is returned

C  EXTEND_Y controls how the interpolation is handled when Y lies
C  outside of the domain of A.  If EXTEND_Y is .TRUE., then the value of
C  A at the edge of the domain is used.  If EXTEND_Y is .FALSE., then
C  BADDATA is returned

C  EXTEND_Z controls how the interpolation is handled when Z lies
C  outside of the domain of A.  If EXTEND_Z is .TRUE., then the value of
C  A at the edge of the domain is used.  If EXTEND_Z is .FALSE., then
C  BADDATA is returned

      IMPLICIT NONE
      REAL LINEAR_INTERPOLATION
      LOGICAL DB,EXTEND_X,EXTEND_Y,EXTEND_Z
      INTEGER MAXX,MAXY,MAXZ,NX,NY,NZ,I,J,K,COUNT,STRICT_CONTROL
      INTEGER IX_FACE(2),IY_FACE(2),IZ_FACE(2)
      REAL XMIN,YMIN,ZMIN,DELX,DELY,DELZ,X,Y,Z,BADDATA
      REAL CORNER_DATA_1D(2),X_FACE(2),Y_FACE(2),Z_FACE(2)
      REAL CORNER_DATA_2D(2,2)
      REAL CORNER_DATA_3D(2,2,2)
      REAL A(MAXX,MAXY,MAXZ)

C  Find the grid point numbers and coordinates that bracket (X,Y,Z).
      CALL FIND_GRIDBOX_1D(X,XMIN,DELX,NX,IX_FACE,X_FACE)
      CALL FIND_GRIDBOX_1D(Y,YMIN,DELY,NY,IY_FACE,Y_FACE)
      CALL FIND_GRIDBOX_1D(Z,ZMIN,DELZ,NZ,IZ_FACE,Z_FACE)

C  Find the field values at the eight surrounding points.
      DO K=1,2
         DO J=1,2
            DO I=1,2
               CORNER_DATA_3D(I,J,K)=A(IX_FACE(I),IY_FACE(J),IZ_FACE(K))
            ENDDO
         ENDDO
      ENDDO

C  Count the number of good values.
      COUNT=0
      DO K=1,2
         DO J=1,2
            DO I=1,2
               IF(CORNER_DATA_3D(I,J,K).NE.BADDATA)THEN
                  COUNT=COUNT+1
               ENDIF
            ENDDO
         ENDDO
      ENDDO

C  If there are enough good values, proceed with the interpolation.
      IF(COUNT.GE.STRICT_CONTROL)THEN

C  Interpolate in the first dimension along the four edges with constant
C  coordinates in the second and third dimensions.
         DO K=1,2
            DO J=1,2
               CORNER_DATA_2D(J,K)=LINEAR_INTERPOLATION(X_FACE(1),
     $         CORNER_DATA_3D(1,J,K),X_FACE(2),CORNER_DATA_3D(2,J,K),X,
     $         BADDATA,DB,EXTEND_X)
            ENDDO
         ENDDO

C  Interpolate in the second dimension along the two resulting edges
C  with constant coordinates in the third dimension.
         DO K=1,2
            CORNER_DATA_1D(K)=LINEAR_INTERPOLATION(Y_FACE(1),
     $      CORNER_DATA_2D(1,K),Y_FACE(2),CORNER_DATA_2D(2,K),Y,BADDATA,
     $      DB,EXTEND_Y)
         ENDDO

C  Interpolate in the third dimension along the one resulting edge.
         LINEAR_INTERPOLATION_3D=LINEAR_INTERPOLATION(Z_FACE(1),
     $   CORNER_DATA_1D(1),Z_FACE(2),CORNER_DATA_1D(2),Z,BADDATA,DB,
     $   EXTEND_Z)
      ELSE

C  There are not enought good values.
         LINEAR_INTERPOLATION_3D=BADDATA
      ENDIF

C  Done.
      RETURN
      END
