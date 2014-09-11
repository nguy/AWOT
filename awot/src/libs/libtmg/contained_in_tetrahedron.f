      LOGICAL FUNCTION CONTAINED_IN_TETRAHEDRON(X_A,Y_A,Z_A,X_B,Y_B,Z_B,
     $X_C,Y_C,Z_C,X_D,Y_D,Z_D,X_F,Y_F,Z_F)

C  Thomas Matejka NOAA/NSSL 28 December 1994

C  This function returns .TRUE. if and only if a test point is contained
C  within a tetrahedron or is on its surface.

C  Input:

C  X_A is a real variable that indicates the first coordinate of the
C  first vertex of the tetrahedron.

C  Y_A is a real variable that indicates the second coordinate of the
C  first vertex of the tetrahedron.

C  Z_A is a real variable that indicates the third coordinate of the
C  first vertex of the tetrahedron.

C  X_B is a real variable that indicates the first coordinate of the
C  second vertex of the tetrahedron.

C  Y_B is a real variable that indicates the second coordinate of the
C  second vertex of the tetrahedron.

C  Z_B is a real variable that indicates the third coordinate of the
C  second vertex of the tetrahedron.

C  X_C is a real variable that indicates the first coordinate of the
C  third vertex of the tetrahedron.

C  Y_C is a real variable that indicates the second coordinate of the
C  third vertex of the tetrahedron.

C  Z_C is a real variable that indicates the third coordinate of the
C  third vertex of the tetrahedron.

C  X_D is a real variable that indicates the first coordinate of the
C  fourth vertex of the tetrahedron.

C  Y_D is a real variable that indicates the second coordinate of the
C  fourth vertex of the tetrahedron.

C  Z_D is a real variable that indicates the third coordinate of the
C  fourth vertex of the tetrahedron.

C  X_F is a real variable that indicates the first coordinate of the
C  test point.

C  Y_F is a real variable that indicates the second coordinate of the
C  test point.

C  Z_F is a real variable that indicates the third coordinate of the
C  test point.

      IMPLICIT NONE
      INTEGER SAME_SIDE_OF_PLANE
      REAL X_A,Y_A,Z_A,X_B,Y_B,Z_B,X_C,Y_C,Z_C,X_D,Y_D,Z_D,X_F,Y_F,Z_F

C  The test point is contained within the tetrahedron if and only if it
C  lies on the same side of each face as the other vertex.
      IF(SAME_SIDE_OF_PLANE(X_A,Y_A,Z_A,X_B,Y_B,Z_B,X_C,Y_C,Z_C,X_D,Y_D,
     $Z_D,X_F,Y_F,Z_F).LT.0)THEN
         CONTAINED_IN_TETRAHEDRON=.FALSE.
         RETURN
      ENDIF
      IF(SAME_SIDE_OF_PLANE(X_B,Y_B,Z_B,X_C,Y_C,Z_C,X_D,Y_D,Z_D,X_A,Y_A,
     $Z_A,X_F,Y_F,Z_F).LT.0)THEN
         CONTAINED_IN_TETRAHEDRON=.FALSE.
         RETURN
      ENDIF
      IF(SAME_SIDE_OF_PLANE(X_C,Y_C,Z_C,X_D,Y_D,Z_D,X_A,Y_A,Z_A,X_B,Y_B,
     $Z_B,X_F,Y_F,Z_F).LT.0)THEN
         CONTAINED_IN_TETRAHEDRON=.FALSE.
         RETURN
      ENDIF
      IF(SAME_SIDE_OF_PLANE(X_D,Y_D,Z_D,X_A,Y_A,Z_A,X_B,Y_B,Z_B,X_C,Y_C,
     $Z_C,X_F,Y_F,Z_F).LT.0)THEN
         CONTAINED_IN_TETRAHEDRON=.FALSE.
         RETURN
      ENDIF

C  Done.
      CONTAINED_IN_TETRAHEDRON=.TRUE.
      RETURN
      END
