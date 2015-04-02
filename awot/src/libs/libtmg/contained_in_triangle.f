      LOGICAL FUNCTION CONTAINED_IN_TRIANGLE(X_A,Y_A,Z_A,X_B,Y_B,Z_B,
     $X_D,Y_D,Z_D,X_F,Y_F,Z_F)

C  Thomas Matejka NOAA/NSSL 28 December 1994

C  This function returns .TRUE. if and only if a test point is contained
C  within a triangle or is on its surface

C  Input:

C  X_A is a real variable that indicates the first coordinate of the
C  first vertex of the triangle.

C  Y_A is a real variable that indicates the second coordinate of the
C  first vertex of the triangle.

C  Z_A is a real variable that indicates the third coordinate of the
C  first vertex of the triangle.

C  X_B is a real variable that indicates the first coordinate of the
C  second vertex of the triangle.

C  Y_B is a real variable that indicates the second coordinate of the
C  second vertex of the triangle.

C  Z_B is a real variable that indicates the third coordinate of the
C  second vertex of the triangle.

C  X_D is a real variable that indicates the first coordinate of the
C  third vertex of the triangle.

C  Y_D is a real variable that indicates the second coordinate of the
C  third vertex of the triangle.

C  Z_D is a real variable that indicates the third coordinate of the
C  third vertex of the triangle.

C  X_F is a real variable that indicates the first coordinate of the
C  test point.

C  Y_F is a real variable that indicates the second coordinate of the
C  test point.

C  Z_F is a real variable that indicates the third coordinate of the
C  test point.

      IMPLICIT NONE
      INTEGER SAME_SIDE_OF_LINE
      REAL X_A,Y_A,Z_A,X_B,Y_B,Z_B,X_D,Y_D,Z_D,X_F,Y_F,Z_F

C  The test point is contained within the triangle if and only if it
C  lies on the same side of each face as the other vertex.
      IF(SAME_SIDE_OF_LINE(X_A,Y_A,Z_A,X_B,Y_B,Z_B,X_D,Y_D,Z_D,X_F,Y_F,
     $Z_F).LT.0)THEN
         CONTAINED_IN_TRIANGLE=.FALSE.
         RETURN
      ENDIF
      IF(SAME_SIDE_OF_LINE(X_B,Y_B,Z_B,X_D,Y_D,Z_D,X_A,Y_A,Z_A,X_F,Y_F,
     $Z_F).LT.0)THEN
         CONTAINED_IN_TRIANGLE=.FALSE.
         RETURN
      ENDIF
      IF(SAME_SIDE_OF_LINE(X_D,Y_D,Z_D,X_A,Y_A,Z_A,X_B,Y_B,Z_B,X_F,Y_F,
     $Z_F).LT.0)THEN
         CONTAINED_IN_TRIANGLE=.FALSE.
         RETURN
      ENDIF

C  Done.
      CONTAINED_IN_TRIANGLE=.TRUE.
      RETURN
      END
