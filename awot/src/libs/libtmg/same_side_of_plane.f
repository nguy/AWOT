      INTEGER FUNCTION SAME_SIDE_OF_PLANE(X_A,Y_A,Z_A,X_B,Y_B,Z_B,X_C,
     $Y_C,Z_C,X_D,Y_D,Z_D,X_F,Y_F,Z_F)

C  Thomas Matejka NOAA/NSSL 13 December 1996

C  This function returns 1 if and only if a test point is on the same
C  side of the plane as a reference point.  The function returns 0 if
C  and only if the test point is in the plane.  The function returns -1
C  if and only if the test point is on the opposite side of the plane as
C  the reference point.  The function returns -2 if and only if the test
C  was indeterminate because of degeneracy.

C  Input:

C  X_A is a real variable that indicates the first coordinate of the
C  first point in the plane.

C  Y_A is a real variable that indicates the second coordinate of the
C  first point in the plane.

C  Z_A is a real variable that indicates the third coordinate of the
C  first point in the plane.

C  X_B is a real variable that indicates the first coordinate of the
C  second point in the plane.

C  Y_B is a real variable that indicates the second coordinate of the
C  second point in the plane.

C  Z_B is a real variable that indicates the third coordinate of the
C  second point in the plane.

C  X_C is a real variable that indicates the first coordinate of the
C  third point in the plane.

C  Y_C is a real variable that indicates the second coordinate of the
C  third point in the plane.

C  Z_C is a real variable that indicates the third coordinate of the
C  third point in the plane.

C  X_D is a real variable that indicates the first coordinate of the
C  reference point.  The reference point should not be in the plane.

C  Y_D is a real variable that indicates the second coordinate of the
C  reference point.  The reference point should not be in the plane.

C  Z_D is a real variable that indicates the third coordinate of the
C  reference point.  The reference point should not be in the plane.

C  X_F is a real variable that indicates the first coordinate of the
C  test point.

C  Y_F is a real variable that indicates the second coordinate of the
C  test point.

C  Z_F is a real variable that indicates the third coordinate of the
C  test point.

      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      REAL::X_A,Y_A,Z_A,X_B,Y_B,Z_B,X_C,Y_C,Z_C,X_D,Y_D,Z_D,X_F,Y_F,Z_F,
     $BX,BY,BZ,CX,CY,CZ,DX,DY,DZ,NX,NY,NZ,MX,MY,MZ,ALPHA,EX,EY,EZ,FX,FY,
     $FZ,Q,VECTORS_DISTINCT

C  Calculate vector components from the first point in the plane.
      BX=X_B-X_A
      BY=Y_B-Y_A
      BZ=Z_B-Z_A
      CX=X_C-X_A
      CY=Y_C-Y_A
      CZ=Z_C-Z_A
      DX=X_D-X_A
      DY=Y_D-Y_A
      DZ=Z_D-Z_A
      FX=X_F-X_A
      FY=Y_F-Y_A
      FZ=Z_F-Z_A

C  Check whether the two vectors in the plane are sufficiently distinct
C  to define the plane.
      VECTORS_DISTINCT=SQRT(((BY*CZ-BZ*CY)**2+(BX*CZ-BZ*CX)**2+
     $(BX*CY-BY*CX)**2)/((BX**2+BY**2+BZ**2)*(CX**2+CY**2+CZ**2)))
      IF(VECTORS_DISTINCT.LT.PLANE_THRESHOLD)THEN
         SAME_SIDE_OF_PLANE=-2
         RETURN
      ENDIF

C  Calculate the components of a vector normal to the plane.
      NX=BY*CZ-BZ*CY
      NY=-BX*CZ+BZ*CX
      NZ=BX*CY-BY*CX

C  Check whether the test point is on the plane.
      MX=BY*FZ-BZ*FY
      MY=-BX*FZ+BZ*FX
      MZ=BX*FY-BY*FX 
      IF(MY*NZ.EQ.MZ*NY.AND.
     $MX*NZ.EQ.MZ*NX.AND.
     $MX*NY.EQ.MY*NX)THEN
         SAME_SIDE_OF_PLANE=0
         RETURN
      ENDIF

C  Calculate the components of a vector in the plane to the point
C  closest to the reference point.
      ALPHA=(DX*NX+DY*NY+DZ*NZ)/(NX**2+NY**2+NZ**2)
      EX=DX-ALPHA*NX
      EY=DY-ALPHA*NY
      EZ=DZ-ALPHA*NZ

C  Check whether the reference point is on the plane.  If it is, the
C  result is indeterminate.
      IF(DX.EQ.EX.AND.
     $DY.EQ.EY.AND.
     $DZ.EQ.EZ)THEN
         SAME_SIDE_OF_PLANE=-2
         RETURN
      ENDIF

C  Compare vectors to the test point and reference point.
      Q=(DX-EX)*(FX-EX)+(DY-EY)*(FY-EY)+(DZ-EZ)*(FZ-EZ)
      IF(Q.GT.0.)THEN
         SAME_SIDE_OF_PLANE=1
      ELSEIF(Q.LT.0.)THEN
         SAME_SIDE_OF_PLANE=-1
      ELSE
         SAME_SIDE_OF_PLANE=0
      ENDIF

C  Done.
      RETURN
      END
