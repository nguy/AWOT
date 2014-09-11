      INTEGER FUNCTION SAME_SIDE_OF_LINE(X_A,Y_A,Z_A,X_B,Y_B,Z_B,X_D,
     $Y_D,Z_D,X_F,Y_F,Z_F)

C  Thomas Matejka NOAA/NSSL 1 September 1995

C  This function returns 1 if and only if a test point is on the same
C  side of the line as a reference point.  The function returns 0 if and
C  only if the test point is on the line.  The function returns -1 if
C  and only if the test point is on the opposite side of the line as the
C  reference point.  The function returns -2 if and only if the test was
C  indeterminate because of degeneracy.

C  Input:

C  X_A is a real variable that indicates the first coordinate of the
C  first point in the line.

C  Y_A is a real variable that indicates the second coordinate of the
C  first point in the line.

C  Z_A is a real variable that indicates the third coordinate of the
C  first point in the line.

C  X_B is a real variable that indicates the first coordinate of the
C  second point in the line.

C  Y_B is a real variable that indicates the second coordinate of the
C  second point in the line.

C  Z_B is a real variable that indicates the third coordinate of the
C  second point in the line.

C  X_D is a real variable that indicates the first coordinate of the
C  reference point.  The reference point should not be on the line.

C  Y_D is a real variable that indicates the second coordinate of the
C  reference point.  The reference point should not be on the line.

C  Z_D is a real variable that indicates the third coordinate of the
C  reference point.  The reference point should not be on the line.

C  X_F is a real variable that indicates the first coordinate of the
C  test point.

C  Y_F is a real variable that indicates the second coordinate of the
C  test point.

C  Z_F is a real variable that indicates the third coordinate of the
C  test point.

      IMPLICIT NONE
      REAL X_A,Y_A,Z_A,X_B,Y_B,Z_B,X_D,Y_D,Z_D,X_F,Y_F,Z_F,BX,BY,BZ,DX,
     $DY,DZ,ALPHA,EX,EY,EZ,FX,FY,FZ,Q

C  Calculate vector components from the first point on the line.
      BX=X_B-X_A
      BY=Y_B-Y_A
      BZ=Z_B-Z_A
      DX=X_D-X_A
      DY=Y_D-Y_A
      DZ=Z_D-Z_A
      FX=X_F-X_A
      FY=Y_F-Y_A
      FZ=Z_F-Z_A

C  Check whether the test point is on the line.
      IF(BY*FZ.EQ.BZ*FY.AND.
     $BX*FZ.EQ.BZ*FX.AND.
     $BX*FY.EQ.BY*FX)THEN
         SAME_SIDE_OF_LINE=0
      ELSE

C  Calculate the components of a vector along the line to the point
C  closest to the reference point.
         ALPHA=(DX*BX+DY*BY+DZ*BZ)/(BX**2+BY**2+BZ**2)
         EX=ALPHA*BX
         EY=ALPHA*BY
         EZ=ALPHA*BZ

C  Check whether the reference point is on the line.  If it is, the
C  result is indeterminate.
         IF(DX.EQ.EX.AND.
     $   DY.EQ.EY.AND.
     $   DZ.EQ.EZ)THEN
            SAME_SIDE_OF_LINE=-2
         ELSE

C  Compare vectors to the test point and reference point.
            Q=(DX-EX)*(FX-EX)+(DY-EY)*(FY-EY)+(DZ-EZ)*(FZ-EZ)
            IF(Q.GT.0.)THEN
               SAME_SIDE_OF_LINE=1
            ELSEIF(Q.LT.0.)THEN
               SAME_SIDE_OF_LINE=-1
            ELSE
               SAME_SIDE_OF_LINE=0
            ENDIF
         ENDIF
      ENDIF

C  Done.
      RETURN
      END
