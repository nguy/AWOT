      SUBROUTINE GRIDBOX_1D(XMIN,DELX,NX,X,
     $IX_FACE,X_FACE)

C  Thomas Matejka NOAA/NSSL 28 April 1998
C=======================================================================
C  This subroutine finds consecutive points in a one-dimensional grid
C  that enclose a specified coordinate.  When the specified coordinate
C  lies beyond the edge of the grid, the consecutive grid points are
C  those at the closer edge of the grid.  When the grid consists of only
C  one point, the consecutive grid points are both this point.

C  Input:

C  XMIN (real) specifies the minimum grid coordinate.

C  DELX (real) specifies the grid point spacing.

C  NX (integer) specifies the number of grid points.

C  X (real) specifies the coordinate for which to find the enclosing
C  grid points.

C  Output:

C  IX_FACE(1) and IX_FACE(2) (integer) return the smaller and larger
C  consecutive grid point numbers that enclose X.

C  X_FACE(1) and X_FACE(2) (real) return the smaller and larger
C  consecutive grid coordinates that enclose X.
C=======================================================================
      IMPLICIT NONE
      INTEGER::NX,I
      INTEGER,DIMENSION(2)::IX_FACE
      REAL::X,XMIN,DELX,XMAX
      REAL,DIMENSION(2)::X_FACE

      IF(NX.GT.1)THEN
         XMAX=XMIN+DELX*FLOAT(NX-1)
         IF(X.LE.XMIN)THEN
            IX_FACE(1)=1
            IX_FACE(2)=2
         ELSEIF(X.GE.XMAX)THEN
            IX_FACE(1)=NX-1
            IX_FACE(2)=NX
         ELSE
            IX_FACE(1)=IFIX((X-XMIN)/DELX)+1
            IX_FACE(2)=IX_FACE(1)+1
         ENDIF
      ELSE
         IX_FACE(1)=1
         IX_FACE(2)=1
      ENDIF
      DO I=1,2
         X_FACE(I)=XMIN+FLOAT(IX_FACE(I)-1)*DELX
      ENDDO

      END SUBROUTINE GRIDBOX_1D
