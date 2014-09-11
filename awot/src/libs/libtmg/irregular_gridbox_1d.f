      SUBROUTINE IRREGULAR_GRIDBOX_1D(N,X_ARRAY,X,
     $IX_FACE,X_FACE)

C  Thomas Matejka NOAA/NSSL 29 April 1998
C=======================================================================
C  This subroutine finds consecutive points in an irregular
C  one-dimensional grid that enclose a specified coordinate.  When the
C  specified coordinate lies beyond the edge of the grid, the
C  consecutive grid points are those at the closer edge of the grid.
C  When the grid consists of only one point, the consecutive grid points
C  are both this point.

C  Input:

C  N (integer) specifies the number of grid points.

C  X_ARRAY (1d real array 1:N).  X_ARRAY(I) specifies the coordinate of
C  the Ith grid point.  X(I) must be monotonically increasing for
C  I=1,...,N.

C  X (real) specifies the coordinate for which to find the enclosing
C  grid points.

C  Output:

C  IX_FACE(1) and IX_FACE(2) (integer) return the smaller and larger
C  consecutive grid point numbers that enclose X.

C  X_FACE(1) and X_FACE(2) (real) return the smaller and larger
C  consecutive grid coordinates that enclose X.
C=======================================================================
      IMPLICIT NONE
      INTEGER::N,I,IX
      INTEGER,DIMENSION(2)::IX_FACE
      REAL::X
      REAL,DIMENSION(2)::X_FACE
      REAL,DIMENSION(N)::X_ARRAY

      IF(N.GT.1)THEN
         IF(X.LE.X_ARRAY(1))THEN
            IX_FACE(1)=1
            IX_FACE(2)=2
         ELSEIF(X.GE.X_ARRAY(N))THEN
            IX_FACE(1)=N-1
            IX_FACE(2)=N
         ELSE
            DO IX=2,N
               IF(X_ARRAY(IX).GT.X)THEN
                  IX_FACE(1)=IX-1
                  IX_FACE(2)=IX
                  EXIT
               ENDIF
            ENDDO
         ENDIF
      ELSE
         IX_FACE(1)=1
         IX_FACE(2)=1
      ENDIF
      DO I=1,2
         X_FACE(I)=X_ARRAY(IX_FACE(I))
      ENDDO

      END SUBROUTINE IRREGULAR_GRIDBOX_1D
