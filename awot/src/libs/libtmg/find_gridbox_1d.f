      SUBROUTINE FIND_GRIDBOX_1D(X,XMIN,DELX,NX,IX_FACE,X_FACE)

C  Thomas Matejka NOAA/NSSL 5 October 1993

C  This subroutine finds the coordinates of the faces of the
C  one-dimensional grid box surrounding the point X.

C  The grid consists of NX points.  XMIN is the minimum grid coordinate.
C  The grid points are spaced DELX.

C  In IX_FACE(1) and IX_FACE(2) are returned the smaller and larger grid
C  point numbers of the box.  In X_FACE(1) and X_FACE(2) are returned
C  the smaller and larger coordinates of the box.

C  When the point coincides with the grid, the box coincides with the
C  grid and has zero width.

C  When the point lies beyond the edge of the grid, the box is
C  positioned at the edge of the grid and has zero width.

      IMPLICIT NONE
      INTEGER NX,I
      INTEGER IX_FACE(2)
      REAL X,XMIN,DELX,XMAX
      REAL X_FACE(2)

C  Find the grid point numbers.
      IF(NX.GT.1)THEN
         XMAX=XMIN+DELX*FLOAT(NX-1)
         IF(X.LE.XMIN)THEN
            IX_FACE(1)=1
            IX_FACE(2)=1
         ELSEIF(X.GE.XMAX)THEN
            IX_FACE(1)=NX
            IX_FACE(2)=NX
         ELSE
            IF(AMOD(X-XMIN,DELX).EQ.0.)THEN
               IX_FACE(1)=IFIX((X-XMIN)/DELX)+1
               IX_FACE(2)=IX_FACE(1)
            ELSE
               IX_FACE(1)=IFIX((X-XMIN)/DELX)+1
               IX_FACE(2)=IX_FACE(1)+1
            ENDIF
         ENDIF
      ELSE
         IX_FACE(1)=1
         IX_FACE(2)=1
      ENDIF

C  Find the coordinates.
      DO 1 I=1,2
         X_FACE(I)=XMIN+FLOAT(IX_FACE(I)-1)*DELX
1     CONTINUE

C  Done.
      RETURN
      END
