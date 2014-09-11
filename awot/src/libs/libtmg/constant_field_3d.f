      SUBROUTINE CONSTANT_FIELD_3D(CONSTANT,MAXX,MAXY,MAXZ,NX,NY,NZ,A)

C  Thomas Matejka NOAA/NSSL 30 March 1994

C  This subroutine fills a three-dimensional data field with a
C  designated constant value at each point.

C  Input:

C  CONSTANT is a real variable that specifies the constant value to be
C  assigned to each point of A.

C  MAXX is an integer variable that specifies the first dimension of A
C  in the calling program.

C  MAXY is an integer variable that specifies the second dimension of A
C  in the calling program.

C  MAXZ is an integer variable that specifies the third dimension of A
C  in the calling program.

C  NX is an integer variable that specifies the number of grid points in
C  the first dimension for A.

C  NY is an integer variable that specifies the number of grid points in
C  the second dimension for A.

C  NZ is an integer variable that specifies the number of grid points in
C  the third dimension for A.

C  Output:

C  A is a three-dimensional real array.  A(I,J,K) returns the value, at
C  the Ith grid point in the first dimension, the Jth grid point in the
C  second dimension, and the Kth grid point in the third dimension, of
C  CONSTANT.

      IMPLICIT NONE
      INTEGER MAXX,MAXY,MAXZ,NX,NY,NZ,IX,IY,IZ
      REAL CONSTANT
      REAL A(MAXX,MAXY,MAXZ)

C  Loop through all the points and assign the constant value.
      DO 1 IZ=1,NZ
         DO 2 IY=1,NY
            DO 3 IX=1,NX
               A(IX,IY,IZ)=CONSTANT
3           CONTINUE
2        CONTINUE
1     CONTINUE

C  Done.
      RETURN
      END
