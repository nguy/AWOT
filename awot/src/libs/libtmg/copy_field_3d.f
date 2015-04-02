      SUBROUTINE COPY_FIELD_3D(A,MAXX_IN,MAXY_IN,MAXZ_IN,NX,NY,NZ,
     $MAXX_OUT,MAXY_OUT,MAXZ_OUT,B)

C  Thomas Matejka NOAA/NSSL 30 March 1994

C  This subroutine copies one three-dimensional data field to another
C  three-dimensional data field.

C  Input:

C  A is a three-dimensional real array.  A(I,J,K) specifies the
C  value, at the Ith grid point in the first dimension, the Jth grid
C  point in the second dimension, and the Kth grid point in the third
C  dimension, of the data field to be copied.

C  MAXX_IN is an integer variable that specifies the first dimension of
C  A in the calling program.

C  MAXY_IN is an integer variable that specifies the second dimension of
C  A in the calling program.

C  MAXZ_IN is an integer variable that specifies the third dimension of
C  A in the calling program.

C  NX is an integer variable that specifies the number of grid points in
C  the first dimension for A and B.

C  NY is an integer variable that specifies the number of grid points in
C  the second dimension for A and B.

C  NZ is an integer variable that specifies the number of grid points in
C  the third dimension for A and B.

C  MAXX_OUT is an integer variable that specifies the first dimension of
C  B in the calling program.

C  MAXY_OUT is an integer variable that specifies the second dimension
C  of B in the calling program.

C  MAXZ_OUT is an integer variable that specifies the third dimension of
C  B and the calling program.

C  Output:

C  B is a three-dimensional real array.  B(I,J,K) returns the value, at
C  the Ith grid point in the first dimension, the Jth grid point in the
C  second dimension, and the Kth grid point in the third dimension, of
C  A(I,J,K).

      IMPLICIT NONE
      INTEGER MAXX_IN,MAXY_IN,MAXZ_IN,MAXX_OUT,MAXY_OUT,MAXZ_OUT,NX,NY,
     $NZ,IX,IY,IZ
      REAL A(MAXX_IN,MAXY_IN,MAXZ_IN)
      REAL B(MAXX_OUT,MAXY_OUT,MAXZ_OUT)

C  Loop through all the points and copy the field.
      DO 1 IZ=1,NZ
         DO 2 IY=1,NY
            DO 3 IX=1,NX
               B(IX,IY,IZ)=A(IX,IY,IZ)
3           CONTINUE
2        CONTINUE
1     CONTINUE

C  Done.
      RETURN
      END
