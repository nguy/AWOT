      SUBROUTINE MAT_MULT(A,A_DIM_1,B,B_DIM_1,M,P,N,C,C_DIM_1)

C  Thomas Matejka NOAA/NSSL 5 January 1994

C  This subroutine multiples two matrices.

C  Input:

C  A is a two-dimensional real array.  A(I,J) specifies the element in
C  the Ith row and Jth column of the first matrix to be multiplied.

C  A_DIM_1 is an integer variable that specifies the first dimension of
C  A in the calling program.

C  B is a two-dimensional real array.  B(I,J) specifies the element in
C  the Ith row and Jth column of the second matrix to be multiplied.

C  B_DIM_1 is an integer variable that specifies the first dimension of
C  B in the calling program.

C  M is an integer variable that specifies the number of rows in A.

C  P is an integer variable that specifies the number of columns in A
C  and the number of rows in B.

C  N is an integer variable that specifies the number of columns in B.

C  C_DIM_1 is an integer variable that specifies the first dimension of
C  C in the calling program.

C  Output:

C  C is a two-dimensional real array.  C(I,J) returns the element in the
C  Ith row and Jth column of the product matrix.  C has M rows and N
C  columns.

      IMPLICIT NONE
      INTEGER A_DIM_1,M,P,B_DIM_1,N,C_DIM_1,I,J,K
      REAL SUM
      REAL A(A_DIM_1,P)
      REAL B(B_DIM_1,N)
      REAL C(C_DIM_1,N)

      DO 1 J=1,N
         DO 2 I=1,M
            SUM=0.
            DO 3 K=1,P
               SUM=SUM+A(I,K)*B(K,J)
3           CONTINUE
            C(I,J)=SUM
2        CONTINUE
1     CONTINUE
      RETURN
      END
