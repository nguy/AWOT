      SUBROUTINE MAT_CHOLESKY_FACTOR(A,A_DIM_1,N,L,L_DIM_1,SUCCESS)

C  Thomas Matejka NOAA/NSSL 5 January 1994

C  This subroutine calculates the Cholesky factor of a symmetric,
C  positive-definite matrix.  The Cholesky factor is a lower-triangular
C  matrix that, when postmultiplied by its transpose, gives the
C  symmetric matrix.

C  Input:

C  A is a two-dimensional real array.  A(I,J) specifies the element in
C  the Ith row and Jth column of the symmetric, positive-definite matrix
C  to be factored.

C  A_DIM_1 is an integer variable that specifies the first dimension of
C  A in the calling program.

C  N is an integer variable that specifies the number of rows and the
C  number of columns in A.

C  L_DIM_1 is an integer variable that specifies the first dimension of
C  L in the calling program.

C  Output:

C  L is a two-dimensional real array.  L(I,J) specifies the element in
C  the Ith row and Jth column of the Cholesky factor of A.  L has N rows
C  and N columns and is lower triangular.

C  SUCCESS is a logical variable that returns .TRUE. if and only if the
C  Cholesky factor was able to be found.

      IMPLICIT NONE
      LOGICAL SUCCESS
      INTEGER A_DIM_1,N,L_DIM_1,I,J,K
      REAL SUM,B
      REAL A(A_DIM_1,N)
      REAL L(L_DIM_1,N)

C  Loop through all the elements of the Cholesky factor.
      DO 1 I=1,N
         DO 2 J=1,N

C  Calculate an element on the main diagonal.
            IF(J.EQ.I)THEN
               SUM=0.
               IF(J.GE.2)THEN
                  DO 3 K=1,J-1
                     SUM=SUM+L(J,K)**2
3                 CONTINUE
               ENDIF
               B=A(J,J)-SUM
               IF(B.LE.0.)THEN
                  SUCCESS=.FALSE.
                  RETURN
               ENDIF
               L(J,J)=SQRT(B)

C  Calculate an element below the main diagonal.
            ELSEIF(J.LT.I)THEN
               SUM=0.
               IF(J.GE.2)THEN
                  DO 4 K=1,J-1
                     SUM=SUM+L(I,K)*L(J,K)
4                 CONTINUE
               ENDIF
               L(I,J)=(A(I,J)-SUM)/L(J,J)

C  Set an element above the main diagonal to zero.
            ELSE
               L(I,J)=0.
            ENDIF
2        CONTINUE
1     CONTINUE

C  Done.
      SUCCESS=.TRUE.
      RETURN
      END
