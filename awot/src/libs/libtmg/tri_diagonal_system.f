      SUBROUTINE TRI_DIAGONAL_SYSTEM(A,A_DIM_1,N,B,SOLUTION,SUCCESS)

C  Thomas Matejka NOAA/NSSL 22 April 1994

C  This subroutine solves a set of simultaneous tri-diagonal linear
C  equations.

C  Input:

C  A is a two-dimensional real array.  A(I,J) specifies the coefficient
C  for the Jth variable in the Ith equation.  A must be tridiagonal and
C  of rank N.  A(1,1) must not be zero.

C  A_DIM_1 is an integer variable that specifies the first dimension of
C  A in the calling program.

C  N is an integer variable that specifies the number of equations,
C  which is also the number of variables.

C  B is a one-dimensional real array.  B(I) specifies the right-hand
C  side of the Ith equation.

C  Output:

C  SOLUTION is a one-dimensional real array.  SOLUTION(I) returns the
C  solution for the Ith variable.

C  SUCCESS is a logical variable that returns .TRUE. if and only if a
C  solution was found.

      IMPLICIT NONE
      LOGICAL SUCCESS
      INTEGER I,A_DIM_1,N
      REAL TERM
      REAL B(N),SOLUTION(N)
      REAL A_RR(N)
      REAL B_RR(N)
      REAL A(A_DIM_1,N)

C  Calculate the terms on the diagonal above the main diagonal and the
C  right hand side of the augmented matrix row reduced to upper diagonal
C  form.
      IF(A(1,1).EQ.0.)THEN
         SUCCESS=.FALSE.
         RETURN
      ENDIF
      A_RR(1)=A(1,2)/A(1,1)
      B_RR(1)=B(1)/A(1,1)
      DO 1 I=2,N-1
         TERM=A(I,I)-A_RR(I-1)*A(I,I-1)
         IF(TERM.EQ.0.)THEN
            SUCCESS=.FALSE.
            RETURN
         ENDIF
         A_RR(I)=A(I,I+1)/TERM
         B_RR(I)=(B(I)-B_RR(I-1)*A(I,I-1))/TERM
1     CONTINUE
      TERM=A(N,N)-A_RR(N-1)*A(N,N-1)
      IF(TERM.EQ.0.)THEN
         SUCCESS=.FALSE.
         RETURN
      ENDIF
      B_RR(N)=(B(N)-B_RR(N-1)*A(N,N-1))/TERM

C  Calculate the solution.
      SOLUTION(N)=B_RR(N)
      DO 3 I=N-1,1,-1
         SOLUTION(I)=B_RR(I)-A_RR(I)*SOLUTION(I+1)
3     CONTINUE

C  Done.
      SUCCESS=.TRUE.
      RETURN
      END
