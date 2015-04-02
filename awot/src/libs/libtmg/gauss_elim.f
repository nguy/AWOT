      SUBROUTINE GAUSS_ELIM(A,A_DIM_1,N,B,B_DIM_1,N_PROBLEMS,
     $CHECK_SINGULAR,SINGULAR_THRESHOLD,WRITE_MODE,SOLUTION,
     $SOLUTION_DIM_1,SUCCESS)

C  Thomas Matejka NOAA/NSSL 13 September 1995

C  This subroutine solves a set of simultaneous linear equations by
C  Gaussian elimination to upper-triangular form and back-substitution.

C  The set of equations is of the form Ax = b, where A is a square
C  matrix of coefficients, x is a column vector of the variables, and b
C  is the right-hand side column vector.

C  The number of equations must equal the number of variables (i.e., no
C  row or column of A can be all zero).

C  This subroutine can simultaneously solve several problems having the
C  same coefficients but different right-hand sides.

C  Input:

C  A is a two-dimensional real array.  A(I,J) specifies the coefficient
C  for the Jth variable in the Ith equation.

C  A_DIM_1 is an integer variable that specifies the first dimension of
C  A in the calling program.

C  N is an integer variable that specifies the number of equations,
C  which is also the number of variables.

C  B is a two-dimensional real array.  B(I,J) specifies the right-hand
C  side of the Ith equation for the Jth problem.

C  B_DIM_1 is an integer variable that specifies the first dimension of
C  B in the calling program.

C  N_PROBLEMS is an integer variable that specifies the number of
C  problems.

C  CHECK_SINGULAR is a logical variable.  Matrices to be inverted are
C  checked for near singularity if and only if CHECK_SINGULAR is .TRUE.

C  SINGULAR_THRESHOLD is a real variable that controls the definition of
C  nearly singular.  SINGULAR_THRESHOLD should be greater than or equal
C  to 0..  The smaller SINGULAR_THRESHOLD is, the more close to exactly
C  singular matrices must be to be considered nearly singular.  It is
C  relevant only when CHECK_SINGULAR is .TRUE..

C  WRITE_MODE is an integer variable that specifies whether matrices are
C  written at various steps of the row reduction.  If WRITE_MODE is 0,
C  nothing is written.  If WRITE_MODE is 1, the initial and row-reduced
C  matrices are written.  If WRITE_MODE is 2, the initial, final, and
C  all intermediate matrices during the row reduction are written.

C  SOLUTION_DIM_1 is an integer variable that specifies the first
C  dimension of SOLUTION in the calling program.

C  Output:

C  SOLUTION is a two-dimensional real array.  The solution for the Ith
C  variable in the Jth problem is returned in SOLUTION(I,J).

C  SUCCESS is a logical variable that returns .TRUE. if and only if a
C  solution was found.

      IMPLICIT NONE
      LOGICAL MAT_SINGULAR
      LOGICAL SUCCESS,CHECK_SINGULAR
      INTEGER I,J,K,A_DIM_1,B_DIM_1,N,N_PROBLEMS,WRITE_MODE,
     $SOLUTION_DIM_1
      REAL SINGULAR_THRESHOLD
      REAL A(A_DIM_1,N)
      REAL B(B_DIM_1,N_PROBLEMS)
      REAL SOLUTION(SOLUTION_DIM_1,N_PROBLEMS)
      DOUBLEPRECISION DP_MATRIX(N,N+N_PROBLEMS)

C  Check whether the matrix is nearly singular if requested.
      IF(CHECK_SINGULAR)THEN
         IF(MAT_SINGULAR(A,A_DIM_1,N,SINGULAR_THRESHOLD))THEN
            SUCCESS=.FALSE.
            RETURN
         ENDIF
      ENDIF

C  Load A into the left side of the double precision matrix.
      DO 1 J=1,N
         DO 2 I=1,N
            DP_MATRIX(I,J)=DBLE(A(I,J))
2        CONTINUE
1     CONTINUE

C  Load B into the right side of the double precision matrix.
      DO 3 J=1,N_PROBLEMS
         DO 4 I=1,N
            DP_MATRIX(I,N+J)=DBLE(B(I,J))
4        CONTINUE
3     CONTINUE

C  Row reduce the double precision matrix to upper triangular form.
      CALL UPTRI_PP(DP_MATRIX,N,N,N+N_PROBLEMS,WRITE_MODE)

C  Check whether the row reduction was successful.
      DO 5 J=1,N
         DO 6 I=J,N
            IF(I.EQ.J)THEN
               IF(DP_MATRIX(I,J).NE.1.D0)THEN
                  SUCCESS=.FALSE.
                  RETURN
               ENDIF
            ELSE
               IF(DP_MATRIX(I,J).NE.0.D0)THEN
                  SUCCESS=.FALSE.
                  RETURN
               ENDIF
            ENDIF
6        CONTINUE
5     CONTINUE

C  Back substitute the right side of the double precision matrix to
C  obtain the solutions.
      IF(N.GE.2)THEN
         DO 7 J=N+1,N+N_PROBLEMS
            DO 8 I=N-1,1,-1
               DO 9 K=I+1,N
                  DP_MATRIX(I,J)=DP_MATRIX(I,J)-
     $            DP_MATRIX(I,K)*DP_MATRIX(K,J)
9              CONTINUE
8           CONTINUE
7        CONTINUE
      ENDIF

C  Load the solution into SOLUTION.
      DO 10 J=1,N_PROBLEMS
         DO 11 I=1,N
            SOLUTION(I,J)=SNGL(DP_MATRIX(I,N+J))
11       CONTINUE
10    CONTINUE

C  Done.
      SUCCESS=.TRUE.
      RETURN
      END
