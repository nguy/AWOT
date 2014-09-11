      SUBROUTINE GAUSS_JORDAN_ELIM_MISS_VARS(A,A_DIM_1,N,
     $B,B_DIM_1,N_PROBLEMS,CHECK_SINGULAR,SINGULAR_THRESHOLD,WRITE_MODE,
     $SOLUTION,SOLUTION_DIM_1,SOLUTION_EXISTS,SUCCESS)

C  Thomas Matejka NOAA/NSSL 13 September 1995

C  This subroutine solves a set of simultaneous linear equations by
C  Gauss-Jordan elimination to row-echelon form.

C  The set of equations is of the form Ax = b, where A is a square
C  matrix of coefficients, x is a column vector of the variables, and b
C  is the right-hand side column vector.

C  The coefficients and right-hand sides of some equations may consist
C  of all zeros.  However, the number of not all-zero equations must
C  equal the number of variables with not all-zero coefficients (i.e.,
C  the number of not all-zero rows and not all-zero columns of A must be
C  equal).

C  This subroutine can solve simultaneously several problems having the
C  same coefficients and variables but different right-hand sides.

C  Input:

C  A is a two-dimensional real array.  A(I,J) specifies the coefficient
C  for the Jth variable in the Ith equation.

C  A_DIM_1 is an integer variable that specifies the first dimension of
C  A in the calling program.

C  N is an integer variable that specifies the number of equations
C  (including all-zero equations), which is also the number of variables
C  (including variables that have only zero coefficients).

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
C  variable in the Jth problem is returned in SOLUTION(I,J).  The
C  original numbering of the variables is preserved regardless of
C  whether some variables had only zero coefficients.

C  SOLUTION_EXISTS is a one-dimensional logical array.
C  SOLUTION_EXISTS(I) is .TRUE. if and only if the Ith variable can be
C  solved for.  If SOLUTION_EXISTS(I) is .FALSE., then all the
C  coefficients for the Ith variable were zero.

C  SUCCESS is a logical variable that returns .TRUE. if and only if a
C  solution was found.

      IMPLICIT NONE
      LOGICAL SUCCESS,CHECK_SINGULAR
      LOGICAL SOLUTION_EXISTS(N)
      LOGICAL NON_ZERO_ROW(N),NON_ZERO_COL(N)
      INTEGER I,J,II,JJ,A_DIM_1,B_DIM_1,N,N_PROBLEMS,WRITE_MODE,
     $SOLUTION_DIM_1,NUM_NON_ZERO_COLS,NUM_NON_ZERO_ROWS
      REAL SINGULAR_THRESHOLD
      REAL A(A_DIM_1,N)
      REAL B(B_DIM_1,N_PROBLEMS)
      REAL A_PACKED(N,N)
      REAL B_PACKED(N,N_PROBLEMS)
      REAL SOLUTION(SOLUTION_DIM_1,N_PROBLEMS)
      REAL SOLUTION_PACKED(N,N_PROBLEMS)

C  Count and identify the number of not-all-zero columns in A.
      NUM_NON_ZERO_COLS=0
      DO 1 J=1,N
         DO 2 I=1,N
            IF(A(I,J).NE.0.)THEN
               NON_ZERO_COL(J)=.TRUE.
               NUM_NON_ZERO_COLS=NUM_NON_ZERO_COLS+1
               GOTO 3
            ENDIF
2        CONTINUE
         NON_ZERO_COL(J)=.FALSE.
3        CONTINUE
1     CONTINUE

C  Count and identify the number of not-all-zero rows in A.
      NUM_NON_ZERO_ROWS=0
      DO 4 I=1,N
         DO 5 J=1,N
            IF(A(I,J).NE.0.)THEN
               NON_ZERO_ROW(I)=.TRUE.
               NUM_NON_ZERO_ROWS=NUM_NON_ZERO_ROWS+1
               GOTO 6
            ENDIF
5        CONTINUE
         NON_ZERO_ROW(I)=.FALSE.
6        CONTINUE
4     CONTINUE

C  If the number of not-all-zero columns is greater than the number of
C  not-all-zero rows in A, the set of equations is underdetermined.
      IF(NUM_NON_ZERO_COLS.GT.NUM_NON_ZERO_ROWS)THEN
         SUCCESS=.FALSE.
         RETURN
      ENDIF

C  If the number of not-all-zero rows is greater than the number of
C  not-all-zero columns in A, the set of equations is either
C  inconsistent or redundant.
      IF(NUM_NON_ZERO_ROWS.GT.NUM_NON_ZERO_COLS)THEN
         SUCCESS=.FALSE.
         RETURN
      ENDIF

C  Check whether the corresponding row of B is all zero whenever a row
C  of A is all zero.  If it is not, the set of equations is
C  inconsistent.
      DO 7 I=1,N
         IF(.NOT.NON_ZERO_ROW(I))THEN
            DO 8 J=1,N_PROBLEMS
               IF(B(I,J).NE.0.)THEN
                  SUCCESS=.FALSE.
                  RETURN
               ENDIF
8           CONTINUE
         ENDIF
7     CONTINUE

C  Eliminate all non-zero rows and columns of A and all non-zero rows of
C  B.
      II=0
      DO 9 I=1,N
         IF(NON_ZERO_ROW(I))THEN
            II=II+1
            JJ=0
            DO 10 J=1,N
               IF(NON_ZERO_COL(J))THEN
                  JJ=JJ+1
                  A_PACKED(II,JJ)=A(I,J)
               ENDIF
10          CONTINUE
            DO 11 J=1,N_PROBLEMS
               B_PACKED(II,J)=B(I,J)
11          CONTINUE
         ENDIF
9     CONTINUE

C  Obtain the solution by Gaussian elimination.
      CALL GAUSS_JORDAN_ELIM(A_PACKED,N,NUM_NON_ZERO_COLS,B_PACKED,N,
     $N_PROBLEMS,CHECK_SINGULAR,SINGULAR_THRESHOLD,WRITE_MODE,
     $SOLUTION_PACKED,N,SUCCESS)
      IF(.NOT.SUCCESS)THEN
         RETURN
      ENDIF

C  Restore the correspondence between the original variable numbers and
C  the solutions.  Indicate whether there is a solution for each
C  variable in SOLUTION_EXISTS.
      II=0
      DO 12 I=1,N
         IF(NON_ZERO_COL(I))THEN
            II=II+1
            SOLUTION_EXISTS(I)=.TRUE.
            DO 13 J=1,N_PROBLEMS
               SOLUTION(I,J)=SOLUTION_PACKED(II,J)
13          CONTINUE
         ELSE
            SOLUTION_EXISTS(I)=.FALSE.
            DO 14 J=1,N_PROBLEMS
               SOLUTION(I,J)=0.
14          CONTINUE
         ENDIF
12    CONTINUE

C  Done.
      SUCCESS=.TRUE.
      RETURN
      END
