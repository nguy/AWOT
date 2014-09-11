      SUBROUTINE MAT_INV_SYM_POSDEF(A,A_DIM_1,N,A_INV,A_INV_DIM_1,
     $CHECK_SINGULAR,SINGULAR_THRESHOLD,SUCCESS)

C  Thomas Matejka NOAA/NSSL 13 September 1995

C  This subroutine calculates the inverse of a square, symmetric,
C  positive definite matrix.

C  Input:

C  A is a two-dimensional real array.  A(I,J) specifies the element in
C  the Ith row and Jth column of the square, symmetric, positive
C  definite matrix whose inverse is sought.

C  A_DIM_1 is an integer variable that specifies the first dimension of
C  A in the calling program.

C  N is an integer variable that specifies the number of rows and the
C  number of columns of A.

C  A_INV_DIM_1 is an integer variable that specifies the first dimension
C  of A_INV in the calling program.

C  CHECK_SINGULAR is a logical variable.  Matrices to be inverted are
C  checked for near singularity if and only if CHECK_SINGULAR is .TRUE.

C  SINGULAR_THRESHOLD is a real variable that controls the definition of
C  nearly singular.  SINGULAR_THRESHOLD should be greater than or equal
C  to 0..  The smaller SINGULAR_THRESHOLD is, the more close to exactly
C  singular matrices must be to be considered nearly singular.  It is
C  relevant only when CHECK_SINGULAR is .TRUE..

C  Output:

C  A_INV is a two-dimensional real array.  A_INV(I,J) returns the
C  element in the Ith row and Jth column of the matrix that is the
C  inverse of A.  A_INV has N rows and N columns and is symmetric.

C  SUCCESS is a logical variable that returns .TRUE. if and only if an
C  inverse could be found.

      IMPLICIT NONE
      LOGICAL MAT_SINGULAR
      LOGICAL SUCCESS,CHECK_SINGULAR
      INTEGER I,J,K,A_DIM_1,A_INV_DIM_1,N
      REAL SINGULAR_THRESHOLD
      REAL A(A_DIM_1,N)
      REAL A_INV(A_INV_DIM_1,N)
      REAL L(N,N),L_INV(N,N)

C  Check whether the matrix is nearly singular if requested.
      IF(CHECK_SINGULAR)THEN
         IF(MAT_SINGULAR(A,A_DIM_1,N,SINGULAR_THRESHOLD))THEN
            SUCCESS=.FALSE.
            RETURN
         ENDIF
      ENDIF

C  Calculate the Cholesky factor of A.
      CALL MAT_CHOLESKY_FACTOR(A,A_DIM_1,N,L,N,SUCCESS)
      IF(.NOT.SUCCESS)THEN
         RETURN
      ENDIF

C  Invert the lower triangular Cholesky factor.
      CALL MAT_INV_LOWTRI(L,N,N,L_INV,N,CHECK_SINGULAR,
     $SINGULAR_THRESHOLD,SUCCESS)
      IF(.NOT.SUCCESS)THEN
         RETURN
      ENDIF

C  Calculate the inverse of A by multiplying the transpose of the
C  inverse of the Cholesky factor by the inverse of the Cholesky factor.
      DO 1 I=1,N
         DO 2 J=1,I
            A_INV(I,J)=0.
            DO 3 K=I,N
               A_INV(I,J)=A_INV(I,J)+L_INV(K,I)*L_INV(K,J)
3           CONTINUE
2        CONTINUE
1     CONTINUE
      IF(N.GE.2)THEN
         DO 4 I=2,N
            DO 5 J=1,I-1
               A_INV(J,I)=A_INV(I,J)
5           CONTINUE
4        CONTINUE
      ENDIF
C  Done.
      SUCCESS=.TRUE.
      RETURN
      END
