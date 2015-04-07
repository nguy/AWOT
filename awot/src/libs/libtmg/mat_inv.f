      SUBROUTINE MAT_INV(A,A_DIM_1,N,A_INV,A_INV_DIM_1,CHECK_SINGULAR,
     $SINGULAR_THRESHOLD,WRITE_MODE,SUCCESS)

C  Thomas Matejka NOAA/NSSL 13 September 1995

C  This subroutine calculates the inverse of a square matrix.

C  Input:

C  A is a two-dimensional real array.  A(I,J) specifies the element in
C  the Ith row and Jth column of the square matrix whose inverse is
C  sought.

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

C  WRITE_MODE is an integer variable that specifies whether matrices are
C  written at various steps of the row reduction.  If WRITE_MODE is 0,
C  nothing is written.  If WRITE_MODE is 1, the initial and row-reduced
C  matrices are written.  If WRITE_MODE is 2, the initial, final, and
C  all intermediate matrices during the row reduction are written.

C  Output:

C  A_INV is a two-dimensional real array.  A_INV(I,J) returns the
C  element in the Ith row and Jth column of the matrix that is the
C  inverse of A.  A_INV has N rows and N columns.

C  SUCCESS is a logical variable that returns .TRUE. if and only if an
C  inverse could be found.

      IMPLICIT NONE
      LOGICAL MAT_SINGULAR
      LOGICAL SUCCESS,CHECK_SINGULAR
      INTEGER I,J,K,A_DIM_1,A_INV_DIM_1,N,WRITE_MODE
      REAL SINGULAR_THRESHOLD
      REAL A(A_DIM_1,N)
      REAL A_INV(A_INV_DIM_1,N)
      DOUBLEPRECISION DP_MATRIX(N,2*N)

C  Check whether the matrix is nearly singular if requested.
      IF(CHECK_SINGULAR)THEN
         IF(MAT_SINGULAR(A,A_DIM_1,N,SINGULAR_THRESHOLD))THEN
            SUCCESS=.FALSE.
            RETURN
         ENDIF
      ENDIF

C  Load A into the left half of the augmented matrix.
      DO 1 J=1,N
         DO 2 I=1,N
            DP_MATRIX(I,J)=DBLE(A(I,J))
2        CONTINUE
1     CONTINUE

C  Load the identity matrix into the right half of the augmented matrix.
      DO 3 J=N+1,2*N
         DO 4 I=1,N
            IF(I.EQ.J-N)THEN
               DP_MATRIX(I,J)=1.D0
            ELSE
               DP_MATRIX(I,J)=0.D0
            ENDIF
4        CONTINUE
3     CONTINUE

C  Row reduce the augmented matrix to upper triangular form.
      CALL UPTRI_PP(DP_MATRIX,N,N,2*N,WRITE_MODE)

C  Check whether an inverse was found.
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

C  Back substitute the right half of the augmented matrix to obtain the
C  inverse.
      IF(N.GE.2)THEN
         DO 7 J=N+1,2*N
            DO 8 I=N-1,1,-1
               DO 9 K=I+1,N
                  DP_MATRIX(I,J)=DP_MATRIX(I,J)-
     $            DP_MATRIX(I,K)*DP_MATRIX(K,J)
9              CONTINUE
8           CONTINUE
7        CONTINUE
      ENDIF

C  Load the right half of the augmented matrix into A_INV.
      DO 10 J=1,N
         DO 11 I=1,N
            A_INV(I,J)=SNGL(DP_MATRIX(I,J+N))
11       CONTINUE
10    CONTINUE

C  Done.
      SUCCESS=.TRUE.
      RETURN
      END
