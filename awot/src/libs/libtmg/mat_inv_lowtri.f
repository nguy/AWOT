      SUBROUTINE MAT_INV_LOWTRI(A,A_DIM_1,N,A_INV,A_INV_DIM_1,
     $CHECK_SINGULAR,SINGULAR_THRESHOLD,SUCCESS)

C  Thomas Matejka NOAA/NSSL 13 September 1995

C  This subroutine calculates the inverse of a square, lower-triangular
C  matrix.

C  Input:

C  A is a two-dimensional real array.  A(I,J) specifies the element in
C  the Ith row and Jth column of the square, lower-triangular matrix
C  whose inverse is sought.

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
C  inverse of A.  A_INV has N rows and N columns and is
C  lower-triangular.

C  SUCCESS is a logical variable that returns .TRUE. if and only if an
C  inverse could be found.

      IMPLICIT NONE
      LOGICAL MAT_SINGULAR
      LOGICAL SUCCESS,CHECK_SINGULAR
      INTEGER I,J,K,A_DIM_1,A_INV_DIM_1,N
      REAL SUM,SINGULAR_THRESHOLD
      REAL A(A_DIM_1,N)
      REAL A_INV(A_INV_DIM_1,N)

C  Check whether the matrix is nearly singular if requested.
      IF(CHECK_SINGULAR)THEN
         IF(MAT_SINGULAR(A,A_DIM_1,N,SINGULAR_THRESHOLD))THEN
            SUCCESS=.FALSE.
            RETURN
         ENDIF
      ENDIF

C  Loop through the columns of the inverse matrix.
      DO 1 J=1,N

C  Set the elements above the main diagonal to zero.
         IF(J.GT.1)THEN
            DO 2 I=1,J-1
               A_INV(I,J)=0.
2           CONTINUE
         ENDIF

C  Calculate the element on the main diagonal.
         A_INV(J,J)=1./A(J,J)

C  Calculate the elements below the main diagonal.
         IF(J.LT.N)THEN
            DO 3 I=J+1,N
               SUM=0.
               DO 4 K=J,I-1
                  SUM=SUM+A(I,K)*A_INV(K,J)
4              CONTINUE
               A_INV(I,J)=-SUM/A(I,I)
3           CONTINUE
         ENDIF
1     CONTINUE

C  Done.
      SUCCESS=.TRUE.
      RETURN
      END
