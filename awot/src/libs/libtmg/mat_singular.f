      LOGICAL FUNCTION MAT_SINGULAR(A,A_DIM_1,N,SINGULAR_THRESHOLD)

C  Thomas Matejka NOAA/NSSL 13 September 1995

C  This function tests whether a square matrix is singular or nearly
C  singular.

C  Input:

C  A is a two-dimensional real array.  A(I,J) specifies the element in
C  the Ith row and Jth column of the square matrix to be tested for near
C  singularity.

C  A_DIM_1 is an integer variable that specifies the first dimension of
C  A in the calling program.

C  N is an integer variable that specifies the number of rows and the
C  number of columns in A.

C  SINGULAR_THRESHOLD is a real variable that controls the definition of
C  nearly singular.  SINGULAR_THRESHOLD should be greater than or equal
C  to 0..  The smaller SINGULAR_THRESHOLD is, the more close to exactly
C  singular matrix A must be to be considered nearly singular.

C  Output:

C  The logical function returns .TRUE. if and only if A is singular or
C  nearly singular.

      IMPLICIT NONE
      REAL DET_SCALED_MAT
      INTEGER A_DIM_1,N
      REAL SINGULAR_THRESHOLD
      REAL A(A_DIM_1,N)

C  Compare the determinant of the scaled matrix to the specified
C  threshold.
      IF(ABS(DET_SCALED_MAT(A,A_DIM_1,N)).LT.SINGULAR_THRESHOLD)
     $THEN
         MAT_SINGULAR=.TRUE.
      ELSE
         MAT_SINGULAR=.FALSE.
      ENDIF

C  Done.
      RETURN
      END
