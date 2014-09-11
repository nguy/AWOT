      SUBROUTINE MAT_WRITE(A,A_DIM_1,M,N,CODE)

C  Thomas Matejka NOAA/NSSL 5 January 1994

C  This subroutine writes a matrix.

C  Input:

C  A is a two-dimensional double-precision array.  A(I,J) specifies the
C  element in the Ith row and Jth column of the matrix to be written.

C  A_DIM_1 is an integer variable that specifies the first dimension of
C  A in the calling program.

C  M is an integer variable that specifies the number of rows in A.

C  N is an integer variable that specifies the number of columns in A.

C  CODE is an integer variable that controls the matrix label.  If CODE
C  is 0, there is no label.  If CODE is 1, the label is "Original
C  matrix:".  If CODE is 2, the label is "Matrix after row-reduction
C  operation:".  If CODE is 3, the label is "Row-reduced matrix:".

      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      INTEGER A_DIM_1,M,N,CODE,I,J
      DOUBLEPRECISION A(A_DIM_1,N)

C  Write the matrix label.
      IF(CODE.EQ.0)THEN
         CONTINUE
      ELSEIF(CODE.EQ.1)THEN
         WRITE(TMMLIB_MESSAGE_UNIT,*)' Original matrix:'
         WRITE(TMMLIB_MESSAGE_UNIT,*)
      ELSEIF(CODE.EQ.2)THEN
         WRITE(TMMLIB_MESSAGE_UNIT,*)' Matrix after row-reduction ',
     $   'operation:'
         WRITE(TMMLIB_MESSAGE_UNIT,*)
      ELSEIF(CODE.EQ.3)THEN
         WRITE(TMMLIB_MESSAGE_UNIT,*)' Row-reduced matrix:'
         WRITE(TMMLIB_MESSAGE_UNIT,*)
      ENDIF

C  Write the matrix.
      DO 1 I=1,M
         WRITE(TMMLIB_MESSAGE_UNIT,2)(A(I,J),J=1,N)
2        FORMAT(1X,8E16.8)
         WRITE(TMMLIB_MESSAGE_UNIT,*)
1     CONTINUE
      WRITE(TMMLIB_MESSAGE_UNIT,*)
      WRITE(TMMLIB_MESSAGE_UNIT,*)

C  Done.
      RETURN
      END
