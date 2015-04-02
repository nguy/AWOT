      SUBROUTINE UPTRI_PP(A,A_DIM_1,NUM_ROWS_A,NUM_COLS_A,WRITE_MODE)

C  Thomas Matejka NOAA/NSSL 2 September 1993

C  This subroutine row-reduces a matrix to upper-triangular form using
C  partial pivoting.

C  Input:

C  A is a two-dimensional double-precision array.  A(I,J) specifies the
C  element in the Ith row and Jth column of the matrix to be
C  row-reduced.

C  A_DIM_1 is an integer variable that specifies the first dimension of
C  A in the calling program.

C  NUM_ROWS_A is an integer variable that specifies the number of rows
C  in A.

C  NUM_COLS_A is an integer variable that specifies the number of
C  columns in A.

C  WRITE_MODE is an integer variable that specifies whether matrices are
C  written at various steps of the row reduction.  If WRITE_MODE is 0,
C  nothing is written.  If WRITE_MODE is 1, the initial and row-reduced
C  matrices are written.  If WRITE_MODE is 2, the initial, final, and
C  all intermediate matrices during the row reduction are written.

C  LU is an integer variable that specifies the logical device number to
C  write to.

C  Output:

C  The row-reduced form replaces the original data in A.

      IMPLICIT NONE
      INTEGER WRITE_MODE,A_DIM_1,NUM_ROWS_A,NUM_COLS_A,PIVOT_ROW,
     $PIVOT_COL,IROW,ICOL,IROW_STORE
      DOUBLEPRECISION W,X,Y,Z,STORE
      DOUBLEPRECISION A(A_DIM_1,NUM_COLS_A)

C  Write the original matrix if requested.
      IF(WRITE_MODE.GE.1)THEN
         CALL MAT_WRITE(A,A_DIM_1,NUM_ROWS_A,NUM_COLS_A,1)
      ENDIF

C  Initialize the pivot row and the pivot column.
      PIVOT_ROW=1
      PIVOT_COL=1

C  Check whether the pivot column is all zero.
1     CONTINUE
      DO 2 IROW=PIVOT_ROW,NUM_ROWS_A
         IF(A(IROW,PIVOT_COL).NE.0.D0)THEN
            GOTO 3
         ENDIF
2     CONTINUE

C  The pivot column is all zero.  Advance the pivot column.
      IF(PIVOT_COL.LT.NUM_COLS_A)THEN
         PIVOT_COL=PIVOT_COL+1
         GOTO 1
      ELSE
         GOTO 10
      ENDIF

C  The pivot column is not all zero.  Interchange rows so that the term
C  in the pivot column and pivot row divided by the largest term in the
C  row is a maximum.
3     CONTINUE
      X=0.D0
      DO 4 IROW=PIVOT_ROW,NUM_ROWS_A
         W=0.D0
         DO 5 ICOL=PIVOT_COL,NUM_COLS_A
            IF(DABS(A(IROW,ICOL)).GT.W)THEN
               W=DABS(A(IROW,ICOL))
            ENDIF
5        CONTINUE
         IF(W.NE.0.D0)THEN
            IF(DABS(A(IROW,PIVOT_COL))/W.GT.X)THEN
               X=DABS(A(IROW,PIVOT_COL)/W)
               IROW_STORE=IROW
            ENDIF
         ENDIF
4     CONTINUE
      IF(IROW_STORE.NE.PIVOT_ROW)THEN
         DO 6 ICOL=1,NUM_COLS_A
            STORE=A(PIVOT_ROW,ICOL)
            A(PIVOT_ROW,ICOL)=A(IROW_STORE,ICOL)
            A(IROW_STORE,ICOL)=STORE
6        CONTINUE
         IF(WRITE_MODE.GE.2)THEN
            CALL MAT_WRITE(A,A_DIM_1,NUM_ROWS_A,NUM_COLS_A,2)
         ENDIF
      ENDIF

C  Divide the terms in the pivot row by the term in the pivot row and
C  column.
      Y=A(PIVOT_ROW,PIVOT_COL)
      A(PIVOT_ROW,PIVOT_COL)=1.D0
      IF(PIVOT_COL+1.LE.NUM_COLS_A)THEN
         DO 7 ICOL=PIVOT_COL+1,NUM_COLS_A
            A(PIVOT_ROW,ICOL)=A(PIVOT_ROW,ICOL)/Y
7        CONTINUE
      ENDIF
      IF(WRITE_MODE.GE.2)THEN
         CALL MAT_WRITE(A,A_DIM_1,NUM_ROWS_A,NUM_COLS_A,2)
      ENDIF

C  Eliminate the terms in the pivot column in rows below the pivot row.
      IF(PIVOT_ROW+1.LE.NUM_ROWS_A)THEN
         DO 8 IROW=PIVOT_ROW+1,NUM_ROWS_A
            Z=A(IROW,PIVOT_COL)
            A(IROW,PIVOT_COL)=0.D0
            IF(PIVOT_COL+1.LE.NUM_COLS_A)THEN
               DO 9 ICOL=PIVOT_COL+1,NUM_COLS_A
                  A(IROW,ICOL)=A(IROW,ICOL)-A(PIVOT_ROW,ICOL)*Z
9              CONTINUE
            ENDIF
8        CONTINUE
         IF(WRITE_MODE.GE.2)THEN
            CALL MAT_WRITE(A,A_DIM_1,NUM_ROWS_A,NUM_COLS_A,2)
         ENDIF
      ENDIF
 
C  Row reduction for the current pivot row and pivot column is complete.
C  Advance the pivot row and pivot column.
      IF(PIVOT_ROW.LT.NUM_ROWS_A.AND.PIVOT_COL.LT.NUM_COLS_A)THEN
         PIVOT_ROW=PIVOT_ROW+1
         PIVOT_COL=PIVOT_COL+1
         GOTO 1
      ENDIF

C  The matrix is row-reduced.  Write the row-reduced matrix if
C  requested.
10    CONTINUE
      IF(WRITE_MODE.EQ.1)THEN
         CALL MAT_WRITE(A,A_DIM_1,NUM_ROWS_A,NUM_COLS_A,3)
      ENDIF

C  Done.
      RETURN
      END
