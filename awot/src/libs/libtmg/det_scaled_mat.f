      REAL FUNCTION DET_SCALED_MAT(A,A_DIM_1,N)

C  Thomas Matejka NOAA/NSSL 7 September 1995

C  This real function scales a matrix so that the largest magnitude in a
C  row or column is 1. and then calculates its determinant.

C  Input:

C  A is a two-dimensional real array.  A(I,J) specifies the element in
C  the Ith row and Jth column of the square matrix whose scaled
C  determinant is sought.

C  A_DIM_1 is an integer variable that specifies the first dimension of
C  A in the calling program.

C  N is an integer variable that specifies the number of rows and the
C  number of columns in A.

C  Output:

C  The real function returns the determinant of the scaled matrix.

      IMPLICIT NONE
      DOUBLEPRECISION DETERMINANT
      INTEGER A_DIM_1,N,I,J
      REAL A(A_DIM_1,N)
      DOUBLEPRECISION X
      DOUBLEPRECISION B(N,N)

C  Copy A into B.
      DO J=1,N
         DO I=1,N
            B(I,J)=DBLE(A(I,J))
         ENDDO
      ENDDO

C  Scale the rows of the matrix so that the largest magnitude in a row
C  is 1..
      DO I=1,N
         X=DABS(B(I,1))
         IF(N.GE.2)THEN
            DO J=2,N
               X=DMAX1(X,DABS(B(I,J)))
            ENDDO
         ENDIF
         IF(X.GT.0.D0)THEN
            DO J=1,N
               B(I,J)=B(I,J)/X
            ENDDO
         ELSE
            DET_SCALED_MAT=0.
            RETURN
         ENDIF
      ENDDO

C  Scale the columns of the matrix so that the largest magnitude in a
C  column is 1..
      DO J=1,N
         X=DABS(B(1,J))
         IF(N.GE.2)THEN
            DO I=2,N
               X=DMAX1(X,DABS(B(I,J)))
            ENDDO
         ENDIF
         IF(X.GT.0.D0)THEN
            DO I=1,N
               B(I,J)=B(I,J)/X
            ENDDO
         ELSE
            DET_SCALED_MAT=0.
            RETURN
         ENDIF
      ENDDO

C  Calculate the determinant of the scaled matrix.
      DET_SCALED_MAT=SNGL(DETERMINANT(B,N,N))

C  Done.
      RETURN
      END
