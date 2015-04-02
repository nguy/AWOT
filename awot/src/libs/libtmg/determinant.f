      DOUBLEPRECISION recursive FUNCTION DETERMINANT(A,A_DIM_1,N) result ( value )

C  Thomas Matejka NOAA/NSSL 18 August 1994
C  This recursive function calculates the determinant of a square matrix.

C  Input:
C  A is a two-dimensional double precision array.  A(I,J) specifies the
C  element in the Ith row and Jth column of the square matrix whose
C  determinant is sought.

C  A_DIM_1 is an integer variable that specifies the first dimension of A in the calling program.
C  N is an integer variable that specifies the number of rows and the number of columns of A.

C  Output:
C  The double precision function returns the value of the determinant of A.

      IMPLICIT NONE
      INTEGER A_DIM_1,N,I,J,JJ,J_MINOR
      DOUBLEPRECISION COFACTOR
      DOUBLEPRECISION A(A_DIM_1,N)
      DOUBLEPRECISION B(N,N)

C  The determinant of a 1 x 1 matrix is the element itself.
      IF(N.EQ.1)THEN
         Value=A(1,1)

C  The determinant of a 2 x 2 or larger matrix is defined as the sum of
C  the products of the elements of the first row times their cofactors.
      ELSE
         Value=0.D0

         DO 1 J=1,N
            J_MINOR=0
            DO 2 JJ=1,N
               IF(JJ.NE.J)THEN
                  J_MINOR=J_MINOR+1
                  DO 3 I=2,N
                     B(I-1,J_MINOR)=A(I,JJ)
3                 CONTINUE
               ENDIF
2           CONTINUE

            COFACTOR=(-1.D0)**(J+1)*DETERMINANT(B,N,N-1)
            Value = Value +A(1,J)*COFACTOR

1        CONTINUE
      ENDIF

      RETURN
      END
