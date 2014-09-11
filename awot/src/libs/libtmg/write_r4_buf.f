      SUBROUTINE WRITE_R4_BUF(R4BUF,N,STRING)

C  Thomas Matejka NOAA/NSSL 9 February 1996

C  This subroutine interprets a one-dimensional array of four-byte
C  floating point numbers as a series of sets of four unformatted bytes.

C  Input:

C  R4BUF is a one-dimensional real*4 array.  R4BUF(I) is the Ith
C  four-byte floating point number to be interpreted as the Ith set of
C  four unformatted bytes.

C  N is an integer variable that specifies the number of elements in
C  R4BUF and the number of sets of four unformatted bytes in STRING.

C  Output:

C  STRING is the character string whose first 4*N unformatted bytes
C  correspond to N four-byte floating-point numbers.

      IMPLICIT NONE
      CHARACTER*4 C4
      CHARACTER*(*) STRING
      REAL*4 S4
      REAL*4 R4BUF(N)
      INTEGER N,K
      EQUIVALENCE(C4,S4)

C  Loop through the sets of four unformatted bytes.
      DO K=1,N

C  Copy the input four-byte floating point number to equivalenced
C  memory.
         S4=R4BUF(K)

C  Copy equivalenced memory to the four output bytes.
         STRING(4*K-3:4*K)=C4(1:4)
      ENDDO

C  Done.
      RETURN
      END
