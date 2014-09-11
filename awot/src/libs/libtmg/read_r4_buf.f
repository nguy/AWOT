      SUBROUTINE READ_R4_BUF(STRING,N,R4BUF)

C  Thomas Matejka NOAA/NSSL 9 February 1996

C  This subroutine interprets a series of sets of four unformatted bytes
C  as a one-dimensional array of four-byte floating-point numbers.

C  Input:

C  STRING is a character string whose first 4*N unformatted bytes are to
C  be interpreted as N four-byte floating-point numbers.

C  N is an integer variable that specifies the number of sets of four
C  unformatted bytes in STRING and the number of elements in R4BUF.

C  Output:

C  R4BUF is a one-dimensional real*4 array.  R4BUF(I) is the Ith
C  four-byte floating-point number that corresponds to the Ith set of
C  four unformatted bytes.

      IMPLICIT NONE
      CHARACTER*4 C4
      CHARACTER*(*) STRING
      REAL*4 S4
      REAL*4 R4BUF(N)
      INTEGER N,K
      EQUIVALENCE(C4,S4)

C  Loop through the sets of four unformatted bytes.
      DO K=1,N

C  Copy four input bytes to equivalenced memory.
         C4(1:4)=STRING(4*K-3:4*K)

C  Copy equivalenced memory to the output four-byte floating-point
C  number.
         R4BUF(K)=S4
      ENDDO

C  Done.
      RETURN
      END
