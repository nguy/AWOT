      SUBROUTINE WRITE_I2_BUF(I2BUF,N,STRING)

C  Thomas Matejka NOAA/NSSL 9 February 1996

C  This subroutine interprets a one-dimensional array of two-byte
C  integers as a series of sets of two unformatted bytes.

C  Input:

C  I2BUF is a one-dimensional integer*2 array.  I2BUF(I) is the Ith
C  two-byte integer to be interpreted as the Ith set of two unformatted
C  bytes.

C  N is an integer variable that specifies the number of elements in
C  I2BUF and the number of sets of two unformatted bytes in STRING.

C  Output:

C  STRING is the character string whose first 2*N unformatted bytes
C  correspond to N two-byte integers.

      IMPLICIT NONE
      CHARACTER*2 C2
      CHARACTER*(*) STRING
      INTEGER*2 J2
      INTEGER*2 I2BUF(N)
      INTEGER N,K
      EQUIVALENCE(C2,J2)

C  Loop through the sets of two unformatted bytes.
      DO K=1,N

C  Copy the input two-byte integer to equivalenced memory.
         J2=I2BUF(K)

C  Copy equivalenced memory to the two output bytes.
         STRING(2*K-1:2*K)=C2(1:2)
      ENDDO

C  Done.
      RETURN
      END
