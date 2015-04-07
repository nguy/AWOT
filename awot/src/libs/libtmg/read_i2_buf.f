      SUBROUTINE READ_I2_BUF(STRING,N,I2BUF)

C  Thomas Matejka NOAA/NSSL 9 February 1996

C  This subroutine interprets a series of sets of two unformatted bytes
C  as a one-dimensional array of two-byte integers.

C  Input:

C  STRING is a character string whose first 2*N unformatted bytes are to
C  be interpreted as N two-byte integers.

C  N is an integer variable that specifies the number of sets of two
C  unformatted bytes in STRING and the number of elements in I2BUF.

C  Output:

C  I2BUF is a one-dimensional integer*2 array.  I2BUF(I) is the Ith
C  two-byte integer that corresponds to the Ith set of two unformatted
C  bytes.

      IMPLICIT NONE
      CHARACTER*2 C2
      CHARACTER*(*) STRING
      INTEGER*2 J2
      INTEGER*2 I2BUF(N)
      INTEGER N,K
      EQUIVALENCE(C2,J2)

C  Loop through the sets of two unformatted bytes.
      DO K=1,N

C  Copy two input bytes to equivalenced memory.
         C2(1:2)=STRING(2*K-1:2*K)

C  Copy equivalenced memory to the output two-byte integer.
         I2BUF(K)=J2
      ENDDO

C  Done.
      RETURN
      END
