      SUBROUTINE READ_I4_BUF(STRING,N,I4BUF)

C  Thomas Matejka NOAA/NSSL 9 February 1996

C  This subroutine interprets a series of sets of four unformatted bytes
C  as a one-dimensional array of four-byte integers.

C  Input:

C  STRING is a character string whose first 4*N unformatted bytes are to
C  be interpreted as N four-byte integers.

C  N is an integer variable that specifies the number of sets of four
C  unformatted bytes in STRING and the number of elements in I4BUF.

C  Output:

C  I4BUF is a one-dimensional integer*4 array.  I4BUF(I) is the Ith
C  four-byte integer that corresponds to the Ith set of four unformatted
C  bytes.

      IMPLICIT NONE
      CHARACTER*4 C4
      CHARACTER*(*) STRING
      INTEGER*4 J4
      INTEGER*4 I4BUF(N)
      INTEGER N,K
      EQUIVALENCE(C4,J4)

C  Loop through the sets of four unformatted bytes.
      DO K=1,N

C  Copy four input bytes to equivalenced memory.
         C4(1:4)=STRING(4*K-3:4*K)

C  Copy equivalenced memory to the output four-byte integer.
         I4BUF(K)=J4
      ENDDO

C  Done.
      RETURN
      END
