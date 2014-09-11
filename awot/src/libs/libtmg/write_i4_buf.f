      SUBROUTINE WRITE_I4_BUF(I4BUF,N,STRING)

C  Thomas Matejka NOAA/NSSL 9 February 1996

C  This subroutine interprets a one-dimensional array of four-byte
C  integers as a series of sets of four unformatted bytes.

C  Input:

C  I4BUF is a one-dimensional integer*4 array.  I4BUF(I) is the Ith
C  four-byte integer to be interpreted as the Ith set of four
C  unformatted bytes.

C  N is an integer variable that specifies the number of elements in
C  I4BUF and the number of sets of four unformatted bytes in STRING.

C  Output:

C  STRING is the character string whose first 4*N unformatted bytes
C  correspond to N four-byte integers.

      IMPLICIT NONE
      CHARACTER*4 C4
      CHARACTER*(*) STRING
      INTEGER*4 J4
      INTEGER*4 I4BUF(N)
      INTEGER N,K
      EQUIVALENCE(C4,J4)

C  Loop through the sets of four unformatted bytes.
      DO K=1,N

C  Copy the input four-byte integer to equivalenced memory.
         J4=I4BUF(K)

C  Copy equivalenced memory to the four output bytes.
         STRING(4*K-3:4*K)=C4(1:4)
      ENDDO

C  Done.
      RETURN
      END
