      SUBROUTINE READ_I4(STRING,I4)

C  Thomas Matejka NOAA/NSSL 25 May 1994

C  This subroutine interprets four unformatted bytes as a four-byte
C  integer.

C  Input:

C  STRING is a character string whose first four unformatted bytes are
C  to be interpreted as a four-byte integer.

C  Output:

C  I4 is the four-byte integer that corresponds to the four unformatted
C  bytes.

      IMPLICIT NONE
      CHARACTER*4 C4
      CHARACTER*(*) STRING
      INTEGER*4 I4,J4
      EQUIVALENCE(C4,J4)

C  Copy four input bytes to equivalenced memory.
      C4(1:4)=STRING(1:4)

C  Copy equivalenced memory to the output four-byte integer.
      I4=J4

C  Done.
      RETURN
      END
