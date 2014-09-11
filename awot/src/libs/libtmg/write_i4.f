      SUBROUTINE WRITE_I4(I4,STRING)

C  Thomas Matejka NOAA/NSSL 9 February 1996

C  This subroutine interprets a four-byte integer as four unformatted
C  bytes.

C  Input:

C  I4 is a four-byte integer to be interpreted as four unformatted
C  bytes.

C  Output:

C  STRING is the character string whose first four unformatted bytes
C  correspond to the four-byte integer.

      IMPLICIT NONE
      CHARACTER*4 C4
      CHARACTER*(*) STRING
      INTEGER*4 I4,J4
      EQUIVALENCE(C4,J4)

C  Copy the input four-byte integer to equivalenced memory.
      J4=I4

C  Copy equivalenced memory to the four output bytes.
      STRING(1:4)=C4(1:4)

C  Done.
      RETURN
      END
