      SUBROUTINE WRITE_I2(I2,STRING)

C  Thomas Matejka NOAA/NSSL 9 February 1996

C  This subroutine interprets a two-byte integer as two unformatted
C  bytes.

C  Input:

C  I2 is a two-byte integer to be interpreted as two unformatted bytes.

C  Output:

C  STRING is the character string whose first two unformatted bytes
C  correspond to the two-byte integer.

      IMPLICIT NONE
      CHARACTER*2 C2
      CHARACTER*(*) STRING
      INTEGER*2 I2,J2
      EQUIVALENCE(C2,J2)

C  Copy the input two-byte integer to equivalenced memory.
      J2=I2

C  Copy equivalenced memory to the two output bytes.
      STRING(1:2)=C2(1:2)

C  Done.
      RETURN
      END
