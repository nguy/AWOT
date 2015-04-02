      SUBROUTINE READ_I2(STRING,I2)

C  Thomas Matejka NOAA/NSSL 25 May 1994

C  This subroutine interprets two unformatted bytes as a two-byte
C  integer.

C  Input:

C  STRING is a character string whose first two unformatted bytes are to
C  be interpreted as a two-byte integer.

C  Output:

C  I2 is the two-byte integer that corresponds to the two unformatted
C  bytes.

      IMPLICIT NONE
      CHARACTER*2 C2
      CHARACTER*(*) STRING
      INTEGER*2 I2,J2
      EQUIVALENCE(C2,J2)

C  Copy two input bytes to equivalenced memory.
      C2(1:2)=STRING(1:2)

C  Copy equivalenced memory to the output two-byte integer.
      I2=J2

C  Done.
      RETURN
      END
