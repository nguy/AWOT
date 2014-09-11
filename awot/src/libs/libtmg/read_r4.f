      SUBROUTINE READ_R4(STRING,R4)

C  Thomas Matejka NOAA/NSSL 25 May 1994

C  This subroutine interprets four unformatted bytes as a four-byte
C  floating-point number.

C  Input:

C  STRING is a character string whose first four unformatted bytes are
C  to be interpreted as a four-byte floating-point number.

C  Output:

C  R4 is the four-byte floating-point number that corresponds to the
C  four unformatted bytes.

      IMPLICIT NONE
      CHARACTER*4 C4
      CHARACTER*(*) STRING
      REAL*4 R4,S4
      EQUIVALENCE(C4,S4)

C  Copy four input bytes to equivalenced memory.
      C4(1:4)=STRING(1:4)

C  Copy equivalenced memory to the output four-byte floating-point
C  number.
      R4=S4

C  Done.
      RETURN
      END
