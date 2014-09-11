      SUBROUTINE WRITE_R4(R4,STRING)

C  Thomas Matejka NOAA/NSSL 9 February 1996

C  This subroutine interprets a four-byte floating point number as four
C  unformatted bytes.

C  Input:

C  R4 is a four-byte floating point number to be interpreted as four
C  unformatted bytes.

C  Output:

C  STRING is the character string whose first four unformatted bytes
C  correspond to the four-byte floating point number.

      IMPLICIT NONE
      CHARACTER*4 C4
      CHARACTER*(*) STRING
      REAL*4 R4,S4
      EQUIVALENCE(C4,S4)

C  Copy the input four-byte integer to equivalenced memory.
      S4=R4

C  Copy equivalenced memory to the four output bytes.
      STRING(1:4)=C4(1:4)

C  Done.
      RETURN
      END
