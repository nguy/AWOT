      SUBROUTINE READ_I1(C,I4)

C  Thomas Matejka NOAA/NSSL 14 November 1996

C  This subroutine interprets one unformatted byte as a four-byte
C  integer.

C  Input:

C  C is a character variable that is to be interpreted as an integer.

C  Output:

C  I4 is the four-byte integer that corresponds to the unformatted byte.

      IMPLICIT NONE
      CHARACTER(LEN=1)::C
      CHARACTER(LEN=4)::C4
      INTEGER*4 I4,J4
      EQUIVALENCE(C4,J4)

C  Copy the input byte to equivalenced memory.
      J4=0
      C4(4:4)=C

C  Copy equivalenced memory to the output four-byte integer.
      I4=J4

C  Done.
      RETURN
      END
