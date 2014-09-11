      LOGICAL FUNCTION COMPARE_STRINGS(STRING1,STRING2)

C  Thomas Matejka NOAA/NSSL 26 March 1993

C  This function returns .TRUE. if and only if the string STRING1 to the
C  last non-blank character and the string STRING2 to the last non-blank
C  character are identical.

      IMPLICIT NONE
      INTEGER S_L
      CHARACTER*(*) STRING1,STRING2
      INTEGER IEND1,IEND2

      IEND1=S_L(STRING1)
      IEND2=S_L(STRING2)
      IF(STRING1(1:IEND1).EQ.STRING2(1:IEND2))THEN
         COMPARE_STRINGS=.TRUE.
      ELSE
         COMPARE_STRINGS=.FALSE.
      ENDIF
      RETURN
      END
