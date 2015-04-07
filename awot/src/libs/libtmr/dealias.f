      SUBROUTINE DEALIAS(V_IN,V_NYQUIST,V_CENTER,V_OUT,IFOLD)

C  Thomas Matejka 10 January 1995

C  This subroutine dealiases a value.

C  Input:

C  V_IN is a real variable that specifies a possibly aliased value.

C  V_NYQUIST is a real variable that specifies one half of the total
C  unambiguous range of V_IN.

C  V_CENTER is a real variable that specifies the center of the
C  unambiguous range to dealias V_IN into.

C  Output:

C  V_OUT is a real variable that returns the dealiased value of V_IN.

C  IFOLD is an integer variable that specifies the number of unambiguous
C  ranges added to V_IN to dealias it.

      IMPLICIT NONE
      INTEGER IFOLD
      REAL V_IN,V_OUT,V_CENTER,V_NYQUIST,X

      X=(V_CENTER-V_IN)/(2.*V_NYQUIST)
      IFOLD=IFIX(X+SIGN(0.5,X))
      IF(AMOD(V_CENTER-V_IN-SIGN(V_NYQUIST,V_CENTER-V_IN),2.*V_NYQUIST)
     $.EQ.0.)THEN
         IFOLD=IFOLD-IFIX(SIGN(1.,X))
      ENDIF
      V_OUT=V_IN+FLOAT(IFOLD)*2.*V_NYQUIST
      RETURN
      END
