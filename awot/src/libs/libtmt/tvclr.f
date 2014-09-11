      REAL FUNCTION TVCLR(Z1,TV1,Z2,TV2,Z)

C  Thomas Matejka NOAA/NSSL 9 March 1993

      IMPLICIT NONE
      INCLUDE 'tmtlib.inc'
      REAL Z1,TV1,Z2,TV2,Z,GAMMA

      IF(Z2.NE.Z1)THEN
         IF(TV2.EQ.TV1)THEN
            TVCLR=TV1
         ELSE
            GAMMA=(TV2-TV1)/(Z2-Z1)
            TVCLR=TV1+GAMMA*(Z-Z1)
         ENDIF
      ELSE
         WRITE(TMTLIB_MESSAGE_UNIT,*)'TVCLR:  Z1 AND Z2 MUST NOT BE ',
     $   'EQUAL.'
         STOP
      ENDIF
      RETURN
      END
