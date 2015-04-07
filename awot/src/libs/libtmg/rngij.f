      SUBROUTINE RNGIJ(INSEED,IRN_MIN,IRN_MAX,IRN)

C  Thomas Matejka NOAA/NSSL 11 November 1993

C  This subroutine generates IRN, an integer random number evenly
C  distributed between IRN_MIN and IRN_MAX.  INSEED is an integer that
C  is used to initialize the algorithm the first time the subroutine is
C  called.

      IMPLICIT NONE
      LOGICAL FIRST_TIME
      INTEGER*4 INSEED,SEED,M,I,IRN_MIN,IRN_MAX,IRN
      REAL A,RN
      DATA FIRST_TIME/.TRUE./
      SAVE FIRST_TIME,SEED,M,A,I

      IF(FIRST_TIME)THEN
         FIRST_TIME=.FALSE.
         M=2**20
         A=FLOAT(M)
         I=2**10+3
         SEED=100001+IFIX((999999.-100001.)*FLOAT(INSEED-1)/(32767.-1.))
         IF(MOD(SEED,2).EQ.0)THEN
            SEED=SEED+1
         ENDIF
      ENDIF
      SEED=MOD(I*SEED,M)
      RN=FLOAT(SEED)*FLOAT(IRN_MAX-IRN_MIN+1)/A+FLOAT(IRN_MIN)-0.5
      IRN=IFIX(RN+SIGN(0.5,RN))
      RETURN
      END

