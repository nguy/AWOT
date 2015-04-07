      SUBROUTINE RNGAB(INSEED,RN_MIN,RN_MAX,RN)

C  Thomas Matejka NOAA/NSSL 11 November 1993

C  This subroutine generates RN, a real random number evenly distributed
C  between RN_MIN and RN_MAX.  INSEED is an integer that is used to
C  initialize the algorithm the first time the subroutine is called.

      IMPLICIT NONE
      LOGICAL FIRST_TIME
      INTEGER*4 INSEED,SEED,M,I
      REAL RN_MIN,RN_MAX,RN,A
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
      RN=FLOAT(SEED)*(RN_MAX-RN_MIN)/A+RN_MIN
      RETURN
      END
