      SUBROUTINE RNGN(INSEED,RN_MEAN,RN_STDEV,RN)

C  Thomas Matejka NOAA/NSSL 11 November 1993

C  This subroutine generates RN, a real random number normally
C  distributed around RN_MEAN with a standard devaition RN_STDEV.
C  INSEED is an integer that is used to initialize the algorithm the
C  first time the subroutine is called.

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      LOGICAL FIRST_TIME
      INTEGER*4 INSEED,SEED,M,I
      REAL RN_MEAN,RN_STDEV,RN,A,X1,X2
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
      X1=FLOAT(SEED)/A
      SEED=MOD(I*SEED,M)
      X2=FLOAT(SEED)/A
      RN=RN_STDEV*SQRT(-2.*ALOG(X1))*COS(2.*PI*X2)+RN_MEAN
      RETURN
      END
