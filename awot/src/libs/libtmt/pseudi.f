      SUBROUTINE PSEUDI(T,P,DTDPI,DQSIDTI,DQSIDPI)

C  Thomas Matejka NOAA/NSSL 18 June 1993

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      INCLUDE 'tmtlib.inc'
      REAL VAPPSI,ALHSUB
      REAL T,P,DTDPI,DQSIDTI,DQSIDPI,ESI,ALS,A,B,C,D

      IF(T.LE.T_FREEZE)THEN
         ESI=VAPPSI(T)
         IF(ESI.GE.P)THEN
            WRITE(TMTLIB_MESSAGE_UNIT,*)'PSEUDI:  ESI IS GREATER THAN ',
     $      'OR EQUAL TO P.'
            STOP
         ENDIF
         ALS=ALHSUB(T)
         A=UNIV_GAS*T/(CP_DRY*(P-ESI)*MW_DRY+CP_VAP*ESI*MW_H2O)
         B=-(ALS+UNIV_GAS*T*(1./MW_H2O-1./MW_DRY))*
     $   ((P-ESI)*MW_DRY+ESI*MW_H2O)/
     $   (CP_DRY*(P-ESI)*MW_DRY+CP_VAP*ESI*MW_H2O)
         C=-MW_H2O*MW_DRY*ESI/((P-ESI)*MW_DRY+ESI*MW_H2O)**2
         D=P*MW_DRY*MW_H2O*ALS/((P-ESI)*MW_DRY+ESI*MW_H2O)**2/T/
     $   (UNIV_GAS*T/MW_H2O/ESI-1./RHO_ICE)
         DTDPI=(A+B*C)/(1.-B*D)
         DQSIDTI=(C+A*D)/(A+B*C)
         DQSIDPI=(C+A*D)/(1.-B*D)
      ELSE
         DTDPI=TMTLIB_BADFLAG
         DQSIDTI=TMTLIB_BADFLAG
         DQSIDPI=TMTLIB_BADFLAG
      ENDIF
      RETURN
      END
