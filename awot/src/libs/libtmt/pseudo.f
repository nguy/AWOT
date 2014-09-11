      SUBROUTINE PSEUDO(T,P,DTDPP,DQSDTP,DQSDPP)

C  Thomas Matejka NOAA/NSSL 18 June 1993

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      INCLUDE 'tmtlib.inc'
      REAL VAPPS,ALHVAP
      REAL T,P,DTDPP,DQSDTP,DQSDPP,ES,ALV,A,B,C,D

      ES=VAPPS(T)
      IF(ES.GE.P)THEN
         WRITE(TMTLIB_MESSAGE_UNIT,*)'PSEUDO:  ES IS GREATER THAN OR ',
     $   'EQUAL TO P.'
         STOP
      ENDIF
      ALV=ALHVAP(T)
      A=UNIV_GAS*T/(CP_DRY*(P-ES)*MW_DRY+CP_VAP*ES*MW_H2O)
      B=-(ALV+UNIV_GAS*T*(1./MW_H2O-1./MW_DRY))*
     $((P-ES)*MW_DRY+ES*MW_H2O)/
     $(CP_DRY*(P-ES)*MW_DRY+CP_VAP*ES*MW_H2O)
      C=-MW_H2O*MW_DRY*ES/((P-ES)*MW_DRY+ES*MW_H2O)**2
      D=P*MW_DRY*MW_H2O*ALV/((P-ES)*MW_DRY+ES*MW_H2O)**2/T/
     $(UNIV_GAS*T/MW_H2O/ES-1./RHO_WAT)
      DTDPP=(A+B*C)/(1.-B*D)
      DQSDTP=(C+A*D)/(A+B*C)
      DQSDPP=(C+A*D)/(1.-B*D)
      RETURN
      END
