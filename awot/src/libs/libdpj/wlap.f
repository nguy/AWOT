      Function WLAP(TEMP,PRESS)
C
      DOUBLE PRECISION F0,F1,F2
C
C  TEMP  = TEMPERATURE                                    (DEG C)
C  PRESS = PRESSURE                                       (MB)
C  WLAP  = MOIST LAPSE RATE                               (DEG C/MB)
      H =  (2500.-2.274*TEMP)*1000.0
      ES = VAPOR(TEMP)
      WS = 0.62198*ES/(PRESS-ES)
      F0 = H/(TEMP+273.16)
      F1 = 286.998+WS*F0
      F2 = 1004.*286.998+0.62198*WS*F0*F0
      WLAP = F1*286.998*(TEMP+273.16)/(PRESS*F2)
      Return
      End
