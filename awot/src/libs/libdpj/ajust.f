      Function AJUST(TEMP,REFL,PALT,TLAP)

C   Returns A TEMPERATURE ADJUSTED TO A REFERENCE LEVEL  IN     (DEG C)
C   TEMP =   TEMPERATURE TO BE ADJUSTED                         (DEG C)
C   REFL =   REFERENCE LEVEL                                    (M)
C   PALT =   PRESSURE ALTITUDE                                  (M)
C   TLAP =   LAPSE RATED OF TEMPERATURE                         (DEG C/M)
C   TLAP =   0 DEFAULT VALUE,0.005577, IS USED

      If (TLAP.eq.0) TLAP=0.005577
      AJUST  = TEMP + TLAP*(PALT-REFL)
      Return
      End
