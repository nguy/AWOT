      Function EPOT3(TEMP,PRESS,WMR)

!     Returns EQUIVALENT POTENTIAL TEMPERATURE               (DEG C)
!     Given termperature (C0, Pressure [mb], and mixing ration [gm/gm]
      DATA M, C1, C2Kelvin /1, 2.64, 273.16/

      R = 621.98*((WMR*PRESS)/621.98)/(PRESS-(WMR*PRESS)/621.98)
      EE=(WMR*PRESS)/621.98

!     TETONS FORMULA

      THETA = (TEMP + C2Kelvin) * ((1000./PRESS)**0.28544)

      DO j = 1, 2
         THETE=THETA*EXP(C1*R/(TLCL2(TEMP,PRESS,WMR)+C2Kelvin))
         C1=2.627+.0003*(THETE-300.)
      End Do

      EPOT3 = THETE - C2Kelvin
      Return

      End
