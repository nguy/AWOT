      Function EPOTT(TEMP,DEWP,PRESS)
      WMR = 6.11*(10.0**(7.5*DEWP/(237.3+DEWP)))
      WMR = 621.98*WMR/(PRESS - WMR)
      TL = TLCL2(TEMP,PRESS,WMR)+273.16
      A = (TEMP+273.16)*((1000.0+1.6078*WMR)/PRESS)**0.28544
      B = EXP((3.1329-0.00237*TL)*(WMR/TL))
      C = EXP(EXP(1.62*ALOG(A*B)+14.3*ALOG(TL)-96.0))
      EPOTT = A*B*C - 273.16
      Return
      End
