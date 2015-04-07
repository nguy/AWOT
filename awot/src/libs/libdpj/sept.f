      Function Sept (Temp,Dp,Press)

c  Routine to compute equivalent potential temperature using
c  the Simpson method (handles tropical air better)
c  Input variables are Temp = Temperature [C], Dp = Dewpoint temperature [C]
c  and Press = Pressure [mb]

      If (Temp .lt. -40.0 .or. Temp .gt. 100.0) Go To 1
      If (Dp   .lt. -40.0 .or. Dp   .gt. 100.0) Go To 1
      If (Press .lt. 100.0 .or. Press .gt. 1100.0) Go To 1

      Q=6.11*(10.0**(7.5*DP/(237.3+DP)))
      WMR=621.98*Q/(PRESS-Q)
      TL = TLCL2(TEMP,PRESS,WMR)+273.16
      A = (TEMP+273.16)*((1000.0+1.6078*WMR)/PRESS)**0.28544
      B = EXP((3.1329-0.00237*TL)*(WMR/TL))
      C = EXP(EXP(1.62*ALOG(A*B)+14.3*ALOG(TL)-96.0))
      Sept  = A*B*C
      Return

    1 Sept = 1.E36
      Return

      End
