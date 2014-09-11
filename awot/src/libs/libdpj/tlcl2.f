      Function TLCL2(TEMP,PRESS,WMR)

      Dimension DT(8)
      DATA DT/1.5,.75,0.37,0.2,0.1,0.05,0.025,0.005/
      DATA RES/1.0/,R1/1.0/,R2/1.0/,I/1/
      T = TEMP + 273.16
      E = WMR*0.001*PRESS/(WMR*0.001+0.62198)
      A = 3.5*ALOG(T)-ALOG(E)-4.805
      TC= 2840.0/A+55.0
      RKM = 1.0/(0.28544*(1.0-0.245*WMR*0.001))
      DO 50 K = 1,500
      R1 = SIGN(R1,RES)
      X   = 6.11*(WMR+621.95)/(WMR*PRESS)
      Y   = (T/TC)**RKM
      Z   = 7.5*(TC-273.15)/(TC-35.85)
      RES = ALOG10(X*Y) + Z
      If (ABS(RES).le.0.005) Go To 100
      R2 = SIGN(R2,RES)
      If (SIGN(R1,R2).ne.SIGN(R2,R1)) I = I+1
      If (I.gt.8) I = 8
      TC = TC -R2*DT(I)
   50 Continue
  100 TLCL2 = TC-273.16
      Return
      End
