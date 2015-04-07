      Function WBULB(TEMP,PRESS,DP,PSTOP)

C   Returns WETBULB POTENTIAL  TEMPERATURE                 (DEG C)
C                                        IF PSTOP = 1000.0 (MB)
C   Returns            WETBULB TEMPERATURE                 (DEG C)
C                                        IF PSTOP = PRESS  (MB)
C   TEMP  = TEMPERATURE                                    (DEG C)
C   PRESS = PRESSURE                                       (MB)
C   DP    = DEWPOINT                                       (DEG C)
C   PSTOP = PRESSURE AT WHICH CALCULATION STOPS            (MB)

      If (ABS(TEMP-DP).le.0.5) Go To 8
      If (DP.gt.TEMP) Go To 3
      W = FMIXD(DP,PRESS)
      THETA = POTM(TEMP,PRESS,W) + 273.16
      TL = TLCL2(TEMP,PRESS,W) +273.16
      If (TL .ge. 1.E36 .OR. W .le .0.0) Go To 3
      PL = 1000.*(TL/THETA)**3.487796
      Go To 9
    8 Continue
      PL = PRESS
      TL = TEMP + 273.16
    9 Continue
      T1 = TL
      P1 = PL
      If (ABS (PL-1000.0) .le.0.5) Go To 4
      If ((PSTOP-PL).gt.0.5) SGN = 1.0
      If ((PL-PSTOP).gt.0.5) SGN =-1.0
      DO 200 L = 1,200
      DIFF = ABS(PSTOP-P1)*SGN
      If (DIFF.le.0.5) Go To 5
      If (DIFF.le.5.0) DDP = DIFF
      If (DIFF.gt.5.0) DDP = SGN*5.
      If (DIFF.gt.10.) DDP = SGN*10.
      If (DIFF.gt.25.) DDP = SGN*25.
      P2 = P1  +  DDP
      A1 = DDP*WLAP(T1,P1)
      B1 = DDP*WLAP(T1+A1/2.,P1+DDP/2.)
      C1 = DDP*WLAP(T1+B1/2.,P1+DDP/2.)
      D1 = DDP*WLAP(T1+C1,P2)
      T2 = T1 + A1/6.+B1/3.+C1/3.+D1/6.
      T1 = T2
      P1 = P2
  200 Continue
    3 WBULB = 1.E36
      Return
    4 T2 = TL
    5 WBULB = T2 -273.16
      Return
      End
