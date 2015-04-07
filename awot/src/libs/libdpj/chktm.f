      Subroutine CHKTM(T)
  
C   Routine to check for time flips over midnight
C   this will add 24 hours if a flip occurs

      DATA OCN /86400.0/, TMJMP /0.0/, TH /0.0/
      Save Ocn, Th, Tmjmp

      IF (T .lt. 0.0) Return

      IF (T-TH .lt. -86300.0) Go To 1
      TH = T
      T = T + TMJMP
      Return

    1 TMJMP = OCN
      TH = T
      T = T + TMJMP

      Return
      End
