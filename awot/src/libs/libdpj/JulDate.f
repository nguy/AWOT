      INTEGER FUNCTION JulDate (YEAR, MONTH, DAY)
C
C---COMPUTES THE JULIAN DATE (JD) GIVEN A GREGORIAN CALENDAR
C   DATE (YEAR,MONTH,DAY) from 1 January 4713 BC
C   From http://aa.usno.navy.mil/faq/docs/JD_Formula.html
C
      INTEGER YEAR, MONTH, DAY, I, J, K

      I= YEAR
      J= MONTH
      K= DAY

      JD= K-32075+1461*(I+4800+(J-14)/12)/4+367*(J-2-(J-14)/12*12) 
     2    /12-3*((I+4900+(J-14)/12)/100)/4

      Juldate = JD

      RETURN
      END

