      Subroutine Tcnvt(Itime, Ih, Im, Is) ! Time: hhmmss->HH MM SS

C Routine to convert time from hhmmss. format to hours, min, & secs

      IH = ITIME/10000
      IM = (ITIME - IH*10000)/100
      IS = ITIME - IH*10000 - IM*100

      Return
      End
