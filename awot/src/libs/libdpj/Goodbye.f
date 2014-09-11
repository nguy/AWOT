C  "<08-Oct-1993 21:54:30><UTC>"
********************************************************************
      Subroutine Goodbye(LuHomeW)
c  Say goodbye and exit.
      Implicit none
      Integer*4 LuhomeW
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      Write (LuhomeW,50001)
      Stop
50001 Format(/' Exit - Hasta luego, Amigo mio! - A bientot, mon ami!'/
     >' Tschuess, mein Freund - Assalaamu ''alaykum - Shalom aleichem'/
     >' Dzai jian - sayonnara - Dag, mijn vriend! - Farvel, min ven!'/
     >' Bless-Bless - Good-bye')
      End ! subroutine goodbye ends
