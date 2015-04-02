      Function TIMZ(IH,IM,IS)          ! Time: IH,IM,IS->seconds

C Function to convert time from 3 integers (IH,IM,IS) to seconds.

      TIMZ = Float(IH)*3600.0 + Float(IM)*60.0 + Float(IS)

      Return

      End
