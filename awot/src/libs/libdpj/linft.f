      Subroutine Linft (X, Y, Npts, A, SigmaA, B, SigmaB, R, Bad)

C Routine TO MAKE A LEAST SQUARES LINEAR FIT OF X AND Y
c to the formula y = A + Bx

      DOUBLE PRECISION Sum, SumX, SumY, SumX2, SumXY, SumY2
      DOUBLE PRECISION X1, Y1, DELTA, VARNCE

      Dimension X(Npts), Y(Npts)

      Sum = 0.0
      SumX = 0.0
      SumY = 0.0
      SumX2 = 0.0
      SumXY = 0.0
      SumY2 = 0.0
      Mpts  = 0

      Do i = 1, Npts
         If (X(i) .eq. Bad .or. Y(i) .eq. Bad) Cycle
         Mpts = Mpts + 1
         X1 = X(I)
         Y1 = Y(I)
         Sum = Sum + 1.0
         SumX = SumX + X1
         SumY = SumY + Y1
         SumX2 = SumX2 + X1*X1
         SumY2 = SumY2 + Y1*Y1
         SumXY = SumXY + X1*Y1
      End Do

      If (Mpts .lt. 3) Then
         A = Bad
         B = Bad
         SigmaA = Bad
         SigmaB = Bad
         R = Bad
      Else
         DELTA = Sum*SumX2 - SumX*SumX
         A = (SumX2*SumY-SumX*SumXY)/DELTA
         B = (SumXY*Sum-SumX*SumY)/DELTA
         
         C = Mpts-2
         VARNCE = (SumY2+A*A*Sum+B*B*SumX2-2.0*(A*SumY+B*SumXY-A*B*SumX))/C
         SIGMAA = Dsqrt(Dabs(VARNCE*SumX2/DELTA))
         SIGMAB = Dsqrt(Dabs(VARNCE*Sum  /DELTA))
         R = (Sum*SumXY-SumX*SumY)/DSQRT(DELTA*(Sum*SumY2-SumY*SumY))
      End If

      Return
      End
