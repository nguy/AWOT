      Subroutine Clot(Z, L1, K, Cint, Sf1, Spi1, Spj1, Ijk, Idash1, Lab)
c
c     Contour program with equally spaced contours given by cinter-set
c     condat(1)=0.0-or assign contour values by placing them in condat.
c     Contour values must be stored in condat in assending order
c     beginning with any value except 0.0. be certain to supply values
c     both larger and smaller then expected max. and min. contours.
c     Contour values should be floating integers because the label code
c     will not plot fractional parts.
c     Cinter(contour interval) must be assigned a value. All contours
c     that are not integer multiples of cinter are tested for contour
c     spacing and if contours are too close that part of the contour is
c     not plotted. To suppress this test set cinter to 1.0.
c     the program can draw contours that are 1/2, 1/3, 1/4, 1/5, 1/6,
c     1/8, 1/10, 1/16, 1/20, 1/25, 1/32, 1/40, or 1/50 of cinter
c     if logarithmic contour spaces are required,fractional parts are
c     usually necessary. The code should be modified to by pass the
c     spacing test. Also contour labels will truncate the fractional parts.
c     Data to be contoured are placed in Z(L,K) array. Data are assumed
c     to be stored in colum beginning with the lower left corner. L is
c     the colume length and K the row(left to right) length.
c     Sf is the scale factor.  WARNING: all values in Z are multiplied
c         by Sf; on normal exit Z values are divided by Sf: this could
c         sligthly change the Z values Sf <> 1.
c     Spi and Spj are the grid interval spacings in inches for the
c     calcomp plotter. For the stromberg carlson 4020 microfilm they may
c     also be thought of as inches. The default value of the maximum
c     grid size for 35 mm frames is 11 inches and Spi and Spj must be
c     selected fo fit your grid into this space
c              Spj=x size in inches/number of x grid points
c              Spi=y size in inches/number of y grid points
c     Spi is the grid interval up and down and spj is left and right.
c     The parameter ijk is the number of interpolated points on a grid
c     square in both directions. The maximum is 15 and the minimum is 2.
c     The program runs much faster for small ijk but the contours are
c     not so well curved.  Experimentation is suggested.
c     This program produces dashed contours controlled by Idash
c     Idash=0 no dashing
c     Idash=1 all dashed
c     Idash=2 zero contour dashed
c     Idash=3 contours lt 0 dashed
c     Idash=4 contours gt 0 dashed
c     Idash=5 contours le 0 dashed
c     Idash=6 contours ge 0 dashed
c     Idash=7 intermediate contours dashed
c     Idash lt 0 or gt 7 all dashed.
c     The ad tolarence must be increased for a 32 bit machine.
c     The plot commands may have to be modified for another system.
c     Klab is the labeling feature (0=no contour labels, 1=labels)
c         Klab=2 means only dashed contours labeled
c         Klab=3 means no dashed contours labeled
c
c     Jan 1988 mods by R. Hueftle
c         Subroutine Pltsg was modified to improve appearance of contours
c             close to labels.
c         Common block /Cont_Pen/Npen(1..3) is used for calls to Newpn
c             Npen(1) is used for negative contours.
c             Npen(2) is used for zero contours.
c             Npen(3) is used for positive contours.
c             Npen(1..3) are initialized in a data statement to all = 1;
c                 they can be changed outside of Clot (but NOT with data
c                 statements).  Npen(3) is used to set Newpn upon exit
c                 from Clot.
c
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      Real Max, Min

      Integer E0, F0, H, A1, B1, X1, Y1, X2, Y2, Side, Aside, Trip

      Common /Clotx/ Es, Fs, Ek, Fk, Spi, Spj, Kount, Ipen, Klab, Cvar,
     . Dnum, Idsh, Pn, Qn, Imumd, Mumd, Numd, Height
      Common /Contr/ Condat(31)
      Common /Numsz/ Hgt
      Common /Cont_Pen/ Npen(3)                                         !RAH
  
      Dimension Z(L1,K)                                                 !RAH
      Dimension Est(4001), Fst(4001)
     > ,  F(4,4), Px(4,15), Py(15,15), Xt(15)                           !RAH
     > ,  Tw(15), Tx(15), Ty(15), Tz(15)                                !RAH

      Save Len, Ni, Adj
  
      Data Len/3996/, Ni/10/, Adj/.0001/
  
      Pn=.1
      Qn=.075
  
      If (Hgt .Gt. 1.0 .Or. Hgt .Le. 0.0) Hgt=0.07
  
      Height=Hgt
      L = L1
      H = K
      Cinter = Cint
      Sf = Sf1
      Spi = Spi1
      Spj = Spj1
      Jk = Ijk
      Idash = Idash1
      Klab = Lab
      Ds = 1./(Float(Jk) - 1.)
      Imumd = 0
      Max = Spi
      If (Spj .Gt. Max) Max = Spj
      If (Abs(Spi - Spj) .gt. 0.2*Max) Imumd = 1
      Min = Z(1,1)*Sf
      Max = Min
  
      Do I = 1, L
         Do J = 1, H
            Z(I,J) = Z(I,J) * Sf
            If (Z(I,J) .Gt. Max) Max = Z(I,J)
            If (Z(I,J) .Lt. Min) Min = Z(I,J)
         End Do
      End Do

      Cva1 = Aint(Min/Cinter)*Cinter
      Icntst = 0
      Ctest = 0.
      Iii = L
      If (H .Gt. L) Iii = H
      Cspace = 200./Float(Iii)
      If (Condat(1) .Ne. 0.) Icntst = 1
      If (Icntst .Eq. 0) Go To 166

      If (Condat(1) - Min .Gt. Cinter) Write (7,161) Min
 161  Format (' Contours Missing About Min=',F8.1)
      If (Condat(31) .Lt. Max - Cinter) Write (7,163) Max
 163  Format (' Contours Missing Around Max=',F8.1)
      If (Condat(1) .Gt. Max) Write (7,164)
 164  Format (26H Contour Values Exceed Max)

      Do I = 3,31
         Istore = I - 1
         If (Condat(I) .Eq. Condat(I-1)) Go To 134
      End Do
 
 134  Do I = Istore, 31
          Condat(I) = Condat(I-1) + Cinter
      End Do

      Do I = 1, 31
         Cva1 = Condat(I)
         N = I
         If (Cva1 .Gt. Min) Go To 165
      End Do

 160  Write (7,162)
 162  Format (' Min Exceeds Contour Values; may be altered!!')       !RAH
      Call NewPn(NPen(3)) ! pen to use on exit from this subroutine  !RAH
      Return

 165  Icon = N

 166  Do I = 1,L
         Do J = 1,H
            Z(I,J) = Z(I,J) - Cva1
         End Do
      End Do

      Iii = Float(Jk)/Sqrt(Spi*Spi + Spj*Spj)
      Jo = Iii*Ni
      I3 = 0
 11   T1 = I3
      Cval = Aint(Min/Cinter)*Cinter + Cinter*T1 - Cva1
      If (Icntst .Eq. 0) Go To 190
      Iii = Icon + I3
      Cvar = Condat(Iii)
      Cval = Cvar - Cva1
      Ctest = 0.
      If (Aint(Cvar/Cinter)*Cinter .Ne. Cvar) Ctest = 1.
      If (Cinter .Eq. 1.0) Ctest = 0.0
      If (Ctest .Lt. 1.) Go To 190
      Sptest = Cinter*Cspace
      Tcvar = Amod(Cvar, Cinter)
      Tcvar = Abs(Tcvar)
      Xc = .03125*Cinter
      Tc = Tcvar/Xc

      Do I = 1, 31, 2
         If (Abs(Tc - Float(I)) .Lt. Adj)    Tcvar = .03125*Cinter
      End Do

      Do I = 2, 30, 4
         If (Abs(Tc - Float(I)) .Lt. Adj)    Tcvar = .0625*Cinter
      End Do

      Do I = 4, 28, 8
         If (Abs(Tc - Float(I)) .Lt. Adj)    Tcvar = .125*Cinter
      End Do

      Do I = 8, 24, 16
         If (Abs(Tc - Float(I)) .Lt. Adj)    Tcvar = .25*Cinter
      End Do

      If (Abs(Tc - 16.) .Lt. Adj)    Tcvar = .5*Cinter
      Xc = .025*Cinter
      Tc = Tcvar/Xc

      Do 179 I = 1,39,2
         If (Mod(I,5) .Eq. 0) Go To 179
         If (Abs(Tc - Float(I)) .Lt. Adj)    Tcvar = .025*Cinter
 179  Continue

      Do 180 I = 2,38,4
         If (Mod(I,5) .Eq. 0) Go To 180
         If (Abs(Tc - Float(I)) .Lt. Adj)    Tcvar = .05*Cinter
 180  Continue

      Do 181 I = 4,36,8
         If (I .Eq. 20) Go To 181
         If (Abs(Tc - Float(I)) .Lt. Adj)    Tcvar = .1*Cinter
 181  Continue

      Do 182 I = 8,32,8
         If (Abs(Tc - Float(I)) .Lt. Adj)    Tcvar = .2*Cinter
 182  Continue

      Xc = .1666666667*Cinter
      Tc = Tcvar/Xc

      Do 183 I = 1,5,4
         If (Abs(Tc - Float(I)) .Lt. Adj)    Tcvar = .1666666667*Cinter
 183  Continue

      Do 185 I = 2,4,2
         If (Abs(Tc - Float(I)) .Lt. Adj)    Tcvar = .3333333333*Cinter
 185  Continue

      Xc = .02*Cinter
      Tc = Tcvar/Xc

      Do 187 I = 1,49,2
         If (Mod(I,5) .Eq. 0) Go To 187
         If (Abs(Tc - Float(I)) .Lt. Adj)    Tcvar = .02*Cinter
 187  Continue

      Do 189 I = 2,48,2
         If (Mod(I,5) .Eq. 0) Go To 189
         If (Abs(Tc - Float(I)) .Lt. Adj)    Tcvar = .04*Cinter
 189  Continue

      If (Tcvar .Gt. .51*Cinter) Tcvar = .5*Cinter
      T1 = Cinter*Cspace
      If (Tcvar .Lt. .51*Cinter) Sptest = T1*.5
      If (Tcvar .Lt. .34*Cinter) Sptest = T1*.3333333333
      If (Tcvar .Lt. .26*Cinter) Sptest = T1*.25
      If (Tcvar .Lt. .21*Cinter) Sptest = T1*.2
      If (Tcvar .Lt. .17*Cinter) Sptest = T1*.1666666667
      If (Tcvar .Lt. .13*Cinter) Sptest = T1*.125
      If (Tcvar .Lt. .11*Cinter) Sptest = T1*.1
      If (Tcvar .Lt. .07*Cinter) Sptest = T1*.0625
      If (Tcvar .Lt. .06*Cinter) Sptest = T1*.05
      If (Tcvar .Lt. .049*Cinter) Sptest = T1*.04
      If (Tcvar .Lt. .04*Cinter) Sptest = T1*.03125
      If (Tcvar .Lt. .03*Cinter) Sptest = T1*.025
      If (Tcvar .Lt. .024*Cinter) Sptest = T1*.02

 190  Ad = Adj*Cval + Adj
      Cvar = Cval + Cva1
      If (I3 .Gt. 60) Go To 140
      If (Cvar .Gt. Max) Go To 5000

      Idsh = 1
      If (Idash .Eq. 0) Idsh = 0
      If (Idash .Eq. 2 .And. Cvar .Ne. 0.) Idsh = 0
      If (Idash .Eq. 3 .And. Cvar .Ge. 0.) Idsh = 0
      If (Idash .Eq. 4 .And. Cvar .Le. 0.) Idsh = 0
      If (Idash .Eq. 5 .And. Cvar .Gt. 0.) Idsh = 0
      If (Idash .Eq. 6 .And. Cvar .Lt. 0.) Idsh = 0
      If (Idash .Eq. 7 .And. Ctest .Lt. 1.) Idsh = 0
      T1 = Abs(Cvar)
  
c  Change pens to match new value of Cvar.                              !RAH
                                                                        !RAH
      If (Cvar.Lt.0.) Then                                              !RAH
          Call NewPn (Npen(1))                                          !RAH
      Else If (Cvar.Eq.0.0) Then                                        !RAH
          Call NewPn (Npen(2))                                          !RAH
      Else                                                              !RAH
          Call NewPn (Npen(3))                                          !RAH
      End If                                                            !RAH
  
      Numd = 5
      Dnum = 0.0
      If (T1 .Le. 1.) Go To 131
      Dnum = Alog10(T1)
      Dnum = Aint(Dnum)
      Numd = Dnum + 5.
 131  Numd = Hgt*Float(Numd)/Ds  + 0.5
      Iii = Hgt/Ds + 0.5
      If (Cval .Lt. 0.) Numd = Numd + Iii + 1
      T1 = Spi
      If (Spj .Lt. Spi) T1 = Spj
      Numd = Float(Numd)/(2.0*T1) + 0.5
      Mumd = Numd
      Dnum = Dnum + 2.
      If (Cvar .Lt. 0.) Dnum = Dnum + 1.
      Dnum = Dnum*Hgt
      I11 = 1
      Go To 9
 10   I3 = I3 + 1
      Go To 11
 9    Jside = 1
      J4 = 1
      Go To 13
 20   J4 = J4 + 1
      If (J4 .Eq. H) Go To 7
 13   Jdet = 1
      E0 = -1
      F0 = J4 - 2
      If (Cval .Gt. Z(E0+2,F0+2) .And. Cval .Gt. Z(E0+2,F0+3)) Go To 20
      If (Cval .Lt. Z(E0+2,F0+2) .And. Cval .Lt. Z(E0+2,F0+3)) Go To 20
      Go To 26
 56   X = F0 + 2
      Do 59 I = 1,Jk
      Xt(I) = X
 59   X = X + Ds
      Jdet = 0
      A1 = 1

      Do 61 J = 2,Jk
         If (Cval .Lt. Py(1,J) .And. Cval .Lt. Py(1,J-1)) Go To 61
         If (Cval .Gt. Py(1,J) .And. Cval .Gt. Py(1,J-1)) Go To 61
         Fs=(Cval-Py(1,J))/(Py(1,J)-Py(1,J-1))*(Xt(J)-Xt(J-1))+Xt(J)
         Es = 1.
         B1 = J - 1
         Side = 1
         Go To 78
 61   Continue

      Go To 20
 7    Jside = 3
      J4 = 1
      Go To 104
 105  J4 = J4 + 1
      If (J4 .Eq. H) Go To 102
 104  Jdet = 3
      E0 = L - 3
      F0 = J4 - 2
      If (Cval .Gt. Z(E0+3,F0+2) .And. Cval .Gt. Z(E0+3,F0+3)) Go To 105
      If (Cval .Lt. Z(E0+3,F0+2) .And. Cval .Lt. Z(E0+3,F0+3)) Go To 105
      Go To 26
 106  X = F0 + 2

      Do I = 1, Jk
         Xt(I) = X
         X = X + Ds
      End Do

      Jdet = 0
      A1 = Jk - 1

      Do 109 J = 2,Jk
         If (Cval .Lt. Py(Jk,J) .And. Cval .Lt. Py(Jk,J-1)) Go To 109
         If (Cval .Gt. Py(Jk,J) .And. Cval .Gt. Py(Jk,J-1)) Go To 109
         Fs=(Cval-Py(Jk,J))/(Py(Jk,J)-Py(Jk,J-1))*(Xt(J)-Xt(J-1))+Xt(J)
         Es = L
         B1 = J - 1
         Side = 3
         Go To 78
  109 Continue

      Go To 105
 102  Jside = 2
      I4 = 1
      Go To 111
 112  I4 = I4 + 1
      If (I4 .Eq. L) Go To 110
 111  Jdet = 2
      E0 = I4 - 2
      F0 = H - 3
      If (Cval .Gt. Z(E0+2,F0+3) .And. Cval .Gt. Z(E0+3,F0+3)) Go To 112
      If (Cval .Lt. Z(E0+2,F0+3) .And. Cval .Lt. Z(E0+3,F0+3)) Go To 112
      Go To 26
 114  X = E0 + 2

      Do 115 I = 1,Jk
         Xt(I) = X
 115  X = X + Ds

      Jdet = 0
      B1 = Jk - 1

      Do 116 J = 2,Jk
         If (Cval .Lt. Py(J,Jk) .And. Cval .Lt. Py(J-1,Jk)) Go To 116
         If (Cval .Gt. Py(J,Jk) .And. Cval .Gt. Py(J-1,Jk)) Go To 116
         Es=(Cval-Py(J,Jk))/(Py(J,Jk)-Py(J-1,Jk))*(Xt(J)- Xt(J-1))+Xt(J)
         Fs = H
         A1 = J - 1
         Side = 2
         Go To 78
 116  Continue

      Go To 112
 110  Jside = 4
      I4 = 2
      J4 = 1
 12   If (Cval .Gt. Z(I4,J4) .And. Cval .Lt. Z(I4-1,J4)) Go To 16
      If (Cval .Lt. Z(I4,J4) .And. Cval .Gt. Z(I4-1,J4)) Go To 16
      Go To 17
 16   Jdet = 4
      E0 = I4 - 3
      F0 = J4 - 2
      Go To 26
 75   X = E0 + 2

      Do 76 I = 1,Jk
         Xt(I) = X
 76      X = X + Ds

      Jdet = 0
      B1 = 1

      Do 77 I = 2,Jk
         If ((Cval .Lt. Py(I,1) .And. Cval .Lt. Py(I-1,1)) .Or.
     .        (Cval .Gt. Py(I,1) .And. Cval .Gt. Py(I-1,1))) Go To 77
         Fs = J4
         Es=(Cval-Py(I,1))/(Py(I,1)-Py(I-1,1))*(Xt(I)-Xt(I-1))+Xt(I)
         A1 = I - 1
         Side = 4
         Go To 78
 77   Continue

 78   Ipen = 2
      I12 = I11 - 2
      If (I12 .Lt. 1) Go To 54

      Do 120 I = 1,I12
         If (Es + Adj .Gt. Est(I) .And. Es - Adj .Lt. Est(I) .And.
     .        Fs + Adj .Gt. Fst(I) .And. Fs - Adj .Lt. Fst(I)) Go To 17
 120  Continue

 54   Kount = 1
      If (Idsh .Eq. 0) Call Plot(Fs*Spj, Es*Spi, 3)
      If (Idsh .Eq. 1) Call Dash (Fs*Spj, Es*Spi, 3, Pn, Qn)
      Fst(I11) = Fs
      Est(I11) = Es
      I11 = I11 + 1
      Ek = Es
      Fk = Fs
      Go To 400
 17   If (Jside .Eq. 1) Go To 20
      If (Jside .Eq. 2) Go To 112
      If (Jside .Eq. 3) Go To 105
      I4 = I4 + 1
      If (I4 .Ne. L + 1) Go To 15
      I4 = 2
      J4 = J4 + 1
 15   If (J4 .Eq. H) Go To 10
      Go To 12
 21   If (Side .Eq. 1) E0 = E0 + 1
      If (Side .Eq. 2) F0 = F0 - 1
      If (Side .Eq. 3) E0 = E0 - 1
      If (Side .Eq. 4) F0 = F0 + 1

 26   Do 27 I = 1,4
         Iii = E0 + I
         Do 27 J = 1,4
            Jjj = F0 + J
            If (Iii .Ne. 0 .Or. Jjj .Ne. 0) Go To 71
            Trip = 1
            Go To 27
 71         If (Iii .Ne. 0 .Or. Jjj .Ne. H + 1) Go To 72
            Trip = 2
            Go To 27
 72         If (Iii .Ne. L + 1 .Or. Jjj .Ne. H + 1) Go To 73
            Trip = 3
            Go To 27
 73         If (Iii .Ne. L + 1 .Or. Jjj .Ne. 0) Go To 74
            Trip = 4
            Go To 27
 74         If (Iii .Eq. 0) Go To 28
            If (Iii .Eq. L + 1) Go To 29
            If (Jjj .Eq. 0) Go To 30
            If (Jjj .Eq. H + 1) Go To 31
            F(I,J) = Z(Iii,Jjj)
            Go To 27
 28         F(I,J) = Z(3,Jjj) + 3.*(Z(1,Jjj) - Z(2,Jjj))
            Go To 27
 29         F(I,J) = Z(L-2,Jjj) - 3.*(Z(L-1,Jjj) - Z(L,Jjj))
            Go To 27
 30         F(I,J) = Z(Iii,3) + 3.*(Z(Iii,1) - Z(Iii,2))
            Go To 27
 31         F(I,J) = Z(Iii,H-2) - 3.*(Z(Iii,H-1) - Z(Iii,H))
 27   Continue

      If (Trip .Eq. 1) F(1,1) = F(1,4) + 3.*(F(1,2) - F(1,3))
      If (Trip .Eq. 2) F(1,4) = F(1,1) + 3.*(F(1,3) - F(1,2))
      If (Trip .Eq. 3) F(4,1) = F(1,1) + 3.*(F(3,1) - F(2,1))
      If (Trip .Eq. 4) F(4,4) = F(1,4) + 3.*(F(3,4) - F(2,4))
      Trip = 0
      X = 0.

      Do 33 J = 1,Jk
         T1 = (X - 1.)*(X - 2.)
         Tw(J) = X*T1*(-.1666666667)
         Tx(J) = (X + 1.)*T1*.5
         T1 = (X + 1.)*X
         Ty(J) = T1*(X - 2.)*(-.5)
         Tz(J) = T1*(X - 1.)*.1666666667
         X = X + Ds

         Do 33 I = 1,4
 33   Px(I,J) = Tw(J)*F(1,I) + Tx(J)*F(2,I) + Ty(J)*F(3,I) +Tz(J)*F(4,I)

      Do 34 J = 1,Jk
         Do 34 I = 1,Jk
            Py(I,J) = Tw(J)*Px(1,I) + Tx(J)*Px(2,I) + Ty(J)*Px(3,I)
     .           + Tz(J)*Px(4,I)
            If (Py(I,J) .Eq. Cval) Py(I,J) = Py(I,J) + Ad
 34   Continue

      If (Jdet .Eq. 1) Go To 56
      If (Jdet .Eq. 2) Go To 114
      If (Jdet .Eq. 3) Go To 106
      If (Jdet .Eq. 4) Go To 75
      If (Side .Eq. 1) A1 = 1
      If (Side .Eq. 2) B1 = Jk - 1
      If (Side .Eq. 3) A1 = Jk - 1
      If (Side .Eq. 4) B1 = 1
 400  If (A1 .Lt. 1 .Or. B1 .Lt. 1 .Or. A1 .Ge. Jk .Or. B1 .Ge.Jk)Goto17
      Ikount = 0
      Space = 0.
      If (Icntst .Eq. 0) Go To 40
      If (Ctest .Lt. 1.) Go To 40
      Space = 1.
 40   If (Ikount .Eq. 1) Kount = 1
      Nside = 0
      Go To (44, 43, 46, 45), Side
 41   Go To (45, 44, 43, 46), Side
 42   Go To (43, 46, 45, 44), Side
 43   X1 = A1
      Y1 = B1
      X2 = A1 + 1
      Y2 = B1
      Nside = Nside + 1
      If (Py(X1,Y1) .Eq. Py(X2,Y2)) Go To (41, 42, 47), Nside
      X0 = Float(X2) + (Cval - Py(X2,Y2))*Float(X1 - X2)/(Py(X1,Y1)
     . - Py(X2,Y2))
      If (Float(X1) .Lt. X0 .And. Float(X2) .Gt. X0) Go To 443
      If (Nside - 2) 41, 42, 17
 44   X1 = A1 + 1
      Y1 = B1
      X2 = A1 + 1
      Y2 = B1 + 1
      Nside = Nside + 1
      If(Py(X1,Y1).Eq.Py(X2,Y2)) Go To(41,42,47),Nside
      X0=Float(Y2)+(Cval-Py(X2,Y2))*Float(Y1-Y2)/(Py(X1,Y1)-Py(X2,Y2))
      If (Float(Y1) .Lt. X0 .And. Float(Y2) .Gt. X0) Go To 444
      If (Nside - 2) 41, 42, 17
 45   X1 = A1
      Y1 = B1 + 1
      X2 = A1 + 1
      Y2 = B1 + 1
      Nside = Nside + 1
      If (Py(X1,Y1) .Eq. Py(X2,Y2)) Go To (41, 42, 47), Nside
      X0=Float(X2)+(Cval-Py(X2,Y2))*Float(X1-X2)/(Py(X1,Y1)-Py(X2,Y2))
      If (Float(X1) .Lt. X0 .And. Float(X2) .Gt. X0) Go To 445
      If (Nside - 2) 41, 42, 17
 46   X1 = A1
      Y1 = B1
      X2 = A1
      Y2 = B1 + 1
      Nside = Nside + 1
      If (Py(X1,Y1) .Eq. Py(X2,Y2)) Go To (41, 42, 47), Nside
      X0=Float(Y2)+(Cval-Py(X2,Y2))*Float(Y1-Y2)/(Py(X1,Y1)-Py(X2,Y2))
      If (Float(Y1) .Lt. X0 .And. Float(Y2) .Gt. X0) Go To 6527
      If (Nside - 2) 41, 42, 17
   47 Write (7,6525)
 6525 Format (  ' Abnormal Exit of Clot Through format 47; Z '         !RAH
     1 ,'may be altered!!')                                             !RAH
      Write (7,*)  Cval, Es, Fs, Ek, Fk, E0, F0, A1, B1, Side, Jside,
     #             Py(x1,y1), Py(x2,y2), x1, y1, x2, y2
      Call NewPn(NPen(3)) ! pen to use on exit from this subroutine     !RAH
      Return
 6527 Es = Float(E0 + 2) + Float(X1 - 1)*Ds
      Fs = Float(F0 + 2) + (X0 - 1.)*Ds
      Aside = 1
      If (Space .Gt. 0.) Go To 450
 446  Call Pltsg
      Side = 3
      If (A1 .Le. 1) Go To 48
      A1 = A1 - 1
      Go To 40
 443  Es = Float(E0 + 2) + (X0 - 1.)*Ds
      Fs = Float(F0 + 2) + Float(Y1 - 1)*Ds
      Aside = 4
      If (Space .Gt. 0.) Go To 452
 447  Call Pltsg
      Side = 2
      If (B1 .Le. 1) Go To 48
      B1 = B1 - 1
      Go To 40
 444  Es = Float(E0 + 2) + Float(X1 - 1)*Ds
      Fs = Float(F0 + 2) + (X0 - 1.)*Ds
      Aside = 3
      If (Space .Gt. 0.) Go To 454
 448  Call Pltsg
      Side = 1
      A1 = A1 + 1
      If (A1 .Ge. Jk) Go To 48
      Go To 40
 445  Es = Float(E0 + 2) + (X0 - 1.)*Ds
      Fs = Float(F0 + 2) + Float(Y1 - 1)*Ds
      Aside = 2
      If (Space .Gt. 0.) Go To 456
      Call Pltsg
 449  Side = 4
      B1 = B1 + 1
      If (B1 .Lt. Jk) Go To 40
 48   If (Kount .Gt. Jo) Kount = 1
      Iii = Fs - Adj
      If (Iii .Le. 0) Go To 14
      Iii = Es - Adj
      If (Iii .Le. 0) Go To 14
      Iii = Fs + Adj
      If (Iii .Ge. H) Go To 14
      Iii = Es + Adj
      If (Iii .Ge. L) Go To 14
      I12 = I11 - 3
      If (I12 .Lt. 1) Go To 52
      Do 51 I = 1,I12
      If (Es + Adj .Gt. Est(I) .And. Es - Adj .Lt. Est(I) .And.
     . Fs + Adj .Gt. Fst(I) .And. Fs - Adj .Lt. Fst(I)) Go To 17
 51   Continue
 52   Est(I11) = Es
      Fst(I11) = Fs
      I11 = I11 + 1
      If (I11 .Gt. Len) Go To 5003
      Go To 21
 14   Est(I11) = Es
      Fst(I11) = Fs
      I11 = I11 + 1
      Est(I11) = Es
      Fst(I11) = Fs
      I11 = I11 + 1
      Est(I11) = Es
      Fst(I11) = Fs
      I11 = I11 + 1
      Go To 17
C     Section Does Spacing Test For Intermediate Contours
 450  If (Side .Eq. 2) Go To 475
      If (Side .Eq. 3) Go To 460
      If (Side .Eq. 4) Go To 470
 452  If (Side .Eq. 1) Go To 470
      If (Side .Eq. 2) Go To 465
      If (Side .Eq. 3) Go To 480
 454  If (Side .Eq. 1) Go To 460
      If (Side .Eq. 2) Go To 485
      If (Side .Eq. 4) Go To 480
 456  If (Side .Eq. 1) Go To 475
      If (Side .Eq. 3) Go To 485
      If (Side .Eq. 4) Go To 465
 460  Dcon = Abs(Cval - Py(A1,B1))
      If (Side .Eq. 1) Go To 461
      Xj = Fs - (Float(F0 + 2) + Float(B1 - 1)*Ds)
      Go To 462
 461  Xj = Fk - (Float(F0 + 2) + Float(B1 - 1)*Ds)
 462  T1 = Abs(Fs - Fk)
      If (Xj .Eq. 0.0) Xj=.000001
      Consp = Dcon*Sqrt(1. + T1)/Xj
      Go To 490
 465  Dcon = Abs(Cval - Py(A1,B1))
      If (Side .Eq. 4) Go To 466
      Xi = Es - (Float(E0 + 2) + Float(A1 - 1)*Ds)
      Go To 467
 466  Xi = Ek - (Float(E0 + 2) + Float(A1 - 1)*Ds)
 467  T1 = Abs(Es - Ek)
      If (Xi .Eq. 0.0) Xi=.000001
      Consp = Dcon*Sqrt(1. + T1)/Xi
      Go To 490
 470  Dcon = Abs(Cval - Py(A1,B1))
      If (Side .Eq. 4) Go To 471
      Xi = Es - (Float(E0 + 2) + Float(A1 - 1)*Ds)
      Xj = Fk - (Float(F0 + 2) + Float(B1 - 1)*Ds)
      Go To 487
 471  Xi = Ek - (Float(E0 + 2) + Float(A1 - 1)*Ds)
      Xj = Fs - (Float(F0 + 2) + Float(B1 - 1)*Ds)
      Go To 487
 475  Dcon = Abs(Cval - Py(A1,B1+1))
      If (Side .Eq. 2) Go To 476
      Xi = Es - (Float(E0 + 2) + Float(A1 - 1)*Ds)
      Xj = Float(F0 + 2) + Float(B1)*Ds - Fk
      Go To 487
 476  Xi = Ek - (Float(E0 + 2) + Float(A1 - 1)*Ds)
      Xj = Float(F0 + 2) + Float(B1)*Ds - Fs
      Go To 487
 480  Dcon = Abs(Cval - Py(A1+1,B1))
      If (Side .Eq. 4) Go To 481
      Xi = Float(E0 + 2) + Float(A1)*Ds - Es
      Xj = Fk - (Float(F0 + 2) + Float(B1 - 1)*Ds)
      Go To 487
 481  Xi = Float(E0 + 2) + Float(A1)*Ds - Ek
      Xj = Fs - (Float(F0 + 2) + Float(B1 - 1)*Ds)
      Go To 487
 485  Dcon = Abs(Cval - Py(A1+1,B1+1))
      If (Side .Eq. 2) Go To 486
      Xi = Float(E0 + 2) + Float(A1)*Ds - Es
      Xj = Float(F0 + 2) + Float(B1)*Ds - Fk
      Go To 487
 486  Xi = Float(E0 + 2) + Float(A1)*Ds - Ek
      Xj = Float(F0 + 2) + Float(B1)*Ds - Fs
  487 If (Xi .Eq. 0.0) Xi=.000001
      If (Xj .Eq. 0.0) Xj=.000001
      Consp = Dcon*Sqrt(Xi*Xi + Xj*Xj)/(Xi*Xj)
 490  Ipen = 2
c     If (Consp .Gt. Sptest) Ipen = 3
      If (Ipen .Eq. 3) Ikount = 1
      If (Kount .Gt. 10 .And. Kount .Lt. Mumd + 11 .And.((Klab .Eq. 1)
     >  .or. ((Klab .Eq. 3).and.(Idsh.eq.0))
     >  .or. ((Klab .Eq. 2).and.(Idsh.eq.1))))
     1 Ipen=3
      If (Ikount .Eq. 1) Kount = 1
      Space = 0.
      Go To (446, 449, 448, 447), Aside
 5003 Write (7,5004)
 5004 Format (47H0increase Est And Fst And If(I11...  Statements)
 5000 Do 233 I=1,L
      Do 233 J=1,H
 233  Z(I,J) = (Z(I,J) + Cva1)/Sf
      Call NewPn(NPen(3)) ! pen to use on exit from this subroutine     !RAH
      Return
  140 Write (7,141) Cval
 141  Format (51H0more Than 30 Contours.Return To Main.Last Cval Was,
     1F10.2)
      Go To 5000
      End
