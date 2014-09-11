
      Function Vt_Adj (dBZ,Hgt,HSnow,HRain)
  
C This function computes mean terminal velocity from the reflectivity
C according to Paul Willis' 2-parameter gamma distribution and
C a snow relationship developed by Heymsfield (1978).
C Reflectivity, in terms of dBZ, must by passed to this routine.
 
C HEIGHT MUST BE IN KM
  
      Z = dBZ
      If (Z .gt. 63.0) z = 63.0
      ZZ=10.0**(Z/10.0)
  
      H = Hgt
      If (H .gt. 20.0) H = 20.0
  
C density correction term (rhoo/rho)*0.4 [rho(Z)=rhoo exp-(z/H), where
C  H is the scale height = 9.58125 from Gray's inner 2 deg composite]
  
      DCOR=EXP(0.4*H/9.58125)
  
C The snow relationship (Atlas et al., 1973) --- VT=0.817*Z**0.063 (m/s)
C     BFS determined value of .500 empirically Dec1,1991
  
      VTS=-DCOR * (0.500*ZZ**0.063)
  
C The rain relationship --- from Willis analytical-gamma distribution
  
C     TERM1=7.331/ZZ**0.010022
C     TERM2=0.14034*ZZ**0.095238
C     VTR=-DCOR * (5.5011E+09/(TERM1+TERM2)**10.5)
  
C The rain relationship (Joss and Waldvogel,1971) -VT=2.6*Z**.107 (m/s)
  
      VTR=-DCOR * (2.6*ZZ**.107)
  
C test if height is in the transition region between SNOW and RAIN
C  defined as HRain km < H < HSnow km
C  if in the transition region do a linear weight of VTR and VTS
 
      VT_Adj=VTR*(HSnow-H)/0.5 + VTS*(H-HRain)/0.5
      If (H.lt.HRain) VT_Adj=VTR
      If (H.gt.HSnow) VT_Adj=VTS
      Return

      End
