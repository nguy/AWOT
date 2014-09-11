      SUBROUTINE SMOTH (YIN,YOUT,K,NWI,Bad)

C  Subroutine to smooth data using Bartlett filter with triangular
C  weighting function with up to 126 points.
C  YIN=Input array
C  YOUT=Output array (can be the same as YIN)
C  K=# of data points
C  NWI=filter width
C  W=weight function
C  SCRACH=temperary function
 
      Dimension YIN(K),YOUT(K)
      DIMENSION W(251),SCRACH(126)

      DATA MAXN /125/

      MAXPTS = 2 * MAXN + 1
      NW = NWI
      IF(NW .GT. MAXPTS) RETURN
      N = NW / 2
      IF(N*2 .EQ. NW) NW = NW + 1
      NP1 = N + 1
      NWP1 = NW + 1

C     CALCULATE BARTLETT WEIGHTS

      SQNP1 = NP1 * NP1

      DO I = 1, N
         W(I) = I / SQNP1
         J = NWP1 - I
         W(J) = W(I)
      End Do

      W(NP1) = 1. / FLOAT(NP1)
 
C     FILTER DATA
 
      IF (K .EQ. 0) RETURN

      DO I = 1, N
         YOUT(I) = Bad
         YOUT(K+1-I) = Bad
         SCRACH(I) = YIN(I)
      End Do
  
      IN    = N
      KMN   = K - N
  
      DO I = NP1, KMN
         IN    = IN + 1
         IF(IN .GT. 126) IN = 1
         SCRACH(IN) = YIN(I)
         IW    = IN
         KFLAG = 0
         IF (YIN(I) .eq. Bad) KFLAG = 1
         SUM   = YIN(I) * W(NP1)

         DO IA = 1,N
            IW    = IW - 1
            IF (IW .LT. 1) IW = 126
            IF (SCRACH(IW).eq. Bad .OR. YIN(I+IA) .eq. Bad) KFLAG = 1
            SUM   = SUM + (SCRACH(IW) + YIN(I+IA)) * W(NP1-IA)
         End Do

         IF (KFLAG .EQ. 1) SUM=Bad
         YOUT(I)= SUM
      End Do

      RETURN
      END
