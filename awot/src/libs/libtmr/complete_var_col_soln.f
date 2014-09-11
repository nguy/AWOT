      SUBROUTINE COMPLETE_VAR_COL_SOLN(Z1,TV1,Z2,TV2,Z3,RHO3,N,Z,DIV,
     $SD_DIV,WP,SD_WP,F,SD_F,W_CONSTRAINT,BADDATA,DIV_ADJ,WP_ADJ,F_ADJ,
     $W,SUCCESS)

C  Thomas Matejka and Diana Bartels NOAA/NSSL 23 September 1994

C  This program calculates variationally adjusted horizontal divergence,
C  variationally adjusted vertical hydrometeor velocity, variationally
C  adjusted hydrometeor terminal fall speed, and vertical air velocity
C  in a column from measurements of the horizontal divergence, vertical
C  hydrometeor velocity, and hydrometeor terminal fall speed at levels
C  in the column.  The solution is found only at and beween the highest
C  set of consecutive levels having complete data, but only if the
C  highest and lowest levels in this set are not identical.

C  The solution is formulated using six conditions: (1) The adjusted
C  horizontal divergences should be as close as possible to the input
C  horizontal divergences.  (2) The adjusted vertical hydrometeor
C  velocities should be as close as possible to the input vertical
C  hydrometeor velocities.  (3) The adjusted hydrometeor terminal fall
C  speeds should be as close as possible to the input hydrometeor
C  terminal fall speeds.  (4) The adjusted horizontal divergences and
C  the vertical air velocities should satisfy anelastic continuity
C  exactly.  (5) The adjusted vertical hydrometeor velocity, the
C  adjusted hydrometeor terminal fall speed, and the vertical air
C  velocity should be exactly related.  (6) The vertical air velocity
C  should equal any input constraints exactly.  All adjustments are
C  performed in terms of mass divergence (not divergence) and vertical
C  air mass flux (not vertical velocity).

C  Input:

C  Z1 is a real variable that specifies an altitude (m).

C  TV1 is a real variable that specifies the virtual temperature (K) at
C  Z1.

C  Z2 is a real variable that specifies an altitude (m).  Z2 must not be
C  equal to Z1.

C  TV2 is a real variable that specifies the virtual temperature (K) at
C  Z2.

C  Z3 is a real variable that specifies an altitude (m).  Z3 may be
C  equal to Z1 or Z2.
 
C  RHO3 is a real variable that specifies the density (kg m**-3) at Z3.
C  Z1, TV1, Z2, TV2, Z3, and RHO3 define a density profile in a
C  constant-lapse-rate atmosphere.

C  N is an integer variable that specifies the number of layers in the
C  column and one less than the number of levels in the column.

C  Z is a one-dimensional real array with elements indexed from 0 to N.
C  Z(I) specifies the altitude (m) of the Ith level above the bottom of
C  the column.

C  DIV is a one-dimensional real array with elements indexed from 0 to
C  N.  DIV(I) specifies the horizontal divergence (s**-1) at the Ith
C  level above the bottom of the column.  If it is missing, it should
C  equal BADDATA.

C  SD_DIV is a one-dimensional real array with elements indexed from 0
C  to N.  SD_DIV(I) specifies the standard deviation of DIV(I) (s**-1).
C  If it is missing, it should equal BADDATA.  SD_DIV(I) must not be
C  missing for DIV(I) to be usable.

C  WP is a one-dimensional real array with elements indexed from 0 to N.
C  WP(I) specifies the mean echo-power-weighted vertical hydrometeor
C  velocity (m s**-1) at the Ith level above the bottom of the column.
C  If it is missing, it should equal BADDATA.  Vertical hydrometeor
C  velocity is positive for upward motion.

C  SD_WP is a one-dimensional real array with elements indexed from 0 to
C  N.  SD_WP(I) specifies the standard deviation of WP(I) (m s**-1).  If
C  it is missing, it should equal BADDATA.  SD_WP(I) must not be missing
C  for WP(I) to be usable.

C  F is a one-dimensional real array with elements indexed from 0 to N.
C  F(I) specifies the mean echo-power-weighted hydrometeor terminal fall
C  speed (m s**-1) at the Ith level above the bottom of the column.  If
C  it is missing, it should equal BADDATA.  Hydrometeor terminal fall
C  speed is greater than or equal to zero.

C  SD_F is a one-dimensional real array with elements indexed from 0 to
C  N.  SD_F(I) specifies the standard deviation of F(I) (m s**-1).  If
C  it is missing, it should equal BADDATA.  SD_F(I) must not be missing
C  for F(I) to be usable.

C  W_CONSTRAINT is a one-dimensional real array with elements indexed
C  from 0 to N.  W_CONSTRAINT(I) specifies a strong vertical air
C  velocity constraint (m s**-1) at the Ith level above the bottom of
C  the column.  If it is missing, it should equal BADDATA.  Vertical air
C  velocity boundary conditions, if desired, are specified in this way.
C  Vertical air velocity is positive for upward motion.

C  BADDATA is a real variable that indicates a missing value as
C  described.

C  Output:

C  DIV_ADJ is a one-dimensional real array with elements indexed from 0
C  to N.  DIV_ADJ(I) returns the adjusted horizontal divergence(s**-1)
C  at the Ith level above the bottom of the column.  If it is not
C  calculated, it is returned as BADDATA.

C  WP_ADJ is a one-dimensional real array with elements indexed from 0
C  to N.  WP_ADJ(I) returns the adjusted mean echo-power-weighted
C  vertical hydrometeor velocity (m s**-1) at the Ith level above the
C  bottom of the column.  If it is not calculated, it is returned as
C  BADDATA.

C  F_ADJ is a one-dimensional real array with elements indexed from 0 to
C  N.  F_ADJ(I) returns the adjusted mean echo-power-weighted terminal
C  fall speed (m s**-1) at the Ith level above the bottom of the column.
C  If it is not calculated, it is returned as BADDATA.

C  W is a one-dimensional real array with elements indexed from 0 to N.
C  W(I) returns the vertical air velocity (m s**-1) at the Ith level
C  above the bottom of the column.  If it is not calculated, it is
C  returned as BADDATA.

C  SUCCESS is a logical variable that is returned .TRUE. if and only if
C  a solution of any sort was able to be found.

      IMPLICIT NONE
      REAL DENCLR,AV_DENCLR
      LOGICAL SUCCESS
      INTEGER N,K,KMIN,KMAX,I
      REAL Z1,TV1,Z2,TV2,Z3,RHO3,BADDATA
      REAL DELTAZ_ABV(0:N-1),RHO_BAR_ABV(0:N-1)
      REAL DELTAZ_BLW(N),RHO_BAR_BLW(N)
      REAL Z(0:N),DIV(0:N),WP(0:N),F(0:N),SD_DIV(0:N),SD_WP(0:N),
     $SD_F(0:N),W_CONSTRAINT(0:N),DIV_ADJ(0:N),WP_ADJ(0:N),F_ADJ(0:N),
     $W(0:N),RHO(0:N),DELTAZ(0:N),RHO_BAR(0:N)
      DOUBLEPRECISION RHS(N),SOLUTION(N),LAMBDA(N)
      DOUBLEPRECISION PHI(0:N),CAP_LAMBDA(0:N)
      DOUBLEPRECISION A(15,0:N)
      DOUBLEPRECISION B(7,0:N)
      DOUBLEPRECISION C(4,N)
      DOUBLEPRECISION COEFF_MATRIX(N,-1:1)

C  Initialize the output as missing.
      DO K=0,N
         DIV_ADJ(K)=BADDATA
         WP_ADJ(K)=BADDATA
         F_ADJ(K)=BADDATA
         W(K)=BADDATA
      ENDDO

C  Find the top of the data column.
      DO K=N,1,-1
         IF(DIV(K).NE.BADDATA.AND.
     $   SD_DIV(K).NE.BADDATA.AND.
     $   WP(K).NE.BADDATA.AND.
     $   SD_WP(K).NE.BADDATA.AND.
     $   F(K).NE.BADDATA.AND.
     $   SD_F(K).NE.BADDATA)THEN
            KMAX=K
            GOTO 1
         ENDIF
      ENDDO
      SUCCESS=.FALSE.
      RETURN
1     CONTINUE

C  Find the bottom of the data column.
      DO K=KMAX-1,0,-1
         IF(DIV(K).EQ.BADDATA.OR.
     $   SD_DIV(K).EQ.BADDATA.OR.
     $   WP(K).EQ.BADDATA.OR.
     $   SD_WP(K).EQ.BADDATA.OR.
     $   F(K).EQ.BADDATA.OR.
     $   SD_F(K).EQ.BADDATA)THEN
            KMIN=K+1
            GOTO 2         
         ENDIF
      ENDDO
      KMIN=0
2     CONTINUE

C  Check whether there is at least one layer.
      IF(KMAX.LE.KMIN)THEN
         SUCCESS=.FALSE.
         RETURN
      ENDIF

C  Set PHI to 0. or 1. according to whether vertical air velocity
C  constraints exist or are missing.
      DO K=KMIN,KMAX
         IF(W_CONSTRAINT(K).EQ.BADDATA)THEN
            PHI(K)=0.D0
         ELSE
            PHI(K)=1.D0
         ENDIF
      ENDDO

C  Calculate the depth of the bottom half-layer, the top half-layer, and
C  the middle layers centered on levels.
      DELTAZ(KMIN)=(Z(KMIN+1)-Z(KMIN))/2.
      DELTAZ(KMAX)=(Z(KMAX)-Z(KMAX-1))/2.
      IF(KMAX.GE.KMIN+2)THEN
         DO K=KMIN+1,KMAX-1
            DELTAZ(K)=(Z(K+1)-Z(K-1))/2.
         ENDDO
      ENDIF

C  Calculate the depth of the half layers below levels.
      DO K=KMIN+1,KMAX
         DELTAZ_BLW(K)=(Z(K)-Z(K-1))/2.
      ENDDO

C  Calculate the depth of the half layers above levels.
      DO K=KMIN,KMAX-1
         DELTAZ_ABV(K)=(Z(K+1)-Z(K))/2.
      ENDDO

C  Calculate the density at levels.
      DO K=KMIN,KMAX
         RHO(K)=DENCLR(Z1,TV1,Z2,TV2,Z3,RHO3,Z(K))
      ENDDO

C  Calculate the mean density in the bottom half-layer, the top
C  half-layer, and the middle layers centered on levels.
      RHO_BAR(KMIN)=AV_DENCLR(Z1,TV1,Z2,TV2,Z3,RHO3,Z(KMIN),
     $(Z(KMIN)+Z(KMIN+1))/2.)
      RHO_BAR(KMAX)=AV_DENCLR(Z1,TV1,Z2,TV2,Z3,RHO3,
     $(Z(KMAX-1)+Z(KMAX))/2.,Z(KMAX))
      IF(KMAX.GE.KMIN+2)THEN
         DO K=KMIN+1,KMAX-1
            RHO_BAR(K)=AV_DENCLR(Z1,TV1,Z2,TV2,Z3,RHO3,(Z(K-1)+Z(K))/2.,
     $      (Z(K)+Z(K+1))/2.)
         ENDDO
      ENDIF

C  Calculate the mean density in the half-layers below levels.
      DO K=KMIN+1,KMAX
         RHO_BAR_BLW(K)=AV_DENCLR(Z1,TV1,Z2,TV2,Z3,RHO3,
     $   (Z(K-1)+Z(K))/2.,Z(K))
      ENDDO

C  Calculate the mean density in the half-layers above levels.
      DO K=KMIN,KMAX-1
         RHO_BAR_ABV(K)=AV_DENCLR(Z1,TV1,Z2,TV2,Z3,RHO3,Z(K),
     $   (Z(K)+Z(K+1))/2.)
      ENDDO
          
C  Calculate coefficients in the equations.
      DO K=KMIN,KMAX
         A(1,K)=DBLE(2.*RHO_BAR(K)*DELTAZ(K)/SD_DIV(K)**2)
         IF(K.EQ.KMIN)THEN
            A(2,K)=0.D0
         ELSE
            A(2,K)=DBLE(RHO_BAR_BLW(K)*DELTAZ_BLW(K))
         ENDIF
         IF(K.EQ.KMAX)THEN
            A(3,K)=0.D0
         ELSE
            A(3,K)=DBLE(RHO_BAR_ABV(K)*DELTAZ_ABV(K))
         ENDIF
         A(4,K)=DBLE(2.*RHO_BAR(K)*DELTAZ(K)*DIV(K)/SD_DIV(K)**2)
         A(5,K)=DBLE(2.*RHO_BAR(K)*DELTAZ(K)/SD_WP(K)**2)
         A(6,K)=DBLE(2.*RHO_BAR(K)*DELTAZ(K)*WP(K)/SD_WP(K)**2)
         A(7,K)=DBLE(2.*RHO_BAR(K)*DELTAZ(K)/SD_F(K)**2)
         A(8,K)=DBLE(2.*RHO_BAR(K)*DELTAZ(K)*F(K)/SD_F(K)**2)
         IF(K.EQ.KMIN)THEN
            A(9,K)=0.D0
         ELSEIF(W_CONSTRAINT(K).EQ.BADDATA)THEN
            A(9,K)=DBLE(RHO(K))
         ELSE
            A(9,K)=0.D0
         ENDIF
         IF(K.EQ.KMAX)THEN
            A(10,K)=0.D0
         ELSEIF(W_CONSTRAINT(K).EQ.BADDATA)THEN
            A(10,K)=DBLE(-RHO(K))
         ELSE
            A(10,K)=0.D0
         ENDIF
         IF(K.EQ.KMIN)THEN
            A(11,K)=0.D0
         ELSE
            A(11,K)=DBLE(RHO_BAR_ABV(K-1)*DELTAZ_ABV(K-1))
         ENDIF
         IF(K.EQ.KMIN)THEN
            A(12,K)=0.D0
         ELSE
            A(12,K)=DBLE(RHO_BAR_BLW(K)*DELTAZ_BLW(K))
         ENDIF
         IF(K.EQ.KMIN)THEN
            A(13,K)=0.D0
         ELSE
            A(13,K)=DBLE(-RHO(K-1))
         ENDIF
         IF(K.EQ.KMIN)THEN
            A(14,K)=0.D0
         ELSE
            A(14,K)=DBLE(RHO(K))
         ENDIF
         IF(W_CONSTRAINT(K).NE.BADDATA)THEN
            A(15,K)=DBLE(W_CONSTRAINT(K))
         ELSE
            A(15,K)=0.D0
         ENDIF
      ENDDO
      DO K=KMIN,KMAX
         IF(W_CONSTRAINT(K).NE.BADDATA)THEN
            B(1,K)=(A(6,K)*A(7,K)+A(8,K)*A(5,K)-A(15,K)*A(5,K)*A(7,K))/
     $      (A(7,K)+A(5,K))
         ELSE
            B(1,K)=0.D0
         ENDIF
      ENDDO
      DO K=KMIN+1,KMAX
         B(2,K)=A(4,K-1)*A(11,K)/A(1,K-1)+A(4,K)*A(12,K)/A(1,K)+
     $   A(6,K-1)*A(13,K)/A(5,K-1)+A(8,K-1)*A(13,K)/A(7,K-1)+
     $   A(6,K)*A(14,K)/A(5,K)+A(8,K)*A(14,K)/A(7,K)
         IF(K.EQ.KMIN+1)THEN
            B(3,K)=0.D0
         ELSE
            B(3,K)=-A(2,K-1)*A(11,K)/A(1,K-1)
         ENDIF
         B(4,K)=-A(3,K-1)*A(11,K)/A(1,K-1)-A(2,K)*A(12,K)/A(1,K)
         IF(K.EQ.KMAX)THEN
            B(5,K)=0.D0
         ELSE
            B(5,K)=-A(3,K)*A(12,K)/A(1,K)
         ENDIF
         B(6,K)=-A(13,K)/A(5,K-1)-A(13,K)/A(7,K-1)
         B(7,K)=-A(14,K)/A(5,K)-A(14,K)/A(7,K)
      ENDDO
      DO K=KMIN+1,KMAX
         C(1,K)=-B(2,K)-PHI(K-1)*B(1,K-1)*B(6,K)-PHI(K)*B(1,K)*B(7,K)
         IF(K.EQ.KMIN+1)THEN
            C(2,K)=0.D0
         ELSE
            C(2,K)=B(3,K)+(1.D0-PHI(K-1))*A(9,K-1)*B(6,K)
         ENDIF
         C(3,K)=B(4,K)+(1.D0-PHI(K-1))*A(10,K-1)*B(6,K)+
     $   (1.D0-PHI(K))*A(9,K)*B(7,K)
         IF(K.EQ.KMAX)THEN
            C(4,K)=0.D0
         ELSE
            C(4,K)=B(5,K)+(1.D0-PHI(K))*A(10,K)*B(7,K)
         ENDIF
      ENDDO

C  Solve for LAMBA.
      COEFF_MATRIX(1,-1)=0.
      COEFF_MATRIX(1,0)=C(3,KMIN+1)
      COEFF_MATRIX(1,1)=C(4,KMIN+1)
      COEFF_MATRIX(KMAX-KMIN,-1)=C(2,KMAX)
      COEFF_MATRIX(KMAX-KMIN,0)=C(3,KMAX)
      COEFF_MATRIX(KMAX-KMIN,1)=0.
      IF(KMAX.GE.KMIN+3)THEN
         DO I=2,KMAX-KMIN-1
            COEFF_MATRIX(I,-1)=C(2,KMIN+I)
            COEFF_MATRIX(I,0)=C(3,KMIN+I)
            COEFF_MATRIX(I,1)=C(4,KMIN+I)
         ENDDO
      ENDIF
      DO I=1,KMAX-KMIN
         RHS(I)=C(1,KMIN+I)
      ENDDO
      CALL TRI_DIAGONAL_EFFICIENT_DP(COEFF_MATRIX,N,KMAX-KMIN,RHS,
     $SOLUTION,SUCCESS)
      IF(.NOT.SUCCESS)THEN
         RETURN
      ENDIF
      DO K=KMIN+1,KMAX
         LAMBDA(K)=SOLUTION(K-KMIN)
      ENDDO

C  Solve for CAP_LAMBDA.
      CAP_LAMBDA(KMIN)=PHI(KMIN)*B(1,KMIN)+
     $LAMBDA(KMIN+1)*(1.D0-PHI(KMIN))*A(10,KMIN)
      CAP_LAMBDA(KMAX)=PHI(KMAX)*B(1,KMAX)+
     $LAMBDA(KMAX)*(1.D0-PHI(KMAX))*A(9,KMAX)
      IF(KMAX.GE.KMIN+2)THEN
         DO K=KMIN+1,KMAX-1
            CAP_LAMBDA(K)=PHI(K)*B(1,K)+LAMBDA(K)*(1.D0-PHI(K))*A(9,K)+
     $      LAMBDA(K+1)*(1.D0-PHI(K))*A(10,K)
         ENDDO
      ENDIF

C  Solve for DIV_ADJ, WP_ADJ, F_ADJ, and W.
      DO K=KMIN,KMAX
         DIV_ADJ(K)=A(4,K)/A(1,K)-SNGL(LAMBDA(K))*A(2,K)/A(1,K)-
     $   SNGL(LAMBDA(K+1))*A(3,K)/A(1,K)
         WP_ADJ(K)=A(6,K)/A(5,K)-SNGL(CAP_LAMBDA(K))/A(5,K)
         F_ADJ(K)=A(8,K)/A(7,K)-SNGL(CAP_LAMBDA(K))/A(7,K)
         W(K)=WP_ADJ(K)+F_ADJ(K)
      ENDDO

C  Done.
      RETURN
      END
