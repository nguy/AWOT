      SUBROUTINE W_ADJUST_LEVELS(Z1,TV1,Z2,TV2,Z3,RHO3,N,Z,DIV,SD_DIV,
     $I_BC1,W_BC1,SD_W_BC1,I_BC2,W_BC2,SD_W_BC2,BADDATA,DIV_ADJ,
     $SD_DIV_ADJ,W,SD_W,SUCCESS)

C  Thomas Matejka NOAA/NSSL 23 September 1994

C  This subroutine calculates variationally adjusted horizontal
C  divergence and its standard deviation and the vertical air velocity
C  and its standard deviation in a column from measurements of the
C  horizontal divergence at levels and vertical air velocity boundary
C  conditions at two levels in the column.  The solution is found only
C  at and between the two boundary levels, but only if the divergence is
C  not missing at or anywhere in between the two boundary levels and
C  only if the two boundary levels are not identical.

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

C  I_BC1 is an integer variable that specifies the level number above
C  the bottom of the column at which the first vertical air velocity
C  boundary condition is specified.

C  W_BC1 is a real variable that specifies the vertical air velocity
C  boundary condition (m s**-1) at the I_BC1th level above the bottom of
C  the column.  If it is missing, it should equal BADDATA.

C  SD_W_BC1 is a real variable that specifies the standard deviation of
C  W_BC1 (m s**-1).  If it is missing, it should equal BADDATA.
C  SD_W_BC1 must not be missing for W_BC1 to be usable.

C  I_BC2 is an integer variable that specifies the level number above
C  the bottom of the column at which the second vertical air velocity
C  boundary condition is specified.

C  W_BC2 is a real variable that specifies the vertical air velocity
C  boundary condition (m s**-1) at the I_BC2th level above the bottom of
C  the column.  If it is missing, it should equal BADDATA.

C  SD_W_BC2 is a real variable that specifies the standard deviation of
C  W_BC2 (m s**-1).  If it is missing, it should equal BADDATA.
C  SD_W_BC2 must not be missing for W_BC2 to be usable.

C  BADDATA is a real variable that indicates a missing value as
C  described.

C  Output:

C  DIV_ADJ is a one-dimensional real array with elements indexed from 0
C  to N.  DIV_ADJ(I) returns the adjusted mean horizontal divergence
C  (s**-1) at the Ith level above the bottom of the column.  If it
C  cannot be calculated, it is returned as BADDATA.

C  SD_DIV_ADJ is a one-dimensional real array with elements indexed from
C  0 to N.  SD_DIV_ADJ(I) returns the standard deviation of DIV_ADJ(I)
C  (s**-1).  If it cannot be calculated, it is returned as BADDATA.

C  W is a one-dimensional real array with elements indexed from 0 to N.
C  W(I) returns the vertical air velocity (m s**-1) at the Ith level
C  above the bottom of the column.  If it cannot be calculated, it is
C  returned as BADDATA.

C  SD_W is a one-dimensional real array with elements indexed from 0 to
C  N.  SD_W(I) returns the standard deviation of W(I) (m s**-1).  If it
C  cannot be calculated, it is returned as BADDATA.

C  SUCCESS is a logical variable that is returned as .TRUE. if and only
C  if a solution was found between the two boundary levels.

      IMPLICIT NONE
      REAL DENCLR,AV_DENCLR
      LOGICAL SUCCESS
      INTEGER I,J,N,I_BC1,I_BC2,I_BCLO,I_BCHI
      REAL Z1,Z2,Z3,TV1,TV2,RHO3,W_BC1,W_BC2,SD_W_BC1,SD_W_BC2,SUM1,
     $SUM2,SUM3,SUM4,SUM5,SUM6,SUM7,SUM8,TERM1,TERM3,TERM6,TERM8,
     $BADDATA,W_BCLO,W_BCHI,SD_W_BCLO,SD_W_BCHI
      REAL RHO_BAR_BELOW(N),DEL_Z_BELOW(N)
      REAL DIV(0:N),SD_DIV(0:N),Z(0:N),DIV_ADJ(0:N),SD_DIV_ADJ(0:N),
     $W(0:N),SD_W(0:N),RHO(0:N),RHO_W(0:N),VAR_RHO_W(0:N),VAR_DIV(0:N),
     $VAR_DIV_ADJ(0:N),A(0:N),B(0:N),C(0:N)
      REAL RHO_BAR_ABOVE(0:N-1),DEL_Z_ABOVE(0:N-1)

C  Initialize the output as missing or as the boundary conditions.
      DO I=0,N
         DIV_ADJ(I)=BADDATA
         SD_DIV_ADJ(I)=BADDATA
         W(I)=BADDATA
         SD_W(I)=BADDATA
      ENDDO
      W(I_BC1)=W_BC1
      SD_W(I_BC1)=SD_W_BC1
      W(I_BC2)=W_BC2
      SD_W(I_BC2)=SD_W_BC2

C  Check whether a vertical air velocity boundary condition is missing.
      IF(W_BC1.EQ.BADDATA.OR.
     $SD_W_BC1.EQ.BADDATA.OR.
     $W_BC2.EQ.BADDATA.OR.
     $SD_W_BC2.EQ.BADDATA.OR.
     $I_BC1.EQ.I_BC2)THEN
         SUCCESS=.FALSE.
         RETURN
      ENDIF

C  Identify the lower and higher boundary levels.
      IF(I_BC1.LT.I_BC2)THEN
         I_BCLO=I_BC1
         W_BCLO=W_BC1
         SD_W_BCLO=SD_W_BC1
         I_BCHI=I_BC2
         W_BCHI=W_BC2
         SD_W_BCHI=SD_W_BC2
      ELSE
         I_BCLO=I_BC2
         W_BCLO=W_BC2
         SD_W_BCLO=SD_W_BC2
         I_BCHI=I_BC1
         W_BCHI=W_BC1
         SD_W_BCHI=SD_W_BC1
      ENDIF

C  Check whether a divergence is missing between the two boundary
C  levels.
      DO I=I_BCLO,I_BCHI
         IF(DIV(I).EQ.BADDATA.OR.
     $   SD_DIV(I).EQ.BADDATA)THEN
            SUCCESS=.FALSE.
            RETURN
         ENDIF
      ENDDO

C  Calculate the density and the variance of the divergence at levels.
      DO I=I_BCLO,I_BCHI
         RHO(I)=DENCLR(Z1,TV1,Z2,TV2,Z3,RHO3,Z(I))
         VAR_DIV(I)=SD_DIV(I)**2
      ENDDO

C  Calculate the depth and the mean density in half-layers below levels.
      DO I=I_BCLO+1,I_BCHI
         DEL_Z_BELOW(I)=(Z(I)-Z(I-1))/2.
         RHO_BAR_BELOW(I)=AV_DENCLR(Z1,TV1,Z2,TV2,Z3,RHO3,
     $   (Z(I-1)+Z(I))/2.,Z(I))
      ENDDO

C  Calculate the depth and the mean density in half-layers above levels.
      DO I=I_BCLO,I_BCHI-1
         DEL_Z_ABOVE(I)=(Z(I+1)-Z(I))/2.
         RHO_BAR_ABOVE(I)=AV_DENCLR(Z1,TV1,Z2,TV2,Z3,RHO3,
     $   Z(I),(Z(I)+Z(I+1))/2.)
      ENDDO

C  Calculate the depth and the mass per horizontal area appropriate for
C  each level.
      A(I_BCLO)=RHO_BAR_ABOVE(I_BCLO)*DEL_Z_ABOVE(I_BCLO)
      IF(I_BCHI.GE.I_BCLO+2)THEN
         DO I=I_BCLO+1,I_BCHI-1
            A(I)=RHO_BAR_ABOVE(I)*DEL_Z_ABOVE(I)+RHO_BAR_BELOW(I)*
     $      DEL_Z_BELOW(I)
         ENDDO
      ENDIF
      A(I_BCHI)=RHO_BAR_BELOW(I_BCHI)*DEL_Z_BELOW(I_BCHI)

C  Calculate the vertical mass flux and its variance at the boundary
C  levels.
      RHO_W(I_BCLO)=RHO(I_BCLO)*W_BCLO
      VAR_RHO_W(I_BCLO)=(RHO(I_BCLO)*SD_W_BCLO)**2
      RHO_W(I_BCHI)=RHO(I_BCHI)*W_BCHI
      VAR_RHO_W(I_BCHI)=(RHO(I_BCHI)*SD_W_BCHI)**2

C  Calculate some sums and terms.
      SUM1=0.
      SUM2=0.
      SUM3=0.
      DO I=I_BCLO,I_BCHI
         SUM1=SUM1+DIV(I)*A(I)
         SUM2=SUM2+VAR_DIV(I)*A(I)
         SUM3=SUM3+VAR_DIV(I)*A(I)**2
      ENDDO
      TERM1=RHO_W(I_BCHI)-RHO_W(I_BCLO)+SUM1
      TERM3=VAR_RHO_W(I_BCHI)+VAR_RHO_W(I_BCLO)+SUM3

C  Calculate the adjusted divergence and its variance at levels.
      DO I=I_BCLO,I_BCHI
         DIV_ADJ(I)=DIV(I)-VAR_DIV(I)*TERM1/SUM2
         VAR_DIV_ADJ(I)=(1.-VAR_DIV(I)*A(I)/SUM2)**2*VAR_DIV(I)+
     $   (VAR_DIV(I)/SUM2)**2*(TERM3-VAR_DIV(I)*A(I)**2)
      ENDDO

C  Calculate the standard deviation of the adjusted divergence at
C  levels.
      DO I=I_BCLO,I_BCHI
         SD_DIV_ADJ(I)=SQRT(VAR_DIV_ADJ(I))
      ENDDO

C  Calculate the vertical mass flux and its variance at the levels
C  between the boundary levels.
      IF(I_BCHI.GE.I_BCLO+2)THEN
         DO I=I_BCLO+1,I_BCHI-1
            B(I_BCLO)=RHO_BAR_ABOVE(I_BCLO)*DEL_Z_ABOVE(I_BCLO)
            IF(I.GE.I_BCLO+2)THEN
               DO J=I_BCLO+1,I-1
                  B(J)=RHO_BAR_ABOVE(J)*DEL_Z_ABOVE(J)+RHO_BAR_BELOW(J)*
     $            DEL_Z_BELOW(J)
               ENDDO
            ENDIF
            B(I)=RHO_BAR_BELOW(I)*DEL_Z_BELOW(I)
            C(I)=RHO_BAR_ABOVE(I)*DEL_Z_ABOVE(I)
            IF(I.LE.I_BCHI-2)THEN
               DO J=I+1,I_BCHI-1
                  C(J)=RHO_BAR_ABOVE(J)*DEL_Z_ABOVE(J)+RHO_BAR_BELOW(J)*
     $            DEL_Z_BELOW(J)
               ENDDO
            ENDIF
            C(I_BCHI)=RHO_BAR_BELOW(I_BCHI)*DEL_Z_BELOW(I_BCHI)
            SUM4=0.
            SUM5=0.
            DO J=I_BCLO,I
               SUM4=SUM4+DIV_ADJ(J)*B(J)
               SUM5=SUM5+VAR_DIV(J)*B(J)
            ENDDO
            SUM6=0.
            DO J=I_BCLO,I-1
               SUM6=SUM6+VAR_DIV(J)*B(J)**2
            ENDDO
            TERM6=VAR_RHO_W(I_BCLO)+SUM6
            SUM7=0.
            DO J=I,I_BCHI
               SUM7=SUM7+VAR_DIV(J)*C(J)
            ENDDO
            SUM8=0.
            DO J=I+1,I_BCHI
               SUM8=SUM8+VAR_DIV(J)*C(J)**2
            ENDDO
            TERM8=VAR_RHO_W(I_BCHI)+SUM8
            RHO_W(I)=RHO_W(I_BCLO)-SUM4
            VAR_RHO_W(I)=(SUM7/SUM2)**2*TERM6+(SUM5/SUM2)**2*TERM8+
     $      (C(I)*SUM5/SUM2-B(I)*SUM7/SUM2)**2*VAR_DIV(I)
         ENDDO
      ENDIF

C  Calculate the vertical air velocity and its standard deviation at
C  levels.
      DO I=I_BCLO,I_BCHI
         W(I)=RHO_W(I)/RHO(I)
         SD_W(I)=SQRT(VAR_RHO_W(I))/RHO(I)
      ENDDO

C  Done.
      SUCCESS=.TRUE.
      RETURN
      END
