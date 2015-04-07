      SUBROUTINE DIV_FROM_W_LAYERS(Z1,TV1,Z2,TV2,Z3,RHO3,N,Z,W,SD_W,
     $BADDATA,DIV,SD_DIV)

C  Thomas Matejka NOAA/NSSL 22 February 1995

C  This program calculates horizontal divergence and its standard
C  deviation in a column from measurements of the vertical velocity at
C  levels in the column.  The horizontal divergence is calculated as the
C  means in layers.

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

C  W is a one-dimensional real array with elements indexed from 0 to N.
C  W(I) specifies the vertical air velocity (m s**-1) at the Ith level
C  above the bottom of the column.  If it is not calculated, it is
C  returned as BADDATA.

C  SD_W is a one-dimensional real array with elements indexed from 0 to
C  N.  SD_W(I) specifies the standard deviation of W(I).  If it is not
C  calculated, it is returned as BADDATA.

C  BADDATA is a real variable that indicates a missing value as
C  described.

C  Output:

C  DIV is a one-dimensional real array with elements indexed from 1 to
C  N.  DIV_ADJ(I) returns the adjusted horizontal divergence(s**-1) in
C  the Ith layer from the bottom of the column.  If it is not
C  calculated, it is returned as BADDATA.

C  SD_DIV is a one-dimensional real array with elements indexed from 1
C  to N.  SD_DIV(I) returns the standard deviation of DIV(I).  If it is
C  not calculated, it is returned as BADDATA.

      IMPLICIT NONE
      REAL DENCLR,AV_DENCLR
      INTEGER I,N
      REAL Z1,Z2,Z3,TV1,TV2,RHO3,BADDATA
      REAL DIV(N),SD_DIV(N),RHO_BAR(N),DEL_Z(N)
      REAL Z(0:N),W(0:N),SD_W(0:N),RHO(0:N),RHO_W(0:N),VAR_RHO_W(0:N)

C  Calculate the density, the vertical air velocity, its standard
C  deviation, the vertical mass flux, and its variance at levels.
      DO I=0,N
         RHO(I)=DENCLR(Z1,TV1,Z2,TV2,Z3,RHO3,Z(I))
         IF(W(I).NE.BADDATA.AND.
     $   SD_W(I).NE.BADDATA)THEN
            RHO_W(I)=RHO(I)*W(I)
            VAR_RHO_W(I)=(RHO(I)*SD_W(I))**2
         ELSE
            RHO_W(I)=BADDATA
            VAR_RHO_W(I)=BADDATA
         ENDIF
      ENDDO

C  Calculate the depth and the mean density in layers.
      DO I=1,N
         DEL_Z(I)=Z(I)-Z(I-1)
         RHO_BAR(I)=AV_DENCLR(Z1,TV1,Z2,TV2,Z3,RHO3,Z(I-1),Z(I))
      ENDDO

C  Calculate the divergence and its standard deviation in layers.
      DO I=1,N
         IF(RHO_W(I-1).NE.BADDATA.AND.
     $   VAR_RHO_W(I-1).NE.BADDATA.AND.
     $   RHO_W(I).NE.BADDATA.AND.
     $   VAR_RHO_W(I).NE.BADDATA)THEN
            DIV(I)=(RHO_W(I-1)-RHO_W(I))/(RHO_BAR(I)*DEL_Z(I))
            SD_DIV(I)=SQRT(VAR_RHO_W(I-1)+VAR_RHO_W(I))/(RHO_BAR(I)*
     $      DEL_Z(I))
         ELSE
            DIV(I)=BADDATA
            SD_DIV(I)=BADDATA
         ENDIF
      ENDDO

C  Done.
      RETURN
      END
