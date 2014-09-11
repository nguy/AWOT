      SUBROUTINE DIV_FROM_W_LEVELS(Z1,TV1,Z2,TV2,Z3,RHO3,N,Z,W,SD_W,
     $BADDATA,DIV,SD_DIV)

C  Thomas Matejka NOAA/NSSL 28 February 1996

C  This program calculates horizontal divergence and its standard
C  deviation at levels in a column from measurements of the vertical air
C  velocity at levels in the column.

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

C  DIV is a one-dimensional real array with elements indexed from 0 to
C  N.  DIV_ADJ(I) returns the adjusted horizontal divergence(s**-1) at
C  the Ith level above the bottom of the column.  If it is not
C  calculated, it is returned as BADDATA.

C  SD_DIV is a one-dimensional real array with elements indexed from 0
C  to N.  SD_DIV(I) returns the standard deviation of DIV(I).  If it is
C  not calculated, it is returned as BADDATA.

      IMPLICIT NONE
      REAL DENCLR,AV_DENCLR
      INTEGER I,N
      REAL Z1,Z2,Z3,TV1,TV2,RHO3,BADDATA
      REAL RHO_BAR_BELOW(N),DEL_Z_BELOW(N),RHO_MIDPOINT(N)
      REAL Z(0:N),W(0:N),SD_W(0:N),DIV(0:N),SD_DIV(0:N),RHO(0:N),
     $VAR_W(0:N),A(0:N)
      REAL RHO_BAR_ABOVE(0:N-1),DEL_Z_ABOVE(0:N-1)

C  Calculate the density and the variance of the vertical air velocity
C  at levels.
      DO I=0,N
         RHO(I)=DENCLR(Z1,TV1,Z2,TV2,Z3,RHO3,Z(I))
         IF(SD_W(I).NE.BADDATA)THEN
            VAR_W(I)=SD_W(I)**2
         ELSE
            VAR_W(I)=BADDATA
         ENDIF
      ENDDO

C  Calculate the depth and the mean density in half-layers below levels
C  and the density at the midpoint of layers.
      DO I=1,N
         DEL_Z_BELOW(I)=(Z(I)-Z(I-1))/2.
         RHO_BAR_BELOW(I)=AV_DENCLR(Z1,TV1,Z2,TV2,Z3,RHO3,
     $   (Z(I-1)+Z(I))/2.,Z(I))
         RHO_MIDPOINT(I)=DENCLR(Z1,TV1,Z2,TV2,Z3,RHO3,(Z(I-1)+Z(I))/2.)
      ENDDO

C  Calculate the depth and the mean density in half-layers above levels.
      DO I=0,N-1
         DEL_Z_ABOVE(I)=(Z(I+1)-Z(I))/2.
         RHO_BAR_ABOVE(I)=AV_DENCLR(Z1,TV1,Z2,TV2,Z3,RHO3,Z(I),
     $   (Z(I)+Z(I+1))/2.)
      ENDDO

C  Calculate the mass per horizontal area appropriate for each level.
      A(0)=RHO_BAR_ABOVE(0)*DEL_Z_ABOVE(0)
      IF(N.GE.2)THEN
         DO I=1,N-1
            A(I)=RHO_BAR_ABOVE(I)*DEL_Z_ABOVE(I)+RHO_BAR_BELOW(I)*
     $      DEL_Z_BELOW(I)
         ENDDO
      ENDIF
      A(N)=RHO_BAR_BELOW(N)*DEL_Z_BELOW(N)

C  Calculate the divergence and its standard deviation at levels.
      IF(W(0).NE.BADDATA.AND.
     $VAR_W(0).NE.BADDATA.AND.
     $W(1).NE.BADDATA.AND.
     $VAR_W(1).NE.BADDATA)THEN
         DIV(0)=(W(0)*(RHO(0)-RHO_MIDPOINT(1)/2.)-W(1)*RHO_MIDPOINT(1)/
     $   2.)/A(0)
         SD_DIV(0)=SQRT(VAR_W(0)*(RHO(0)-RHO_MIDPOINT(1)/2.)**2+
     $   VAR_W(1)*(RHO_MIDPOINT(1)/2.)**2)/A(0)
      ELSE
         DIV(0)=BADDATA
         SD_DIV(0)=BADDATA
      ENDIF
      IF(N.GE.2)THEN
         DO I=1,N-1
            IF(W(I-1).NE.BADDATA.AND.
     $      VAR_W(I-1).NE.BADDATA.AND.
     $      W(I).NE.BADDATA.AND.
     $      VAR_W(I).NE.BADDATA.AND.
     $      W(I+1).NE.BADDATA.AND.
     $      VAR_W(I+1).NE.BADDATA)THEN
               DIV(I)=(W(I-1)*(RHO_MIDPOINT(I)/2.)+W(I)*
     $         (RHO_MIDPOINT(I)/2.-RHO_MIDPOINT(I+1)/2.)-W(I+1)*
     $         (RHO_MIDPOINT(I+1)/2.))/A(I)
               SD_DIV(I)=SQRT(VAR_W(I-1)*(RHO_MIDPOINT(I)/2.)**2+
     $         VAR_W(I)*(RHO_MIDPOINT(I)/2.-RHO_MIDPOINT(I+1)/2.)**2+
     $         VAR_W(I+1)*(RHO_MIDPOINT(I+1)/2.)**2)/A(I)
            ELSE
               DIV(I)=BADDATA
               SD_DIV(I)=BADDATA
            ENDIF
         ENDDO
      ENDIF
      IF(W(N-1).NE.BADDATA.AND.
     $VAR_W(N-1).NE.BADDATA.AND.
     $W(N).NE.BADDATA.AND.
     $VAR_W(N-1).NE.BADDATA)THEN
         DIV(N)=(W(N-1)*RHO_MIDPOINT(N)/2.+W(N)*(RHO_MIDPOINT(N)/2.-
     $   RHO(N)))/A(N)
         SD_DIV(N)=SQRT(VAR_W(N-1)*(RHO_MIDPOINT(N)/2.)**2+VAR_W(N)*
     $   (RHO_MIDPOINT(N)/2.-RHO(N))**2)/A(N)
      ELSE
         DIV(N)=BADDATA
         SD_DIV(N)=BADDATA
      ENDIF

C  Done.
      RETURN
      END
