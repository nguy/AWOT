      SUBROUTINE W_UPWARD_DOWNWARD_LAYERS(Z1,TV1,Z2,TV2,Z3,RHO3,N,Z,DIV,
     $SD_DIV,I_BC,W_BC,SD_W_BC,BADDATA,W,SD_W)

C  Thomas Matejka NOAA/NSSL 23 September 1994

C  This subroutine calculates the vertical air velocity and its standard
C  deviation in a column from measurements of the mean horizontal
C  divergence in layers and a specified vertical air velocity at one
C  level in the column.

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

C  DIV is a one-dimensional real array with elements indexed from 1 to
C  N.  DIV(I) specifies the mean horizontal divergence (s**-1) in the
C  Ith layer from the bottom of the column.  If it is missing, it should
C  equal BADDATA.

C  SD_DIV is a one-dimensional real array with elements indexed from 1
C  to N.  SD_DIV(I) specifies the standard deviation of DIV(I) (s**-1).
C  If it is missing, it should equal BADDATA.  SD_DIV(I) must not be
C  missing for DIV(I) to be usable.

C  I_BC is an integer variable that specifies the level number above the
C  bottom of the column at which a vertical air velocity boundary
C  condition is specified.

C  W_BC is a real variable that specifies the vertical air velocity
C  boundary condition (m s**-1) at the I_BCth level above the bottom of
C  the column.  If it is missing, it should equal BADDATA.

C  SD_W_BC is a real variable that specifies the standard deviation of
C  W_BC (m s**-1).  If it is missing, it should equal BADDATA.  SD_W_BC
C  must not be missing for W_BC to be usable.

C  BADDATA is a real variable that indicates a missing value as
C  described.

C  Output:

C  W is a one-dimensional real array with elements indexed from 0 to N.
C  W(I) returns the vertical air velocity (m s**-1) at the Ith level
C  above the bottom of the column.  If it cannot be calculated, it is
C  returned as BADDATA.

C  SD_W is a one-dimensional real array with elements indexed from 0 to
C  N.  SD_W(I) returns the standard deviation of W(I) (m s**-1).  If it
C  cannot be calculated, it is returned as BADDATA.

      IMPLICIT NONE
      REAL DENCLR,AV_DENCLR
      INTEGER I,J,N,I_BC
      REAL Z1,Z2,Z3,TV1,TV2,RHO3,SUM1,SUM2,W_BC,SD_W_BC,BADDATA
      REAL DIV(N),SD_DIV(N),RHO_BAR(N),DEL_Z(N),VAR_DIV(N)
      REAL Z(0:N),W(0:N),SD_W(0:N),RHO(0:N),RHO_W(0:N),VAR_RHO_W(0:N)

C  Check whether the vertical air velocity boundary condition is
C  missing.
      IF(W_BC.EQ.BADDATA.OR.
     $SD_W_BC.EQ.BADDATA)THEN
         DO I=0,N
            W(I)=BADDATA
            SD_W(I)=BADDATA
         ENDDO
         RETURN
      ENDIF

C  Calculate the density at levels.
      DO I=0,N
         RHO(I)=DENCLR(Z1,TV1,Z2,TV2,Z3,RHO3,Z(I))
      ENDDO

C  Calculate the depth, the mean density, and the variance of the
C  divergence in layers.
      DO I=1,N
         DEL_Z(I)=Z(I)-Z(I-1)
         RHO_BAR(I)=AV_DENCLR(Z1,TV1,Z2,TV2,Z3,RHO3,Z(I-1),Z(I))
         IF(SD_DIV(I).NE.BADDATA)THEN
            VAR_DIV(I)=SD_DIV(I)**2
         ELSE
            VAR_DIV(I)=BADDATA
         ENDIF
      ENDDO

C  Calculate the vertical mass flux and its variance at the boundary
C  level.
      RHO_W(I_BC)=RHO(I_BC)*W_BC
      VAR_RHO_W(I_BC)=(RHO(I_BC)*SD_W_BC)**2

C  Calculate the vertical mass flux and its variance at higher levels.
      IF(I_BC.LT.N)THEN
         DO I=I_BC+1,N
            SUM1=0.
            SUM2=0.
            DO J=I_BC+1,I
               IF(DIV(J).NE.BADDATA.AND.
     $         VAR_DIV(J).NE.BADDATA)THEN
                  SUM1=SUM1+DIV(J)*RHO_BAR(J)*DEL_Z(J)
                  SUM2=SUM2+VAR_DIV(J)*(RHO_BAR(J)*DEL_Z(J))**2
               ELSE
                  RHO_W(I)=BADDATA
                  VAR_RHO_W(I)=BADDATA
                  GOTO 1
               ENDIF
            ENDDO
            RHO_W(I)=RHO_W(I_BC)-SUM1
            VAR_RHO_W(I)=VAR_RHO_W(I_BC)+SUM2
1           CONTINUE
         ENDDO
      ENDIF

C  Calculate the vertical mass flux and its variance at lower levels.
      IF(I_BC.GT.0)THEN
         DO I=0,I_BC-1
            SUM1=0.
            SUM2=0.
            DO J=I+1,I_BC
               IF(DIV(J).NE.BADDATA.AND.
     $         VAR_DIV(J).NE.BADDATA)THEN
                  SUM1=SUM1+DIV(J)*RHO_BAR(J)*DEL_Z(J)
                  SUM2=SUM2+VAR_DIV(J)*(RHO_BAR(J)*DEL_Z(J))**2
               ELSE
                  RHO_W(I)=BADDATA
                  VAR_RHO_W(I)=BADDATA
                  GOTO 2
               ENDIF
            ENDDO
            RHO_W(I)=RHO_W(I_BC)+SUM1
            VAR_RHO_W(I)=VAR_RHO_W(I_BC)+SUM2
2           CONTINUE
         ENDDO
      ENDIF

C  Calculate the vertical air velocity and its standard deviation at
C  levels.
      DO I=0,N
         IF(RHO_W(I).NE.BADDATA.AND.
     $   VAR_RHO_W(I).NE.BADDATA)THEN
            W(I)=RHO_W(I)/RHO(I)
            SD_W(I)=SQRT(VAR_RHO_W(I))/RHO(I)
         ELSE
            W(I)=BADDATA
            SD_W(I)=BADDATA
         ENDIF
      ENDDO

C  Done.
      RETURN
      END
