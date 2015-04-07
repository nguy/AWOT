      SUBROUTINE W_UP(Z1,TV1,Z2,TV2,Z3,RHO3,
     $N,Z,DIV,
     $I_BC,W_BC,
     $BADDATA,
     $W,
     $SUCCESS)

C  Thomas Matejka NOAA/NSSL 29 November 2001

C  This subroutine calculates vertical air velocity in a column from
C  measurements of horizontal divergence and one specified vertical air
C  velocity boundary condition by upward integration of the mass
C  continuity equation.

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

C  N is an integer variable that specifies the number of data levels in
C  the column.

C  Z is a one-dimensional real array with elements indexed from 1 to N.
C  Z(I) specifies the altitude (m) of the Ith level, counting from the
C  bottom of the column.

C  DIV is a one-dimensional real array with elements indexed from 1 to
C  N.  DIV(I) specifies the horizontal divergence (s**-1) at Z(I).  If
C  it is missing, it should equal BADDATA.

C  I_BC is an integer variable that specifies that the boundary is at
C  Z(I_BC).

C  W_BC is a real variable that specifies the vertical air velocity
C  boundary condition (m s**-1) at Z(I_BC).  If it is missing, it should
C  equal BADDATA.

C  BADDATA is a real variable that indicates a missing value as
C  described.

C  Output:

C  W is a one-dimensional real array with elements indexed from 1 to N.
C  W(I) returns the vertical air velocity (m s**-1) at Z(I).  If it is
C  missing, it returns BADDATA.

C  SUCCESS is a logical variable that returns .TRUE. if and only if a
C  vertical air velocity could be calculated at at least one level.

      IMPLICIT NONE
      REAL,EXTERNAL::DENCLR,AV_DENCLR
      LOGICAL::SUCCESS
      INTEGER::I,N,I_BC
      REAL::Z1,Z2,Z3,TV1,TV2,RHO3,W_BC,BADDATA
      REAL,DIMENSION(N)::RHO_BAR_BELOW,RHO_BAR_ABOVE,DEL_Z_BELOW,
     $DEL_Z_ABOVE,DIV,Z,W,RHO,RHO_W

C  Calculate the density at levels.
      DO I=1,N
         RHO(I)=DENCLR(Z1,TV1,Z2,TV2,Z3,RHO3,Z(I))
      ENDDO

C  Calculate the depth of and the mean density in half-layers below
C  levels.
      IF(N.GE.2)THEN
         DO I=2,N
            DEL_Z_BELOW(I)=(Z(I)-Z(I-1))/2.
            RHO_BAR_BELOW(I)=AV_DENCLR(Z1,TV1,Z2,TV2,Z3,RHO3,
     $      (Z(I-1)+Z(I))/2.,Z(I))
         ENDDO
      ENDIF

C  Calculate the depth of and the mean density in half-layers above
C  levels.
      IF(N.GE.2)THEN
         DO I=1,N-1
            DEL_Z_ABOVE(I)=(Z(I+1)-Z(I))/2.
            RHO_BAR_ABOVE(I)=AV_DENCLR(Z1,TV1,Z2,TV2,Z3,RHO3,Z(I),
     $      (Z(I)+Z(I+1))/2.)
         ENDDO
      ENDIF

C  Initialize the vertical mass flux to missing.
      DO I=1,N
         RHO_W(I)=BADDATA
      ENDDO

C  Calculate the vertical mass flux at the boundary.
      IF(W_BC.NE.BADDATA)THEN
         RHO_W(I_BC)=RHO(I_BC)*W_BC
      ENDIF

C  Calculate the vertical mass flux above the boundary.
      SUCCESS=.FALSE.
      IF(I_BC.LT.N)THEN
         DO I=I_BC+1,N
            IF(RHO_W(I-1).EQ.BADDATA.OR.
     $      DIV(I-1).EQ.BADDATA.OR.
     $      DIV(I).EQ.BADDATA)THEN
               EXIT
            ENDIF
            RHO_W(I)=RHO_W(I-1)-
     $      DIV(I-1)*RHO_BAR_ABOVE(I-1)*DEL_Z_ABOVE(I-1)-
     $      DIV(I)*RHO_BAR_BELOW(I)*DEL_Z_BELOW(I)
            SUCCESS=.TRUE.
         ENDDO
      ENDIF

C  Calculate the vertical air velocity.
      DO I=1,N
         IF(RHO_W(I).NE.BADDATA)THEN
            W(I)=RHO_W(I)/RHO(I)
         ELSE
            W(I)=BADDATA
         ENDIF
      ENDDO

      END SUBROUTINE W_UP
