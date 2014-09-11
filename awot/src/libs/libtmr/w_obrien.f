      SUBROUTINE W_OBRIEN(Z1,TV1,Z2,TV2,Z3,RHO3,
     $N,Z,DIV,SD_DIV,
     $I_BCLO,W_BCLO,I_BCHI,W_BCHI,
     $BADDATA,
     $DIV_OUT,W,
     $SUCCESS)

C  Thomas Matejka NOAA/NSSL 29 November 2001

C  This subroutine calculates adjusted horizontal divergence and
C  vertical air velocity in a column from measurements of horizontal
C  divergence and two specified vertical air velocity boundary
C  conditions with a version of the O'Brien variational formulation.
C  The solution is found at and between the boundaries only if the
C  divergence is not missing at or between the boundaries.

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

C  SD_DIV is a one-dimensional real array with elements indexed from 1
C  to N.  SD_DIV(I) specifies the rms error of DIV(I) (s**-1).  If it is
C  missing, it should equal BADDATA.

C  I_BCLO is an integer variable that specifies that the lower boundary
C  is at Z(I_BCLO).

C  W_BCLO is a real variable that specifies the vertical air velocity
C  boundary condition (m s**-1) at Z(I_BCLO).  If it is missing, it
C  should equal BADDATA.

C  I_BCHI is an integer variable that specifies that the upper boundary
C  is at Z(I_BCHI).

C  W_BCHI is a real variable that specifies the vertical air velocity
C  boundary condition (m s**-1) at Z(I_BCHI).  If it is missing, it
C  should equal BADDATA.

C  BADDATA is a real variable that indicates a missing value as
C  described.

C  Output:

C  DIV_OUT is a one-dimensional real array with elements indexed from 1
C  to N.  DIV_OUT(I) returns the adjusted horizontal divergence (s**-1)
C  at Z(I).  If it is missing, it returns BADDATA.  DIV_OUT may be the
C  same as DIV, in which case DIV is overwritten.

C  W is a one-dimensional real array with elements indexed from 1 to N.
C  W(I) returns the vertical air velocity (m s**-1) at the Z(I).  If it
C  is missing, it returns BADDATA.

C  SUCCESS is a logical variable that returns .TRUE. if and only if a
C  vertical air velocity could be calculated at at least one level.

      IMPLICIT NONE
      INCLUDE 'tmrlib.inc'
      REAL,EXTERNAL::DENCLR,AV_DENCLR
      LOGICAL::SUCCESS
      INTEGER::I,N,I_BCLO,I_BCHI
      REAL::Z1,Z2,Z3,TV1,TV2,RHO3,SUM1,SUM2,BADDATA,W_BCLO,W_BCHI
      REAL,DIMENSION(N)::RHO_BAR_BELOW,RHO_BAR_ABOVE,DEL_Z_BELOW,
     $DEL_Z_ABOVE,Z,RHO,DIV,SD_DIV,DIV_OUT,W,DIV_TEMP

C  Copy the input array in case it is overwritten.
      DO I=1,N
         DIV_TEMP(I)=DIV(I)
      ENDDO

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

C  Initialize the adjusted divergence and the vertical air velocity to
C  missing.
      DO I=1,N
         DIV_OUT(I)=BADDATA
         W(I)=BADDATA
      ENDDO

C  Check whether a vertical air velocity boundary condition is missing.
      IF(W_BCLO.EQ.BADDATA.OR.
     $W_BCHI.EQ.BADDATA.OR.
     $I_BCLO.GE.I_BCHI)THEN
         SUCCESS=.FALSE.
         RETURN
      ENDIF

C  Check whether a divergence is missing between the boundaries.
      DO I=I_BCLO,I_BCHI
         IF(DIV_TEMP(I).EQ.BADDATA.OR.
     $   SD_DIV(I).EQ.BADDATA)THEN
            SUCCESS=.FALSE.
            RETURN
         ENDIF
      ENDDO

C  Calculate sums.
      SUM1=0.
      SUM2=0.
      DO I=I_BCLO+1,I_BCHI
         SUM1=SUM1+
     $   DIV_TEMP(I-1)*RHO_BAR_ABOVE(I-1)*DEL_Z_ABOVE(I-1)+
     $   DIV_TEMP(I)*RHO_BAR_BELOW(I)*DEL_Z_BELOW(I)
         SUM2=SUM2+
     $   SD_DIV(I-1)**2*RHO_BAR_ABOVE(I-1)*DEL_Z_ABOVE(I-1)+
     $   SD_DIV(I)**2*RHO_BAR_BELOW(I)*DEL_Z_BELOW(I)
      ENDDO

C  Calculate the adjusted divergence.
      DO I=I_BCLO,I_BCHI
         DIV_OUT(I)=DIV_TEMP(I)-
     $   SD_DIV(I)**2*
     $   (RHO(I_BCHI)*W_BCHI-RHO(I_BCLO)*W_BCLO+SUM1)/
     $   SUM2
      ENDDO

C  Calculate the vertical air velocity.
      CALL W_DOWN(Z1,TV1,Z2,TV2,Z3,RHO,
     $N,Z,DIV_OUT,
     $I_BCHI,W_BCHI,
     $BADDATA,
     $W,
     $SUCCESS)
      IF(.NOT.SUCCESS)THEN
         WRITE(TMRLIB_MESSAGE_UNIT,*)'W_OBRIEN:  UNEXPECTED FALSE ',
     $   'SUCCESS from W_DOWN.'
         STOP
      ENDIF

C  Use the specified vertical air velocity boundary condition at the
C  lower boundary instead of the calculated value to eliminate
C  computational errors.
      W(I_BCLO)=W_BCLO

C  Done.
      SUCCESS=.TRUE.

      END SUBROUTINE W_OBRIEN
