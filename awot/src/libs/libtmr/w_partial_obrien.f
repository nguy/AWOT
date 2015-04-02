      SUBROUTINE W_PARTIAL_OBRIEN(Z1,TV1,Z2,TV2,Z3,RHO3,
     $N,Z,DIV,SD_DIV,
     $I_BCLO,W_BCLO,I_BCHI,W_BCHI,DELZ_MAX,
     $BADDATA,
     $DIV_OUT,W,
     $SUCCESS)

C  Thomas Matejka NOAA/NSSL 29 November 2001

C  This subroutine calculates adjusted horizontal divergence and
C  vertical air velocity in a column from measurements of horizontal
C  divergence and two specified vertical air velocity boundary
C  conditions with a version of the O'Brien variational formulation.
C  Divergences that are missing at and immediately above the lower
C  boundary through a depth that is less than a specified threshold are
C  filled.  Otherwise, the solution is found at and between the
C  boundaries only if the divergence is not missing at or between the
C  boundaries.  The missing divergences toward the bottom of the column
C  are estimated as a linear combination of the divergence at the lowest
C  data level and the divergence necessary to mass balance the column.
C  The former is weighted more heavily when the depth of the layer of
C  missing data is smaller.  The solution is in effect a partial
C  adjustment of a downward integration solution toward the O'Brien
C  solution.  When the depth of the layer of missing data is zero, the
C  solution is identical to the O'Brien solution, since the bottom
C  boundary is contiguous with the divergence data as required.  As the
C  depth of the layer of missing data becomes larger (but still less
C  than the threshold), the solution increasingly resembles a downward
C  integration, because the bottom boundary becomes more detatched from
C  the divergence data and, hence, less usable.

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

C  DELZ_MAX is a real variable that specifies the depth (m) of the layer
C  of missing divergences at and above the lower boundary that must not
C  be equaled or exceeded.

C  BADDATA is a real variable that indicates a missing value as
C  described.

C  Output:

C  DIV_OUT is a one-dimensional real array with elements indexed from 1
C  to N.  DIV_OUT(I) returns the possibly filled and adjusted horizontal
C  divergence (s**-1) at Z(I).  If it is missing, it returns BADDATA.
C  DIV_OUT may be the same as DIV, in which case DIV is overwritten.

C  W is a one-dimensional real array with elements indexed from 1 to N.
C  W(I) returns the vertical air velocity (m s**-1) at the Z(I).  If it
C  is missing, it returns BADDATA.

C  SUCCESS is a logical variable that returns .TRUE. if and only if a
C  vertical air velocity could be calculated at at least one level.

      IMPLICIT NONE
      INCLUDE 'tmrlib.inc'
      REAL,EXTERNAL::DENCLR,AV_DENCLR
      LOGICAL::FOUND,SUCCESS
      INTEGER::I,N,I_BCLO,I_BCHI,I_DIVLO
      REAL::Z1,Z2,Z3,TV1,TV2,RHO3,DELZ_MAX,SUM,BADDATA,W_BCLO,W_BCHI,
     $DIV_MB,ALPHA
      REAL,DIMENSION(N)::RHO_BAR_BELOW,RHO_BAR_ABOVE,DEL_Z_BELOW,
     $DEL_Z_ABOVE,Z,RHO,DIV,SD_DIV,DIV_OUT,W,DIV_STORE,SD_DIV_STORE,
     $W_DUM,DIV_TEMP

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

C  Find the lowest level of divergence data.
      FOUND=.FALSE.
      DO I=I_BCLO,I_BCHI
         IF(DIV_TEMP(I).NE.BADDATA.AND.
     $   SD_DIV(I).NE.BADDATA)THEN
            FOUND=.TRUE.
            I_DIVLO=I
            EXIT
         ENDIF
      ENDDO
      IF(.NOT.FOUND)THEN
         SUCCESS=.FALSE.
         RETURN
      ENDIF

C  Check whether the distance between the lower boundary and the lowest
C  level of divergence data equals or exceeds the threshold.
      IF(Z(I_DIVLO)-Z(I_BCLO).GE.DELZ_MAX)THEN
         SUCCESS=.FALSE.
         RETURN
      ENDIF

C  Check whether a divergence is missing between the lowest level of
C  divergence data and the upper boundary.
      DO I=I_DIVLO,I_BCHI
         IF(DIV_TEMP(I).EQ.BADDATA.OR.
     $   SD_DIV(I).EQ.BADDATA)THEN
            SUCCESS=.FALSE.
            RETURN
         ENDIF
      ENDDO

C  Calculate the divergence that would mass balance the column if used
C  for the missing divergences.
      IF(I_DIVLO.GT.I_BCLO)THEN
         CALL W_DOWN(Z1,TV1,Z2,TV2,Z3,RHO3,
     $   N,Z,DIV_TEMP,
     $   I_BCHI,W_BCHI,
     $   BADDATA,
     $   W_DUM,
     $   SUCCESS)
         IF(.NOT.SUCCESS)THEN
            WRITE(TMRLIB_MESSAGE_UNIT,*)'W_PARTIAL_OBRIEN:  ',
     $      'UNEXPECTED FALSE SUCCESS FROM W_DOWN.'
            STOP
         ENDIF
         SUM=RHO_BAR_ABOVE(I_BCLO)*DEL_Z_ABOVE(I_BCLO)
         IF(I_DIVLO-1.GE.I_BCLO+1)THEN
            DO I=I_BCLO+1,I_DIVLO-1
               SUM=SUM+
     $         RHO_BAR_BELOW(I)*DEL_Z_BELOW(I)+
     $         RHO_BAR_ABOVE(I)*DEL_Z_ABOVE(I)
            ENDDO
         ENDIF
         DIV_MB=(RHO(I_BCLO)*W_BCLO-
     $   RHO(I_DIVLO)*W_DUM(I_DIVLO)-
     $   RHO_BAR_BELOW(I_DIVLO)*DEL_Z_BELOW(I_DIVLO)*DIV_TEMP(I_DIVLO))/
     $   SUM
      ENDIF

C  The missing divergences can be estimated as a linear combination of
C  the lowest divergence datum and the divergence needed to mass balance
C  the column.  The proportionality is determined by the depth of the
C  layer of missing data.
      IF(I_DIVLO.GT.I_BCLO)THEN
         ALPHA=(DELZ_MAX-(Z(I_DIVLO)-Z(I_BCLO)))/DELZ_MAX
         DO I=I_BCLO,I_DIVLO-1
            DIV_STORE(I)=ALPHA*DIV_TEMP(I_DIVLO)+
     $      (1.-ALPHA)*DIV_MB
            SD_DIV_STORE(I)=SD_DIV(I_DIVLO)
         ENDDO
      ENDIF
      DO I=I_DIVLO,I_BCHI
         DIV_STORE(I)=DIV_TEMP(I)
         SD_DIV_STORE(I)=SD_DIV(I)
      ENDDO

C  Calculate a variational solution from the filled column.
      CALL W_OBRIEN(Z1,TV1,Z2,TV2,Z3,RHO3,
     $N,Z,DIV_STORE,SD_DIV_STORE,
     $I_BCLO,W_BCLO,I_BCHI,W_BCHI,
     $BADDATA,
     $DIV_OUT,W,
     $SUCCESS)
      IF(.NOT.SUCCESS)THEN
         WRITE(TMRLIB_MESSAGE_UNIT,*)'W_PARTIAL_OBRIEN:  UNEXPECTED ',
     $   'FALSE SUCCESS FROM W_OBRIEN.'
         STOP
      ENDIF

C  Done.
      SUCCESS=.TRUE.

      END SUBROUTINE W_PARTIAL_OBRIEN
