      SUBROUTINE W_PARTIAL_OBRIEN_OR_DOWN(Z1,TV1,Z2,TV2,Z3,RHO3,
     $N,Z,DIV,SD_DIV,
     $I_BCLO,W_BCLO,I_BCHI,W_BCHI,DELZ_MAX,
     $BADDATA,
     $DIV_OUT,W,
     $DID_PARTIAL_OBRIEN,DID_DOWN,SUCCESS)

C  Thomas Matejka NOAA/NSSL 29 November 2001

C  This subroutine calculates adjusted horizontal divergence and
C  vertical air velocity in a column with W_PARTIAL_OBRIEN.  If that is
C  unsuccessful, the subroutine calculates vertical air velocity in the
C  column with W_DOWN.

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
C  to N.  DIV_OUT(I) returns the horizontal divergence (s**-1) at Z(I),
C  possibly filled and adjusted.  If it is missing, it returns BADDATA.
C  DIV_OUT may be the same as DIV, in which case DIV is overwritten.

C  W is a one-dimensional real array with elements indexed from 1 to N.
C  W(I) returns the vertical air velocity (m s**-1) at the Z(I).  If it
C  is missing, it returns BADDATA.

C  DID_PARTIAL_OBRIEN is a logical variable that returns .TRUE. if and
C  only if W_PARTIAL_OBRIEN was used.

C  DID_DOWN is a logical variable that returns .TRUE. if and only if
C  W_DOWN was used.

C  SUCCESS is a logical variable that returns .TRUE. if and only if a
C  vertical air velocity could be calculated at at least one level.

      IMPLICIT NONE
      LOGICAL::DID_PARTIAL_OBRIEN,DID_DOWN,SUCCESS
      INTEGER::I,N,I_BCLO,I_BCHI
      REAL::Z1,Z2,Z3,TV1,TV2,RHO3,DELZ_MAX,BADDATA,W_BCLO,W_BCHI
      REAL,DIMENSION(N)::Z,DIV,SD_DIV,DIV_OUT,W,DIV_TEMP

C  Copy the input array in case it is overwritten.
      DO I=1,N
         DIV_TEMP(I)=DIV(I)
      ENDDO

C  Initialize method used flags.
      DID_DOWN=.FALSE.
      DID_PARTIAL_OBRIEN=.FALSE.

C  Partial O'Brien solution.
      CALL W_PARTIAL_OBRIEN(Z1,TV1,Z2,TV2,Z3,RHO3,
     $N,Z,DIV_TEMP,SD_DIV,
     $I_BCLO,W_BCLO,I_BCHI,W_BCHI,DELZ_MAX,
     $BADDATA,
     $DIV_OUT,W,
     $SUCCESS)
      IF(SUCCESS)THEN
         DID_PARTIAL_OBRIEN=.TRUE.

C  If the partial O'Brien solution was not successful, perform a
C  downward integration.  In this case, the divergence is unchanged.
      ELSE
         CALL W_DOWN(Z1,TV1,Z2,TV2,Z3,RHO3,
     $   N,Z,DIV_TEMP,
     $   I_BCHI,W_BCHI,
     $   BADDATA,
     $   W,
     $   SUCCESS)
         IF(SUCCESS)THEN
            DID_DOWN=.TRUE.
         ENDIF
         DO I=1,N
            DIV_OUT(I)=DIV_TEMP(I)
         ENDDO
      ENDIF

      END SUBROUTINE W_PARTIAL_OBRIEN_OR_DOWN
