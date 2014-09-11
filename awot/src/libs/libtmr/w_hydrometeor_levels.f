      SUBROUTINE W_HYDROMETEOR_LEVELS(Z1,TV1,Z2,TV2,Z3,RHO3,N,Z,WP,
     $SD_WP,F,SD_F,BADDATA,DIV,SD_DIV,W,SD_W)

C  Thomas Matejka NOAA/NSSL 22 February 1995

C  This program calculates vertical air velocity, horizontal divergence,
C  and their standard deviations in a column from measurements of the
C  vertical hydrometeor velocity and hydrometeor terminal fall speed at
C  levels in the column.  The horizontal divergence is calculated at
C  levels.

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

C  W is a one-dimensional real array with elements indexed from 0 to N.
C  W(I) returns the vertical air velocity (m s**-1) at the Ith level
C  above the bottom of the column.  If it is not calculated, it is
C  returned as BADDATA.

C  SD_W is a one-dimensional real array with elements indexed from 0 to
C  N.  SD_W(I) returns the standard deviation of W(I).  If it is not
C  calculated, it is returned as BADDATA.

      IMPLICIT NONE
      INTEGER I,N
      REAL Z1,Z2,Z3,TV1,TV2,RHO3,BADDATA
      REAL Z(0:N),WP(0:N),SD_WP(0:N),F(0:N),SD_F(0:N),W(0:N),SD_W(0:N),
     $DIV(0:N),SD_DIV(0:N)

C  Calculate the the vertical air velocity and its variance and standard
C  deviation at levels.
      DO I=0,N
         IF(WP(I).NE.BADDATA.AND.
     $   SD_WP(I).NE.BADDATA.AND.
     $   F(I).NE.BADDATA.AND.
     $   SD_F(I).NE.BADDATA)THEN
            W(I)=WP(I)+F(I)
            SD_W(I)=SQRT(SD_WP(I)**2+SD_F(I)**2)
         ELSE
            W(I)=BADDATA
            SD_W(I)=BADDATA
         ENDIF
      ENDDO

C  Calculate the divergence and its standard deviation.
      CALL DIV_FROM_W_LEVELS(Z1,TV1,Z2,TV2,Z3,RHO3,N,Z,W,SD_W,BADDATA,
     $DIV,SD_DIV)

C  Done.
      RETURN
      END
