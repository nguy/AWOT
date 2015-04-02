      SUBROUTINE COMPLETE_VAR_COL_SOLN_DRIVER(Z1,TV1,Z2,TV2,Z3,RHO3,N,Z,
     $DIV,SD_DIV,WP,SD_WP,F,SD_F,W_CONSTRAINT,SD_W_CONSTRAINT,BADDATA,
     $DIV_ADJ,SD_DIV_ADJ,WP_ADJ,SD_WP_ADJ,F_ADJ,SD_F_ADJ,W,SD_W,SUCCESS)

C  Thomas Matejka NOAA/NSSL 23 September 1994

C  This program calculates variationally adjusted horizontal divergence,
C  variationally adjusted vertical hydrometeor velocity, variationally
C  adjusted hydrometeor terminal fall speed, vertical air velocity, and
C  their standard deviations in a column from measurements of the
C  horizontal divergence, vertical hydrometeor velocity, and hydrometeor
C  terminal fall speed at levels in the column.  The solution is found
C  only at and beween the highest set of consecutive levels having
C  complete data, but only if the highest and lowest levels in this set
C  are not identical.

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

C  SD_W_CONSTRAINT is a one-dimensional real array with elements indexed
C  from 0 to N.  SD_W_CONSTRAINT(I) specifies the standard deviation of
C  W_CONSTRAINT(I) (m s**-1).  If it is missing, it should equal
C  BADDATA.  SD_W_CONSTRAINT(I) must not be missing for W_CONSTRAINT(I)
C  to be usable.

C  BADDATA is a real variable that indicates a missing value as
C  described.

C  Output:

C  DIV_ADJ is a one-dimensional real array with elements indexed from 0
C  to N.  DIV_ADJ(I) returns the adjusted horizontal divergence(s**-1)
C  at the Ith level above the bottom of the column.  If it is not
C  calculated, it is returned as BADDATA.

C  SD_DIV_ADJ is a one-dimensional real array with elements indexed from
C  0 to N.  SD_DIV_ADJ(I) returns the standard deviation of the adjusted
C  horizontal divergence(s**-1) at the Ith level above the bottom of the
C  column.  If it is not calculated, it is returned as BADDATA.

C  WP_ADJ is a one-dimensional real array with elements indexed from 0
C  to N.  WP_ADJ(I) returns the adjusted mean echo-power-weighted
C  vertical hydrometeor velocity (m s**-1) at the Ith level above the
C  bottom of the column.  If it is not calculated, it is returned as
C  BADDATA.

C  SD_WP_ADJ is a one-dimensional real array with elements indexed from
C  0 to N.  SD_WP_ADJ(I) returns the standard deviation of the adjusted
C  mean echo-power-weighted vertical hydrometeor velocity (m s**-1) at
C  the Ith level above the bottom of the column. If it is not
C  calculated, it is returned as BADDATA.

C  F_ADJ is a one-dimensional real array with elements indexed from 0 to
C  N.  F_ADJ(I) returns the adjusted mean echo-power-weighted terminal
C  fall speed (m s**-1) at the Ith level above the bottom of the column.
C  If it is not calculated, it is returned as BADDATA.

C  SD_F_ADJ is a one-dimensional real array with elements indexed from 0
C  to N.  F_ADJ(I) returns the standard deviation of adjusted mean
C  echo-power-weighted terminal fall speed (m s**-1) at the Ith level
C  above the bottom of the column.  If it is not calculated, it is
C  returned as BADDATA.

C  W is a one-dimensional real array with elements indexed from 0 to N.
C  W(I) returns the vertical air velocity (m s**-1) at the Ith level
C  above the bottom of the column.  If it is not calculated, it is
C  returned as BADDATA.

C  SD_W is a one-dimensional real array with elements indexed from 0 to
C  N.  W(I) returns the standard deviation of vertical air velocity (m
C  s**-1) at the Ith level above the bottom of the column.  If it is not
C  calculated, it is returned as BADDATA.

C  SUCCESS is a logical variable that is returned .TRUE. if and only if
C  a solution of any sort was able to be found.

      IMPLICIT NONE
      LOGICAL SUCCESS
      INTEGER N,K
      REAL Z1,TV1,Z2,TV2,Z3,RHO3,BADDATA
      REAL Z(0:N),DIV(0:N),WP(0:N),F(0:N),SD_DIV(0:N),SD_WP(0:N),
     $SD_F(0:N),W_CONSTRAINT(0:N),SD_W_CONSTRAINT(0:N),DIV_ADJ(0:N),
     $SD_DIV_ADJ(0:N),WP_ADJ(0:N),SD_WP_ADJ(0:N),F_ADJ(0:N),
     $SD_F_ADJ(0:N),W(0:N),SD_W(0:N)

C  Calculate the profiles of adjusted horizontal divergence, adjusted
C  vertical hydrometeor velocity, adjusted hydrometeor terminal fall
C  speed, and vertical air velocity.
      CALL COMPLETE_VAR_COL_SOLN(Z1,TV1,Z2,TV2,Z3,RHO3,N,Z,DIV,SD_DIV,
     $WP,SD_WP,F,SD_F,W_CONSTRAINT,BADDATA,DIV_ADJ,WP_ADJ,F_ADJ,W,
     $SUCCESS)

C  Calculate the standard deviations of the profiles of adjusted
C  horizontal divergence, adjusted vertical hydrometeor velocity,
C  adjusted hydrometeor terminal fall speed, and vertical air velocity.
      IF(SUCCESS)THEN
         CALL COMPLETE_VAR_COL_SOLN_SD(Z1,TV1,Z2,TV2,Z3,RHO3,N,Z,DIV,
     $   SD_DIV,WP,SD_WP,F,SD_F,W_CONSTRAINT,SD_W_CONSTRAINT,BADDATA,
     $   SD_DIV_ADJ,SD_WP_ADJ,SD_F_ADJ,SD_W)
      ELSE
         DO K=0,N
            SD_DIV_ADJ(K)=BADDATA
            SD_WP_ADJ(K)=BADDATA
            SD_F_ADJ(K)=BADDATA
            SD_W(K)=BADDATA
         ENDDO
      ENDIF

C  Done.
      RETURN
      END
