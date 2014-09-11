      SUBROUTINE COMPLETE_VAR_COL_SOLN_SD(Z1,TV1,Z2,TV2,Z3,RHO3,N,Z,DIV,
     $SD_DIV,WP,SD_WP,F,SD_F,W_CONSTRAINT,SD_W_CONSTRAINT,BADDATA,
     $SD_DIV_ADJ,SD_WP_ADJ,SD_F_ADJ,SD_W)

C  Diana Bartels and Thomas Matejka NOAA/NSSL 23 September 1994

C  This subroutine calculates the standard deviations of adjusted
C  horizontal divergence, adjusted vertical hydrometeor velocity,
C  adjusted terminal fall speed, and vertical air velocity that are
C  derived from the complete variational column solution.  A Monte Carlo
C  approach is used.  This subroutine should be called with the same
C  input as that used in COMPLETE_VAR_COL_SOLN only if SUCCESS was
C  returned as .TRUE. from COMPLETE_VAR_COL_SOLN.  This subroutine will
C  stop execution if it is run with data for which SUCCESS was returned
C  as .FALSE. from COMPLETE_VAR_COL_SOLN.

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
C  equal BADDATA.  The only use made of DIV here is to indicate whether
C  or not a value of divergence is missing.

C  SD_DIV is a one-dimensional real array with elements indexed from 0
C  to N.  SD_DIV(I) specifies the standard deviation of DIV(I) (s**-1).
C  If it is missing, it should equal BADDATA.  SD_DIV(I) must not be
C  missing for DIV(I) to be usable.

C  WP is a one-dimensional real array with elements indexed from 0 to N.
C  WP(I) specifies the mean echo-power-weighted vertical hydrometeor
C  velocity (m s**-1) at the Ith level above the bottom of the column.
C  If it is missing, it should equal BADDATA.  Vertical hydrometeor
C  velocity is positive for upward motion.  The only use made of WP here
C  is to indicate whether or not a value of hydrometeor vertical
C  velocity is missing.

C  SD_WP is a one-dimensional real array with elements indexed from 0 to
C  N.  SD_WP(I) specifies the standard deviation of WP(I) (m s**-1).  If
C  it is missing, it should equal BADDATA.  SD_WP(I) must not be missing
C  for WP(I) to be usable.

C  F is a one-dimensional real array with elements indexed from 0 to N.
C  F(I) specifies the mean echo-power-weighted hydrometeor terminal fall
C  speed (m s**-1) at the Ith level above the bottom of the column.  If
C  it is missing, it should equal BADDATA.  Hydrometeor terminal fall
C  speed is greater than or equal to zero.  The only use made of F here
C  is to indicate whether or not a value of hydrometeor terminal fall
C  speed is missing.

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

C  SD_DIV_ADJ is a one-dimensional real array with elements indexed from
C  0 to N.  SD_DIV_ADJ(I) returns the standard deviation of the adjusted
C  horizontal divergence(s**-1) at the Ith level above the bottom of the
C  column.  If it is not calculated, it is returned as BADDATA.

C  SD_WP_ADJ is a one-dimensional real array with elements indexed from
C  0 to N.  SD_WP_ADJ(I) returns the standard deviation of the adjusted
C  mean echo-power-weighted vertical hydrometeor velocity (m s**-1) at
C  the Ith level above the bottom of the column. If it is not
C  calculated, it is returned as BADDATA.

C  SD_F_ADJ is a one-dimensional real array with elements indexed from 0
C  to N.  F_ADJ(I) returns the standard deviation of adjusted mean
C  echo-power-weighted terminal fall speed (m s**-1) at the Ith level
C  above the bottom of the column.  If it is not calculated, it is
C  returned as BADDATA.

C  SD_W is a one-dimensional real array with elements indexed from 0 to
C  N.  W(I) returns the standard deviation of vertical air velocity (m
C  s**-1) at the Ith level above the bottom of the column.  If it is not
C  calculated, it is returned as BADDATA.

      IMPLICIT NONE
      INCLUDE 'tmrlib.inc'
      INTEGER INSEED
      PARAMETER(INSEED=88)
      LOGICAL SUCCESS
      INTEGER N,J,I
      INTEGER N_DIV_ADJ(0:N),N_WP_ADJ(0:N),N_F_ADJ(0:N),N_W(0:N)
      REAL Z1,TV1,Z2,TV2,Z3,RHO3,BADDATA
      REAL DIV(0:N),WP(0:N),F(0:N),Z(0:N),W_CONSTRAINT(0:N),
     $SD_W_CONSTRAINT(0:N),SD_DIV(0:N),SD_WP(0:N),SD_F(0:N),
     $DIV_TRIAL(0:N),F_TRIAL(0:N),WP_TRIAL(0:N),DIV_ADJ_TRIAL(0:N),
     $WP_ADJ_TRIAL(0:N),F_ADJ_TRIAL(0:N),W_TRIAL(0:N),TOT_SQ_ERR_W(0:N),
     $TOT_SQ_ERR_DIV_ADJ(0:N),TOT_SQ_ERR_WP_ADJ(0:N),
     $TOT_SQ_ERR_F_ADJ(0:N),SD_W(0:N),SD_DIV_ADJ(0:N),SD_WP_ADJ(0:N),
     $SD_F_ADJ(0:N),W_CONSTRAINT_TRIAL(0:N)
      
C  Initialize the total squared error and a counter for each variable
C  and at each level to zero.
      DO I=0,N
         TOT_SQ_ERR_W(I)=0.
         TOT_SQ_ERR_DIV_ADJ(I)=0.
         TOT_SQ_ERR_WP_ADJ(I)=0.
         TOT_SQ_ERR_F_ADJ(I)=0.
         N_DIV_ADJ(I)=0
         N_WP_ADJ(I)=0
         N_F_ADJ(I)=0
         N_W(I)=0
      ENDDO

C  Perform many trials.
      DO J=1,COMPLETE_VAR_COL_SOLN_NTRIALS

C  At each level, generate a random value for divergence, vertical
C  hydrometeor velocity, and hydrometeor terminal fall speed.  These
C  values are drawn from normal distributions with the input standard
C  deviations.  A random value is generated only if the original datum
C  is not missing.  These random values are considered, without loss of
C  generality, as perturbations on a state having no divergence or
C  vertical motion.
         DO I=0,N
            IF(DIV(I).NE.BADDATA.AND.
     $      SD_DIV(I).NE.BADDATA)THEN
               CALL RNGN(INSEED,0.,SD_DIV(I),DIV_TRIAL(I))
            ELSE
               DIV_TRIAL(I)=BADDATA
            ENDIF
            IF(WP(I).NE.BADDATA.AND.
     $      SD_WP(I).NE.BADDATA)THEN
               CALL RNGN(INSEED,0.,SD_WP(I),WP_TRIAL(I))
            ELSE
               WP_TRIAL(I)=BADDATA
            ENDIF
            IF(F(I).NE.BADDATA.AND.
     $      SD_F(I).NE.BADDATA)THEN
               CALL RNGN(INSEED,0.,SD_F(I),F_TRIAL(I))
            ELSE
               F_TRIAL(I)=BADDATA
            ENDIF
            IF(W_CONSTRAINT(I).NE.BADDATA)THEN
               CALL RNGN(INSEED,0.,SD_W_CONSTRAINT(I),
     $         W_CONSTRAINT_TRIAL(I))
            ELSE
               W_CONSTRAINT_TRIAL(I)=BADDATA
            ENDIF
         ENDDO

C  Calculate the adjusted horizontal divergence, the adjusted vertical
C  hydrometeor velocity, the adjusted hydrometeor terminal fall speed,
C  and the vertical air velocity from the random input.
         CALL COMPLETE_VAR_COL_SOLN(Z1,TV1,Z2,TV2,Z3,RHO3,N,Z,DIV_TRIAL,
     $   SD_DIV,WP_TRIAL,SD_WP,F_TRIAL,SD_F,W_CONSTRAINT,BADDATA,
     $   DIV_ADJ_TRIAL,WP_ADJ_TRIAL,F_ADJ_TRIAL,W_TRIAL,SUCCESS)
         IF(SUCCESS.EQV..FALSE.)THEN
            WRITE(TMRLIB_MESSAGE_UNIT,*)'COMPLETE_VAR_COL_SOLN_SD:  ',
     $      'TRIAL FAILED.'
            STOP
         ENDIF

C  Calculate the total squared errors for each level and variable.  The
C  errors are equal to the variable values, since the state perturbed
C  with errors was one with no divergence or vertical motion.
         DO I=0,N
            IF(DIV_ADJ_TRIAL(I).NE.BADDATA)THEN
               TOT_SQ_ERR_DIV_ADJ(I)=TOT_SQ_ERR_DIV_ADJ(I)+
     $         DIV_ADJ_TRIAL(I)**2
               N_DIV_ADJ(I)=N_DIV_ADJ(I)+1
            ENDIF
            IF(WP_ADJ_TRIAL(I).NE.BADDATA)THEN
               TOT_SQ_ERR_WP_ADJ(I)=TOT_SQ_ERR_WP_ADJ(I)+
     $         WP_ADJ_TRIAL(I)**2
               N_WP_ADJ(I)=N_WP_ADJ(I)+1
            ENDIF
            IF(F_ADJ_TRIAL(I).NE.BADDATA)THEN
               TOT_SQ_ERR_F_ADJ(I)=TOT_SQ_ERR_F_ADJ(I)+F_ADJ_TRIAL(I)**2
               N_F_ADJ(I)=N_F_ADJ(I)+1
            ENDIF
            IF(W_TRIAL(I).NE.BADDATA)THEN
               TOT_SQ_ERR_W(I)=TOT_SQ_ERR_W(I)+W_TRIAL(I)**2
               N_W(I)=N_W(I)+1
            ENDIF
         ENDDO
      ENDDO

C  Calculate the standard deviation for each level and variable.
      DO I=0,N
         IF(N_DIV_ADJ(I).NE.0)THEN
            SD_DIV_ADJ(I)=SQRT(TOT_SQ_ERR_DIV_ADJ(I)/
     $      FLOAT(N_DIV_ADJ(I)))
         ELSE
            SD_DIV_ADJ(I)=BADDATA
         ENDIF
         IF(N_WP_ADJ(I).NE.0)THEN
            SD_WP_ADJ(I)=SQRT(TOT_SQ_ERR_WP_ADJ(I)/FLOAT(N_WP_ADJ(I)))
         ELSE
            SD_WP_ADJ(I)=BADDATA
         ENDIF
         IF(N_F_ADJ(I).NE.0)THEN
            SD_F_ADJ(I)=SQRT(TOT_SQ_ERR_F_ADJ(I)/FLOAT(N_F_ADJ(I)))
         ELSE
            SD_F_ADJ(I)=BADDATA
         ENDIF
         IF(N_W(I).NE.0)THEN
            SD_W(I)=SQRT(TOT_SQ_ERR_W(I)/FLOAT(N_W(I)))
         ELSE
            SD_W(I)=BADDATA
         ENDIF
      ENDDO

C  Done.
      RETURN
      END
