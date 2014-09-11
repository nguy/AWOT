      SUBROUTINE MULTI_STEP_BINOMIAL_FILTER_WEIGHTS(N_STEPS,
     $WT_OUT_MAX_REACH,WT_OUT,REACH_OUT)

C  Thomas Matejka NOAA/NSSL 22 April 1994

C  This subroutine calculates the effective one-step filter weights for
C  a multiple application of a three-point binomial filter.

C  Input:

C  N_STEPS is an integer variable that specifies the number of binomial
C  filter steps to use.

C  WT_OUT_MAX_REACH is an integer variable that specifies that WT_OUT is
C  dimensioned from -WT_OUT_MAX_REACH to WT_OUT_MAX_REACH in the calling
C  program.

C  Output:

C  WT_OUT is a one-dimensional real array.  WT_OUT(I) returns the
C  effective one-step filter weight at grid point I relative to the
C  point of application.

C  REACH_OUT is an integer variable that returns that non-zero weights
C  in WT_OUT extend from -REACH_OUT to REACH_OUT grid points relative to
C  the point of application.

      IMPLICIT NONE
      INTEGER N_STEPS,WT_OUT_MAX_REACH,REACH_OUT,I_STEP,ARRAY_REACH
      REAL WT_OUT(-WT_OUT_MAX_REACH:WT_OUT_MAX_REACH)
      REAL WT_ARRAY(N_STEPS+1,-1:1)

C  Record the reach of the filter weight array.
      IF(N_STEPS.GE.1)THEN
         ARRAY_REACH=1
      ELSE
         ARRAY_REACH=0
      ENDIF

C  Construct the multiple-step binomial filter weight array.
      IF(N_STEPS.GE.1)THEN
         DO 1 I_STEP=1,N_STEPS
            WT_ARRAY(I_STEP,-1)=0.25
            WT_ARRAY(I_STEP,0)=0.5
            WT_ARRAY(I_STEP,1)=0.25
1        CONTINUE
      ENDIF

C  Calculate the equivalent one-step filter weights.
      CALL FILTER_WEIGHT_CALCULATOR(WT_ARRAY,N_STEPS+1,1,N_STEPS,
     $ARRAY_REACH,WT_OUT_MAX_REACH,WT_OUT,REACH_OUT)

C  Done.
      RETURN
      END
