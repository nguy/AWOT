      SUBROUTINE MULTI_STEP_LEISE_FILTER_WEIGHTS(N_STEPS,
     $WT_OUT_MAX_REACH,WT_OUT,REACH_OUT)

C  Thomas Matejka NOAA/NSSL 22 April 1994

C  This subroutine calculates the effective one-step filter weights for
C  a multiple application of a Leise filter.

C  Input:

C  N_STEPS is an integer variable that specifies the number of Leise
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
      INTEGER N_STEPS,WT_OUT_MAX_REACH,REACH_OUT,I_STEP,J_STEP,
     $ARRAY_REACH,K,N_STEPS_ACTUAL
      REAL WT_OUT(-WT_OUT_MAX_REACH:WT_OUT_MAX_REACH)
      REAL WT_ARRAY(2*N_STEPS+1,-(2**N_STEPS):2**N_STEPS)

C  Calculate the reach of the filter weight array.
      IF(N_STEPS.GE.1)THEN
         ARRAY_REACH=2**N_STEPS
      ELSE
         ARRAY_REACH=0
      ENDIF

C  Calculate the actual number of filter steps.
      IF(N_STEPS.GE.1)THEN
         N_STEPS_ACTUAL=2*N_STEPS-1
      ELSE
         N_STEPS_ACTUAL=0
      ENDIF

C  Initialize all the relevant filter weight array elements to zero.
      IF(N_STEPS_ACTUAL.GE.1)THEN
         DO 1 J_STEP=1,N_STEPS_ACTUAL
            DO 2 K=-ARRAY_REACH,ARRAY_REACH
               WT_ARRAY(J_STEP,K)=0.
2           CONTINUE
1        CONTINUE
      ENDIF

C  Construct the multiple-step Leise filter weight array.
      IF(N_STEPS.GE.1)THEN
         J_STEP=0
         DO 3 I_STEP=1,N_STEPS
            K=2**(I_STEP-1)
            J_STEP=J_STEP+1
            WT_ARRAY(J_STEP,-2*K)=-0.0625
            WT_ARRAY(J_STEP,-K)=0.25
            WT_ARRAY(J_STEP,0)=0.625
            WT_ARRAY(J_STEP,K)=0.25
            WT_ARRAY(J_STEP,2*K)=-0.0625
            IF(I_STEP.NE.N_STEPS)THEN
               J_STEP=J_STEP+1
               WT_ARRAY(J_STEP,-2*K)=-0.0625
               WT_ARRAY(J_STEP,-K)=0.25
               WT_ARRAY(J_STEP,0)=0.625
               WT_ARRAY(J_STEP,K)=0.25
               WT_ARRAY(J_STEP,2*K)=-0.0625
            ENDIF
3        CONTINUE
      ENDIF

C  Calculate the equivalent one-step filter weights.
      CALL FILTER_WEIGHT_CALCULATOR(WT_ARRAY,2*N_STEPS+1,2**N_STEPS,
     $N_STEPS_ACTUAL,ARRAY_REACH,WT_OUT_MAX_REACH,WT_OUT,REACH_OUT)

C  Done.
      RETURN
      END
