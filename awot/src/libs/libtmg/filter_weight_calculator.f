      SUBROUTINE FILTER_WEIGHT_CALCULATOR(WT_ARRAY,WT_ARRAY_DIM_1,
     $WT_ARRAY_MAX_REACH,N_STEPS,ARRAY_REACH,WT_OUT_MAX_REACH,WT_OUT,
     $REACH_OUT)

C  Thomas Matejka NOAA/NSSL 22 April 1994

C  This subroutine calculates the effective one-step filter weight that
C  is equivalent to several successive filter steps.  The filter is
C  one-dimensional.

C  Input:

C  WT_ARRAY is a two-dimensional real array.  WT_ARRAY(I,J) specifies
C  the weight at grid point J for the Ith filter step.

C  WT_ARRAY_DIM_1 is an integer variable that specifies the first
C  dimension of WT_ARRAY in the calling program.

C  N_STEPS is an integer variable that specifies how many filter steps
C  are specified in WT_ARRAY.

C  ARRAY_REACH is an integer variable that specifies that weights in
C  WT_ARRAY run from -ARRAY_REACH to ARRAY_REACH grid points relative to
C  the point of application.

C  WT_OUT_MAX_REACH is an integer variable that specifies that WT_OUT is
C  dimensioned from -WT_OUT_MAX_REACH to WT_OUT_MAX_REACH in the calling
C  program.

C  Output:

C  WT_OUT is a one-dimensional real array.  WT_OUT(J) returns the
C  effective one-step filter weight at grid point J relative to the
C  point of application.

C  REACH_OUT is an integer variable that returns that non-zero weights
C  in WT_OUT run from -REACH_OUT to REACH_OUT grid points relative to
C  the point of application.

      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      INTEGER WT_ARRAY_DIM_1,N_STEPS,ARRAY_REACH,WT_OUT_MAX_REACH,
     $REACH_OUT,I_STEP,J,K,L,WT_ARRAY_MAX_REACH
      REAL WT_ARRAY(WT_ARRAY_DIM_1,
     $-WT_ARRAY_MAX_REACH:WT_ARRAY_MAX_REACH)
      REAL WT_OUT(-WT_OUT_MAX_REACH:WT_OUT_MAX_REACH)
      REAL WT_SUM(-WT_OUT_MAX_REACH:WT_OUT_MAX_REACH)

C  Initialize the effective one-pass filter weights after 0 passes.
      DO 1 K=-WT_OUT_MAX_REACH,WT_OUT_MAX_REACH
         IF(K.EQ.0)THEN
            WT_OUT(K)=1.
         ELSE
            WT_OUT(K)=0.
         ENDIF
1     CONTINUE
      REACH_OUT=0

C  Apply the filter weights step by step.
      IF(N_STEPS.GE.1)THEN
         DO 2 I_STEP=1,N_STEPS
            DO 3 K=-WT_OUT_MAX_REACH,WT_OUT_MAX_REACH
               WT_SUM(K)=0.
3           CONTINUE
            DO 4 K=-WT_OUT_MAX_REACH,WT_OUT_MAX_REACH
               IF(WT_OUT(K).NE.0.)THEN
                  DO 5 J=-ARRAY_REACH,ARRAY_REACH
                     IF(WT_ARRAY(I_STEP,J).NE.0.)THEN
                        L=K+J
                        IF(IABS(L).GT.WT_OUT_MAX_REACH)THEN
                           WRITE(TMMLIB_MESSAGE_UNIT,*)
     $                     'FILTER_WEIGHT_CALCULATOR:  ',
     $                     'WT_OUT_MAX_REACH IS TOO SMALL.'
                           STOP
                        ENDIF
                        WT_SUM(L)=WT_SUM(L)+WT_OUT(K)*WT_ARRAY(I_STEP,J)
                        REACH_OUT=MAX0(REACH_OUT,IABS(L))
                     ENDIF
5                 CONTINUE
               ENDIF
4           CONTINUE
            DO 6 K=-REACH_OUT,REACH_OUT
               WT_OUT(K)=WT_SUM(K)
6           CONTINUE
2        CONTINUE
      ENDIF

C  Done.
      RETURN
      END
