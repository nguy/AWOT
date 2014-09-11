C  tmrlib.inc

C  Thomas Matejka NOAA/NSSL 13 September 1995

C  This INCLUDE file defines constants and memory sizes for TMRLIB.

C  Logical device number for TMRLIB messages.
      INTEGER TMRLIB_MESSAGE_UNIT
      PARAMETER(TMRLIB_MESSAGE_UNIT=7)

C  Threshold for determining nearly singular matrices for the DOP*
C  family of routines.
      REAL DOP_SINGULAR_THRESHOLD
      PARAMETER(DOP_SINGULAR_THRESHOLD=0.001)

C  Threshold for determining nearly singular matrices for the ESTIMATE*
C  family of routines.
      REAL ESTIMATE_SINGULAR_THRESHOLD
      PARAMETER(ESTIMATE_SINGULAR_THRESHOLD=0.001)

C  Whether matrices are written at various steps of the row reduction.
C  If MAT_WRITE_MODE is 0, nothing is written.  If MAT_WRITE_MODE is 1,
C  the initial and row-reduced matrices are written.  If MAT_WRITE_MODE
C  is 2, the initial, final, and all intermediate matrices during the
C  row reduction are written.
      INTEGER MAT_WRITE_MODE
      PARAMETER(MAT_WRITE_MODE=0)

C  Number of trials to use when computing, by Monte Carlo method, the
C  standard deviations of the results from the complete variational
C  column solution.
      INTEGER COMPLETE_VAR_COL_SOLN_NTRIALS
      PARAMETER(COMPLETE_VAR_COL_SOLN_NTRIALS=100)

C  Maximum number of parameters in a sequence step in W_DRIVER.
      INTEGER MAX_W_DRIVER_PARAMETERS
      PARAMETER(MAX_W_DRIVER_PARAMETERS=5)

C  Maximum character string length.
      INTEGER MAX_STRING
      PARAMETER(MAX_STRING=200)
