      SUBROUTINE LLS(NPARAMS,NDATA,X,X_DIM_1,Y,GROUPED,GROUP_SIZE_IN,
     $ESTIMATE_SSQ,SSQ_IN,CHECK_SINGULAR,SINGULAR_THRESHOLD,WRITE_MODE,
     $A_INV_DIM_1,EFFECTIVE_NDATA,C,VAR_C,SSQ,A_INV,SUCCESS)

C  Thomas Matejka NOAA/NSSL 13 September 1995

C  This subroutine performs a linear least squares regression involving
C  any number of independent variables.  The response variable data are
C  all assumed to have the same variance about a potential regression
C  and to be uncorrelated.  Data points may be grouped so that each
C  group affects the regression as a single point.  The variance of the
C  response variable data about the regression may be either known a
C  priori or estimated from the regression under the assumption that the
C  regression model is correct.

C  Input:

C  NPARAMS is an integer variable that specifies the number of predictor
C  variable (including the constant term if used and counting separately
C  terms involving different powers of the same variable) in the
C  regression model.  Each predictor variable has a parameter to be
C  estimated.

C  NDATA is an integer variable that specifies the number of data used
C  to develop the regression.

C  X is a two-dimensional real array.  X(I,J) specifies the value of the
C  Jth predictor variable for the Ith datum.

C  X_DIM_1 is an integer variable that specifies the first dimension of
C  X in the calling program.

C  Y is a one-dimensional real array.  Y(I) specifies the value of the
C  response variable for the Ith datum.

C  GROUPED is a logical variable.  If GROUPED is .TRUE., data points are
C  grouped as specified in GROUP_SIZE_IN so that each group contributes
C  to the regression as one data point.

C  GROUP_SIZE_IN is a one-dimensional integer array that is relevant
C  only when GROUPED is .TRUE..  GROUP_SIZE(I) specifies the total
C  number of data points in the group to which the Ith data point
C  belongs.  Data points in a group contribute to the regression as if
C  they are one point.

C  ESTIMATE_SSQ is a logical variable.  If ESTIMATE_SSQ is .TRUE., then
C  the variance of the response variable data about the regression will
C  be estimated from the data under the assumption that the regression
C  model is correct.  If ESTIMATE_SSQ is .FALSE., then the variance of
C  the response variable data about a potential regression is specified
C  using SSQ_IN as described.

C  SSQ_IN is a real variable that is relevant only when ESTIMATE_SSQ is
C  .FALSE..  SSQ_IN should specify the variance of all the response
C  variable data about a potential regression.

C  CHECK_SINGULAR is a logical variable.  Matrices to be inverted are
C  checked for near-singularity if and only if CHECK_SINGULAR is .TRUE.

C  SINGULAR_THRESHOLD is a real variable that controls the definition of
C  nearly singular.  SINGULAR_THRESHOLD should be greater than or equal
C  to 0..  The smaller SINGULAR_THRESHOLD is, the more close to exactly
C  singular matrices must be to be considered nearly singular.  It is
C  relevant only when CHECK_SINGULAR is .TRUE..

C  WRITE_MODE is an integer variable that specifies whether matrices are
C  written at various steps of the row reduction.  If WRITE_MODE is 0,
C  nothing is written.  If WRITE_MODE is 1, the initial and row-reduced
C  matrices are written.  If WRITE_MODE is 2, the initial, final, and
C  all intermediate matrices during the row reduction are written.

C  A_INV_DIM_1 is an integer variable that specifies the first dimension
C  of A_INV in the calling program.  It should be greater than or equal
C  to NPARAMS.

C  Output:

C  EFFECTIVE_NDATA is an integer variable that returns the effective
C  number of data points.  If GROUPED is .FALSE., it equals NDATA.  If
C  GROUPED is .TRUE., it equals the number of groups of data.

C  C is a one-dimensonal real array.  C(J) returns the estimate of the
C  parameter for the Jth predictor variable.

C  VAR_C is a one-dimensional real array.  VAR_C(J) returns the variance
C  of the estimate of the Jth parameter.  If ESTIMATE_SSQ is true, then
C  VAR_C is calculated under the assumption that the regression model is
C  correct.

C  SSQ is a real variable that returns the variance of all the response
C  variable data about the regression.  If ESTIMATE_SSQ is .TRUE., it is
C  estimated from the data under the assumption that the regression
C  model is correct.  If ESTIMATE_SSQ is .FALSE., then it equals SSQ_IN.

C  A_INV is a two-dimensional real array.  A_INV returns the inverse of
C  A.  A(I,J) is the sum of products of the Ith and Jth predictor
C  variable data.

C  SUCCESS is a logical variable that returns .TRUE. if and only if the
C  regression could be performed.

      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      LOGICAL SUCCESS,ESTIMATE_SSQ,INVERSE_SUCCESS,CHECK_SINGULAR,
     $GROUPED
      INTEGER NDATA,NPARAMS,IDATUM,IPARAM,JPARAM,EFFECTIVE_NDATA,
     $X_DIM_1,WRITE_MODE,A_INV_DIM_1
      INTEGER GROUP_SIZE_IN(NDATA)
      INTEGER GROUP_SIZE(NDATA)
      REAL SSQ,SSQ_IN,SSTOT,SSREG,SSRES,SUM,SINGULAR_THRESHOLD
      REAL C(NPARAMS),VAR_C(NPARAMS)
      REAL Y(NDATA)
      REAL B(NPARAMS)
      REAL X(X_DIM_1,NPARAMS)
      REAL A_INV(A_INV_DIM_1,NPARAMS)
      REAL A(NPARAMS,NPARAMS)

C  Check memory size.
      IF(A_INV_DIM_1.LT.NPARAMS)THEN
         WRITE(TMMLIB_MESSAGE_UNIT,*)'LLS:  A_INV_DIM_1 MUST BE ',
     $   'GREATER THAN OR EQUAL TO NPARAMS.'
         STOP
      ENDIF

C  Set the group sizes and calculate the number of groups.
      IF(GROUPED)THEN
         SUM=0.
         DO 1 IDATUM=1,NDATA
            GROUP_SIZE(IDATUM)=GROUP_SIZE_IN(IDATUM)
            SUM=SUM+1./FLOAT(GROUP_SIZE(IDATUM))
1        CONTINUE
         EFFECTIVE_NDATA=IFIX(SUM+0.5)
      ELSE
         DO 2 IDATUM=1,NDATA
            GROUP_SIZE(IDATUM)=1
2        CONTINUE
         EFFECTIVE_NDATA=NDATA
      ENDIF

C  Check whether there are enough data.  If SSQ is being estimated, need
C  at least one more datum than parameters.  Otherwise, need at least as
C  many data as parameters.
      IF(ESTIMATE_SSQ)THEN
         IF(EFFECTIVE_NDATA.LE.NPARAMS)THEN
            SUCCESS=.FALSE.
            RETURN
         ENDIF
      ELSE
         IF(EFFECTIVE_NDATA.LT.NPARAMS)THEN
            SUCCESS=.FALSE.
            RETURN
         ENDIF
      ENDIF

C  Calculate A, the matrix of the sum of products of the predictor
C  variable data.
      DO 3 IPARAM=1,NPARAMS
         DO 4 JPARAM=1,IPARAM
            A(IPARAM,JPARAM)=0.
            DO 5 IDATUM=1,NDATA
               A(IPARAM,JPARAM)=A(IPARAM,JPARAM)+X(IDATUM,IPARAM)*
     $         X(IDATUM,JPARAM)/FLOAT(GROUP_SIZE(IDATUM))
5           CONTINUE
4        CONTINUE
3     CONTINUE
      IF(NPARAMS.GE.2)THEN
         DO 6 IPARAM=2,NPARAMS
            DO 7 JPARAM=1,IPARAM-1
               A(JPARAM,IPARAM)=A(IPARAM,JPARAM)
7           CONTINUE
6        CONTINUE
      ENDIF

C  Calculate the vector of the sum of products of the predictor variable
C  data and the response variable data.
      DO 8 IPARAM=1,NPARAMS
         B(IPARAM)=0.
         DO 9 IDATUM=1,NDATA
            B(IPARAM)=B(IPARAM)+X(IDATUM,IPARAM)*Y(IDATUM)/
     $      FLOAT(GROUP_SIZE(IDATUM))
9        CONTINUE
8     CONTINUE

C  Find the inverse of A.
      CALL MAT_INV(A,NPARAMS,NPARAMS,A_INV,A_INV_DIM_1,CHECK_SINGULAR,
     $SINGULAR_THRESHOLD,WRITE_MODE,INVERSE_SUCCESS)
      IF(.NOT.INVERSE_SUCCESS)THEN
         SUCCESS=.FALSE.
         RETURN
      ENDIF

C  Calculate the estimated parameters.
      DO 10 IPARAM=1,NPARAMS
         C(IPARAM)=0.
         DO 11 JPARAM=1,NPARAMS
            C(IPARAM)=C(IPARAM)+A_INV(IPARAM,JPARAM)*B(JPARAM)
11       CONTINUE
10    CONTINUE

C  Calculate the total sum of squares of the response variable data.
      SSTOT=0.
      DO 12 IDATUM=1,NDATA
         SSTOT=SSTOT+Y(IDATUM)**2/FLOAT(GROUP_SIZE(IDATUM))
12    CONTINUE

C  Calculate the sum of squares of the response variable data due to the
C  regression.
      SSREG=0.
      DO 13 IPARAM=1,NPARAMS
         SSREG=SSREG+C(IPARAM)*B(IPARAM)
13    CONTINUE

C  Calculate the residual sum of squares of the response variable data.
      SSRES=SSTOT-SSREG
      IF(SSRES.LT.0.)THEN
         SSRES=0.
      ENDIF

C  If requested, estimate the variance of the response variable data
C  about the regression.  Otherwise, use the value specified a priori.
      IF(ESTIMATE_SSQ)THEN
         SSQ=SSRES/FLOAT(EFFECTIVE_NDATA-NPARAMS)
      ELSE
         SSQ=SSQ_IN
      ENDIF

C  Calculate the variances of the estimated parameters.
      DO 14 IPARAM=1,NPARAMS
         VAR_C(IPARAM)=SSQ*A_INV(IPARAM,IPARAM)
14    CONTINUE

C  Done.
      SUCCESS=.TRUE.
      RETURN
      END
