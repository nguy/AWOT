      SUBROUTINE REGRESSION_VAR_2(NPARAMS,NDATA,X,X_DIM_1,Y,VAR_Y,
     $ESTIMATE_SSQ,SSQ_IN,
     $BADDATA,
     $CHECK_SINGULAR,SINGULAR_THRESHOLD,
     $C,VAR_C)

C  Thomas Matejka NOAA/NSSL 24 July 2002

C  This subroutine performs a linear least squares regression involving
C  any number of independent variables. The response variable data are
C  allowed to have different variances but to be uncorrelated. To
C  correct for the variable variances, the regression is performed with
C  transformed data. The variances of the parameters are also estimated.

C  Input:

C  NPARAMS (integer) specifies the number of predictor variables
C  (including the constant term if used and counting separately terms
C  involving different powers of the same variable) in the regression
C  model. Each predictor variable has a parameter to be estimated.

C  NDATA (integer) specifies the number of data used to perform the
C  regression.

C  X (2D real array 1:X_DIM_1,1:NPARAMS). X(I,J) specifies the value of
C  the Jth predictor variable for the Ith datum.

C  X_DIM_1 (integer) specifies the first dimension of X in the calling
C  program.

C  Y (1D real array 1:NDATA). Y(I) specifies the value of the response
C  variable for the Ith datum.

C  VAR_Y (1D real array 1:NDATA). VAR_Y(I) specifies a value
C  proportional to the variance of Y(I). VAR_Y(I) must be greater than
C  0..

C  ESTIMATE_SSQ (logical). .TRUE. indicates that the variance of the
C  response variable data will be estimated from the scatter about the
C  regression under the assumption that the regression model is correct.
C  .FALSE. indicates that the variance of the response variable data is
C  specified using SSQ_IN as described.

C  SSQ_IN (real) specifies the value that, when multiplied by VAR_Y(I),
C  gives the variance of Y(I). It is relevant only when ESTIMATE_SSQ is
C  .FALSE..

C  BADDATA (real) indicates a missing value as described.

C  CHECK_SINGULAR (logical). .TRUE. indicates that matrices to be
C  inverted are checked for near-singularity. .FALSE. indicates that
C  they are not checked.

C  SINGULAR_THRESHOLD (real) controls the definition of
C  near-singularity. SINGULAR_THRESHOLD should be greater than or equal
C  to 0.. The smaller SINGULAR_THRESHOLD is, the more close to exactly
C  singular matrices must be to be considered nearly singular. It is
C  relevant only when CHECK_SINGULAR is .TRUE..

C  Output:

C  C (1D real array 1:NPARAMS). C(I) returns an estimate of the
C  parameter for the Ith predictor variable. If it is missing, it is
C  returned as BADDATA.

C  VAR_C (1D real array 1:NPARAMS). VAR_C(I) returns an estimate of the
C  variance of C(I). If it is missing, it is returned as BADDATA.

      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      LOGICAL::CHECK_SINGULAR,INVERSE_SUCCESS,ESTIMATE_SSQ
      INTEGER::NDATA,NPARAMS,IDATUM,IPARAM,JPARAM,X_DIM_1
      REAL::BADDATA,SSQ,SSTOT,SSREG,SSRES,SINGULAR_THRESHOLD,SSQ_IN
      REAL,DIMENSION(NPARAMS)::B,C,VAR_C
      REAL,DIMENSION(NDATA)::Y,VAR_Y
      REAL,DIMENSION(X_DIM_1,NPARAMS)::X
      REAL,DIMENSION(NPARAMS,NPARAMS)::A,A_INV

C  There must be at least as many data as parameters.
      IF(NDATA.LT.NPARAMS)THEN
         DO IPARAM=1,NPARAMS
            C(IPARAM)=BADDATA
            VAR_C(IPARAM)=BADDATA
         ENDDO
         RETURN
      ENDIF

C  Calculate the matrix of the sum of products of the transformed
C  predictor variable data.
      DO IPARAM=1,NPARAMS
         DO JPARAM=1,IPARAM
            A(IPARAM,JPARAM)=0.
            DO IDATUM=1,NDATA
               A(IPARAM,JPARAM)=A(IPARAM,JPARAM)+X(IDATUM,IPARAM)*
     $         X(IDATUM,JPARAM)/VAR_Y(IDATUM)
            ENDDO
         ENDDO
      ENDDO
      IF(NPARAMS.GE.2)THEN
         DO IPARAM=2,NPARAMS
            DO JPARAM=1,IPARAM-1
               A(JPARAM,IPARAM)=A(IPARAM,JPARAM)
            ENDDO
         ENDDO
      ENDIF

C  Calculate the vector of the sum of products of the transformed
C  predictor variable data and the transformed response variable data.
      DO IPARAM=1,NPARAMS
         B(IPARAM)=0.
         DO IDATUM=1,NDATA
            B(IPARAM)=B(IPARAM)+X(IDATUM,IPARAM)*Y(IDATUM)/VAR_Y(IDATUM)
         ENDDO
      ENDDO

C  Find the inverse of the matrix of the sum of products of the
C  transformed predictor variable data.
      CALL MAT_INV(A,NPARAMS,NPARAMS,A_INV,NPARAMS,CHECK_SINGULAR,
     $SINGULAR_THRESHOLD,0,INVERSE_SUCCESS)
      IF(.NOT.INVERSE_SUCCESS)THEN
         DO IPARAM=1,NPARAMS
            C(IPARAM)=BADDATA
            VAR_C(IPARAM)=BADDATA
         ENDDO
         RETURN
      ENDIF

C  Calculate the parameters.
      DO IPARAM=1,NPARAMS
         C(IPARAM)=0.
         DO JPARAM=1,NPARAMS
            C(IPARAM)=C(IPARAM)+A_INV(IPARAM,JPARAM)*B(JPARAM)
         ENDDO
      ENDDO

C  Estimate the variances of the parameters.

C  The variance of the transformed response variable data will be
C  estimated from the scatter about the regression.
      IF(ESTIMATE_SSQ)THEN

C  There must be more data than parameters.
         IF(NDATA.LE.NPARAMS)THEN
            DO IPARAM=1,NPARAMS
               VAR_C(IPARAM)=BADDATA
            ENDDO
            RETURN
         ENDIF

C  Calculate the total sum of squares of the transformed response
C  variable data.
         SSTOT=0.
         DO IDATUM=1,NDATA
            SSTOT=SSTOT+Y(IDATUM)**2/VAR_Y(IDATUM)
         ENDDO

C  Calculate the sum of squares of the transformed response variable
C  data due to the regression.
         SSREG=0.
         DO IPARAM=1,NPARAMS
            SSREG=SSREG+C(IPARAM)*B(IPARAM)
         ENDDO

C  Calculate the residual sum of squares of the response variable data.
         SSRES=SSTOT-SSREG
         IF(SSRES.LT.0.)THEN
            SSRES=0.
         ENDIF

C  Estimate the variance of the transformed response variable data.
         SSQ=SSRES/FLOAT(NDATA-NPARAMS)

C  Use the specified variance of the transformed response variable data.
      ELSE
         SSQ=SSQ_IN
      ENDIF

C  Estimate the variances of the parameters.
      DO IPARAM=1,NPARAMS
         VAR_C(IPARAM)=SSQ*A_INV(IPARAM,IPARAM)
      ENDDO

      END SUBROUTINE REGRESSION_VAR_2
