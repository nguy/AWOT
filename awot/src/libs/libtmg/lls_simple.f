      SUBROUTINE LLS_SIMPLE(NPARAMS,NDATA,X,X_DIM_1,Y,CHECK_SINGULAR,
     $SINGULAR_THRESHOLD,C,VAR_C,SUCCESS)

C  Thomas Matejka NOAA/NSSL 13 September 1995

C  This subroutine performs a linear least squares regression involving
C  any number of independent variables.  The response variable data are
C  all assumed to have the same variance about a potential regression
C  and to be uncorrelated.  The variance of the response variable data
C  about the regression is estimated from the regression under the
C  assumption that the regression model is correct.

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

C  CHECK_SINGULAR is a logical variable.  Matrices to be inverted are
C  checked for near-singularity if and only if CHECK_SINGULAR is .TRUE.

C  SINGULAR_THRESHOLD is a real variable that controls the definition of
C  nearly singular.  SINGULAR_THRESHOLD should be greater than or equal
C  to 0..  The smaller SINGULAR_THRESHOLD is, the more close to exactly
C  singular matrices must be to be considered nearly singular.  It is
C  relevant only when CHECK_SINGULAR is .TRUE..

C  Output:

C  C is a one-dimensonal real array.  C(J) returns the estimate of the
C  parameter for the Jth predictor variable.

C  VAR_C is a one-dimensional real array.  VAR_C(J) returns the variance
C  of the estimate of the Jth parameter.  If ESTIMATE_SSQ is true, then
C  VAR_C is calculated under the assumption that the regression model is
C  correct.

C  SUCCESS is a logical variable that returns .TRUE. if and only if the
C  regression could be performed.

      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      LOGICAL SUCCESS,CHECK_SINGULAR,INVERSE_SUCCESS
      INTEGER NDATA,NPARAMS,IDATUM,IPARAM,JPARAM,X_DIM_1
      REAL SSQ,SSTOT,SSREG,SSRES,SINGULAR_THRESHOLD
      REAL C(NPARAMS),VAR_C(NPARAMS)
      REAL Y(NDATA)
      REAL B(NPARAMS)
      REAL X(X_DIM_1,NPARAMS)
      REAL A(NPARAMS,NPARAMS),A_INV(NPARAMS,NPARAMS)

C  Check whether there are enough data.  Need at least one more datum
C  than parameters.
      IF(NDATA.LE.NPARAMS)THEN
         SUCCESS=.FALSE.
         RETURN
      ENDIF

C  Calculate A, the matrix of the sum of products of the predictor
C  variable data.
      DO IPARAM=1,NPARAMS
         DO JPARAM=1,IPARAM
            A(IPARAM,JPARAM)=0.
            DO IDATUM=1,NDATA
               A(IPARAM,JPARAM)=A(IPARAM,JPARAM)+X(IDATUM,IPARAM)*
     $         X(IDATUM,JPARAM)
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

C  Calculate the vector of the sum of products of the predictor variable
C  data and the response variable data.
      DO IPARAM=1,NPARAMS
         B(IPARAM)=0.
         DO IDATUM=1,NDATA
            B(IPARAM)=B(IPARAM)+X(IDATUM,IPARAM)*Y(IDATUM)
         ENDDO
      ENDDO

C  Find the inverse of A.
      CALL MAT_INV(A,NPARAMS,NPARAMS,A_INV,NPARAMS,CHECK_SINGULAR,
     $SINGULAR_THRESHOLD,0,INVERSE_SUCCESS)
      IF(.NOT.INVERSE_SUCCESS)THEN
         SUCCESS=.FALSE.
         RETURN
      ENDIF

C  Calculate the estimated parameters.
      DO IPARAM=1,NPARAMS
         C(IPARAM)=0.
         DO JPARAM=1,NPARAMS
            C(IPARAM)=C(IPARAM)+A_INV(IPARAM,JPARAM)*B(JPARAM)
         ENDDO
      ENDDO

C  Calculate the total sum of squares of the response variable data.
      SSTOT=0.
      DO IDATUM=1,NDATA
         SSTOT=SSTOT+Y(IDATUM)**2
      ENDDO

C  Calculate the sum of squares of the response variable data due to the
C  regression.
      SSREG=0.
      DO IPARAM=1,NPARAMS
         SSREG=SSREG+C(IPARAM)*B(IPARAM)
      ENDDO

C  Calculate the residual sum of squares of the response variable data.
      SSRES=SSTOT-SSREG
      IF(SSRES.LT.0.)THEN
         SSRES=0.
      ENDIF

C  Estimate the variance of the response variable data about the
C  regression.
      SSQ=SSRES/FLOAT(NDATA-NPARAMS)

C  Calculate the variances of the estimated parameters.
      DO IPARAM=1,NPARAMS
         VAR_C(IPARAM)=SSQ*A_INV(IPARAM,IPARAM)
      ENDDO

C  Done.
      SUCCESS=.TRUE.
      RETURN
      END
