      SUBROUTINE REGRESSION_1(NPARAMS,NDATA,X,X_DIM_1,Y,
     $BADDATA,
     $CHECK_SINGULAR,SINGULAR_THRESHOLD,
     $C)

C  Thomas Matejka NOAA/NSSL 24 July 2002

C  This subroutine performs a linear least squares regression involving
C  any number of independent variables. The response variable data are
C  all assumed to have the same variance and to be uncorrelated.

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

      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      LOGICAL::CHECK_SINGULAR,INVERSE_SUCCESS
      INTEGER::NDATA,NPARAMS,IDATUM,IPARAM,JPARAM,X_DIM_1
      REAL::BADDATA,SINGULAR_THRESHOLD
      REAL,DIMENSION(NPARAMS)::B,C
      REAL,DIMENSION(NDATA)::Y
      REAL,DIMENSION(X_DIM_1,NPARAMS)::X
      REAL,DIMENSION(NPARAMS,NPARAMS)::A,A_INV

C  There must be at least as many data as parameters.
      IF(NDATA.LT.NPARAMS)THEN
         DO IPARAM=1,NPARAMS
            C(IPARAM)=BADDATA
         ENDDO
         RETURN
      ENDIF

C  Calculate the matrix of the sum of products of the predictor variable
C  data.
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

C  Find the inverse of the matrix of the sum of products of the
C  predictor variable data.
      CALL MAT_INV(A,NPARAMS,NPARAMS,A_INV,NPARAMS,CHECK_SINGULAR,
     $SINGULAR_THRESHOLD,0,INVERSE_SUCCESS)
      IF(.NOT.INVERSE_SUCCESS)THEN
         DO IPARAM=1,NPARAMS
            C(IPARAM)=BADDATA
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

      END SUBROUTINE REGRESSION_1
