      SUBROUTINE LLS_VAR_RESPONSE(NPARAMS,AINV,AINV_DIM_1,X0,SSQ,VAR_Y0,
     $SUCCESS)

C  Thomas Matejka NOAA/NSSL 12 January 1994

C  This subroutine determines the variance of the response variable
C  calculated from a linear least-squares regression equation involving
C  any number of parameters.

C  Input:

C  NPARAMS is an integer variable that specifies the number of predictor
C  variables (including the constant term if used) in the regression
C  equation.

C  AINV is a two-dimensional real array.  AINV specifies the inverse of
C  A.  A(I,J) is the covariance of the possibly transformed Ith and Jth
C  predictor variables from the development of the regression.

C  AINV_DIM_1 is an integer variable that specifies the first dimension
C  of AINV in the calling program.

C  X0 is a one-dimensional real array.  X0(I) specifies the value of the
C  Ith predictor variable for the point of interest.

C  SSQ is a real variable that specifies the variance of the possibly
C  transformed response variable data about the regression.

C  Output:

C  VAR_Y0 is a real variable that returns the variance of the response
C  variable calculated from the regression at the point of interest.

C  SUCCESS is a logical variable that returns .TRUE. if and only if the
C  variance of the response variable calculated at the point of interest
C  can be calculated.  When the regression was developed with correlated
C  predictor variable data, VAR_Y0 may occasionally be negative.  This
C  is a spurious result.

      IMPLICIT NONE
      LOGICAL SUCCESS
      INTEGER NPARAMS,AINV_DIM_1,I,J
      REAL SSQ,VAR_Y0,SUM
      REAL X0(NPARAMS)
      REAL AINV(AINV_DIM_1,NPARAMS)

C  Calculate the variance of the response variable at the specified
C  point.
      SUM=0.
      DO 1 J=1,NPARAMS
         DO 2 I=1,NPARAMS
            SUM=SUM+X0(I)*X0(J)*AINV(I,J)
2        CONTINUE
1     CONTINUE
      VAR_Y0=SSQ*SUM
      IF(VAR_Y0.LT.0)THEN
         SUCCESS=.FALSE.
         RETURN
      ENDIF

C  Done.
      SUCCESS=.TRUE.
      RETURN
      END
