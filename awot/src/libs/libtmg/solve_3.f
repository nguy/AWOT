      SUBROUTINE SOLVE_3(Y0,X1,X2,X_TOL,Y_TOL,MAX_ITERS,X0,SUCCESS)

C  Thomas Matejka NOAA/NSSL 6 May 1998
C=======================================================================
C  This subroutine solves F(X0)=Y0 for X0 by using the average of the
C  method of false position and the method of bisection.  The function
C  to be solved should be defined in a subroutine named SOLVE_FUNCTION.
C=======================================================================
      IMPLICIT NONE
      REAL,EXTERNAL::SOLVE_FUNCTION
      LOGICAL::SUCCESS
      INTEGER::MAX_ITERS,I_ITER
      REAL::X0,X1,X2,X_GUESS_1,X_GUESS_2,X_GUESS,X_MINUS,X_PLUS,
     $Y0,Y1,Y2,Y_GUESS,Y_MINUS,Y_PLUS,
     $X_TOL,Y_TOL

C  Find the function values at the endpoints of the solution domain.
      Y1=SOLVE_FUNCTION(X1)
      Y2=SOLVE_FUNCTION(X2)

C  Check if the solution is an endpoint of the solution domain.
      IF(Y0.EQ.Y1)THEN
         X0=X1
         SUCCESS=.TRUE.
         RETURN
      ENDIF
      IF(Y0.EQ.Y2)THEN
         X0=X2
         SUCCESS=.TRUE.
         RETURN
      ENDIF

C  Check that the value of the function at the endpoints of the solution
C  domain encloses the value of the function at the solution.
      IF(SIGN(1.,Y1-Y0).EQ.SIGN(1.,Y2-Y0))THEN
         SUCCESS=.FALSE.
         RETURN
      ENDIF

C  Initialize the solution domain.
      IF(Y2.GT.Y1)THEN
         X_PLUS=X2
         X_MINUS=X1
         Y_PLUS=Y2
         Y_MINUS=Y1
      ELSE
         X_PLUS=X1
         X_MINUS=X2
         Y_PLUS=Y1
         Y_MINUS=Y2
      ENDIF

C  Iterate to a solution.
      I_ITER=0
      DO
         I_ITER=I_ITER+1

C  Check whether the maximum number of iterations has been exceeded.
         IF(I_ITER.GT.MAX_ITERS)THEN
            SUCCESS=.FALSE.
            RETURN
         ENDIF

C  Estimate a solution.
         X_GUESS_1=X_MINUS+
     $   (Y0-Y_MINUS)*(X_PLUS-X_MINUS)/(Y_PLUS-Y_MINUS)
         X_GUESS_2=(X_MINUS+X_PLUS)/2.
         X_GUESS=(X_GUESS_1+X_GUESS_2)/2.
         Y_GUESS=SOLVE_FUNCTION(X_GUESS)

C  Check whether the solution has converged.
         IF(ABS(Y_GUESS-Y0).LE.Y_TOL.OR.
     $   ABS(X_PLUS-X_MINUS).LE.X_TOL)THEN
            X0=X_GUESS
            SUCCESS=.TRUE.
            RETURN
         ENDIF

C  Narrow the solution domain.
         IF(Y_GUESS.LT.Y0)THEN
            X_MINUS=X_GUESS
            Y_MINUS=Y_GUESS
         ELSE
            X_PLUS=X_GUESS
            Y_PLUS=Y_GUESS
         ENDIF
      ENDDO

      END SUBROUTINE SOLVE_3
