      SUBROUTINE SOR(RELAX,THRESH,MAX_ITERS,A,A_DIM1,N,B,X,MAX_CHANGE,
     $ITER,SUCCESS)

C  Thomas Matejka NOAA/NSSL 15 September 1994

C  This subroutine solves a linear system of equations by successive
C  overrelaxation.

C  Input:

C  RELAX is a real variable that specifies the overrelaxation factor.
C  If RELAX is 1., the system is relaxed to a solution.  If RELAX is
C  greater than 1., the system is overrelaxed to a solution.

C  THRESH is a real variable that specifies the criterion for
C  convergence of the solution.  When changes in the solution are not
C  greater than THRESH, it is assumed that the solution has converged.

C  MAX_ITERS is an integer variable that specifies the maximum number of
C  iterations to perform.

C  A is a two-dimensional real array.  A(I,J) specifies the coefficient
C  of the Jth variable in the Ith equation.  A should have a strong main
C  diagonal to guarantee convergence of the solution.

C  A_DIM_1 is an integer variable that specifies the first dimension of
C  A in the calling program.

C  N is an integer variable that specifies the number of equations and
C  variables in the system of equations.

C  B is a one-dimensional real array.  B(I) specifies the
C  right-hand-side of the Ith equation.

C  Output:

C  X is a one-dimensional real array.  X(I) returns the solution of the
C  Ith variable only if SUCCESS is .TRUE..

C  MAX_CHANGE is a real variable that returns the maximum change in a
C  variable in the last iteration.

C  ITER is an integer variable that returns the number of iterations
C  performed.

C  SUCCESS is a logical variable.  SUCCESS returns .TRUE. if and only if
C  a solution was found within the specified convergence criteria.

      IMPLICIT NONE
      LOGICAL SUCCESS
      INTEGER I,N,ITER,J,MAX_ITERS,A_DIM1
      REAL SUM,X_OLD,MAX_CHANGE,RELAX,THRESH
      REAL B(N),X(N)
      REAL A(A_DIM1,N)

C  Check that there are no zeros on the main diagonal.
      DO I=1,N
         IF(A(I,I).EQ.0.)THEN
            SUCCESS=.FALSE.
            RETURN
         ENDIF
      ENDDO

C  Make an initial guess for the solution.
      DO I=1,N
         X(I)=0.
      ENDDO

C  Iterate to a solution.
      ITER=0
      DOWHILE(ITER.LT.MAX_ITERS)
         ITER=ITER+1
         MAX_CHANGE=0.
         DO I=1,N
            SUM=0.
            DO J=1,N
                  SUM=SUM+A(I,J)*X(J)      
            ENDDO
            X_OLD=X(I)
            X(I)=X_OLD+RELAX*(B(I)-SUM)/A(I,I)
            MAX_CHANGE=AMAX1(MAX_CHANGE,ABS(X(I)-X_OLD))
         ENDDO

C  Check whether the solution has converged.
         IF(MAX_CHANGE.LE.THRESH)THEN
            SUCCESS=.TRUE.
            RETURN
         ENDIF
      ENDDO

C  The solution did not converge within the maximum allowed number of
C  iterations.
      SUCCESS=.FALSE.
      RETURN
      END       
