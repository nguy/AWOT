      SUBROUTINE SOR_MULTI_DIAGONAL_EFFICIENT(RELAX,THRESH,MAX_ITERS,
     $N_DIAGS,AP,AP_DIM1,N,B,X,MAX_CHANGE,ITER,SUCCESS)

C  Thomas Matejka NOAA/NSSL 19 September 1994

C  This subroutine solves a multi-diagonal linear system of equations by
C  successive overrelaxation.  The diagonal coefficients are stored
C  efficiently.

C  Input:

C  RELAX is a real variable that specifies the overrelaxation factor.
C  If RELAX is 1., the system is relaxed to a solution.  If RELAX is
C  greater than 1., the system is overrelaxed to a solution.

C  THRESH is a real variable that specifies the criterion for
C  convergence of the solution.  When changes in the solution are not
C  greater than THRESH, it is assumed that the solution has converged.

C  MAX_ITERS is an integer variable that specifies the maximum number of
C  iterations to perform.

C  N_DIAGS is an integer variable that specifies the number of non-zero
C  diagonals in the linear system.  N_DIAGS must be odd.

C  AP is a two-dimensional real array.  The second dimension of AP is
C  indexed from -N_DIAGS/2 to N_DIAGS/2.  AP(I,J) specifies the
C  coefficient of the (I+J)th variable in the Ith equation.  AP should
C  have a strong main diagonal to guarantee convergence of the solution.
C  AP(I,J) should be 0. when I+J is less than 1 or greater than N.

C  AP_DIM_1 is an integer variable that specifies the first dimension of
C  AP in the calling program.

C  N is an integer variable that specifies the number of equations and
C  variables in the system of equations.

C  B is a one-dimensional real array.  B(I) specifies the
C  right-hand-side of the Ith equation.

C  Output:

C  X is a one-dimensional real array.  X(I) returns the solution for the
C  Ith variable only if SUCCESS is .TRUE..

C  MAX_CHANGE is a real variable that returns the maximum change in a
C  variable in the last iteration.

C  ITER is an integer variable that returns the number of iterations
C  performed.

C  SUCCESS is a logical variable.  SUCCESS returns .TRUE. if and only if
C  a solution was found within the specified convergence criteria.

      IMPLICIT NONE
      LOGICAL SUCCESS
      INTEGER I,N,ITER,MAX_ITERS,AP_DIM1,N_DIAGS,J
      REAL SUM,X_OLD,MAX_CHANGE,RELAX,THRESH
      REAL B(N),X(N)
      REAL AP(AP_DIM1,-N_DIAGS/2:N_DIAGS/2)

C  Check that there are no zeros on the main diagonal.
      DO I=1,N
         IF(AP(I,0).EQ.0.)THEN
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
            DO J=-N_DIAGS/2,N_DIAGS/2
               SUM=SUM+AP(I,J)*X(I+J)
            ENDDO
            X_OLD=X(I)
            X(I)=X_OLD+RELAX*(B(I)-SUM)/AP(I,0)
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
