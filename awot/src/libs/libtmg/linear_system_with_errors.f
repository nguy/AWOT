      SUBROUTINE LINEAR_SYSTEM_WITH_ERRORS(N,A,A_DIM_1,B,RMSE_B,
     $CHECK_SINGULAR,SINGULAR_THRESHOLD,
     $X,RMSE_X,SUCCESS)

C  Thomas Matejka NOAA/NSSL 1 September 1998

C  This subroutine solves a linear system and calculates the rms errors
C  of the solution.

C  Input:

C  N (integer) specifies the number of equations and the number of
C  unknown variables.

C  A (2d real array 1:A_DIM_1,N).  A(I,J) specifies the coefficient of
C  the Jth variable in the Ith equation.

C  A_DIM_1 (integer) specifies the first dimension of A in the calling
C  program.

C  B (1d real array 1:N).  B(I) specifies the right hand side of the Ith
C  equation.

C  RMSE_B (1d real array 1:N).  RMSE_B(I) specifies the rms error
C  variance of B(I).

C  CHECK_SINGULAR (logical) controls whether the matrix to be inverted
C  is checked for near-singularity.  If it is .TRUE., the matrix is
C  checked.  If it is .FALSE., it is not checked.

C  SINGULAR_THRESHOLD (real) controls the definition of nearly singular.
C  The smaller SINGULAR_THRESHOLD is, the more close to exactly singular
C  matrices must be to be considered nearly singular.  It must be
C  greater than or equal to 0..  It is relevant only when CHECK_SINGULAR
C  is .TRUE..

C  Output:

C  X (1d real array 1:N).  X(J) returns the solution for the Jth
C  variable.  It is meaningful only if SUCCESS is .TRUE..

C  RMSE_X (1d real array 1:N).  RMSE_X(J) returns the rms error of the
C  solution for the Jth variable.  It is meaningful only if SUCCESS is
C  .TRUE..

C  SUCCESS (logical) indicates whether a solution was found.  It returns
C  .TRUE. if a solution was found.  It returns .FALSE. if it was not
C  found.

      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      LOGICAL::SUCCESS,CHECK_SINGULAR,INVERSE_SUCCESS
      INTEGER::N,A_DIM_1,I,J,II,JJ
      REAL::SINGULAR_THRESHOLD
      REAL,DIMENSION(N)::B,RMSE_B,BB,Q,X,RMSE_X
      REAL,DIMENSION(A_DIM_1,N)::A
      REAL,DIMENSION(N,N)::AA,P,P_INV

C  Scale the equations by the rms error of the right hand side.
      DO I=1,N
         DO J=1,N
            AA(I,J)=A(I,J)/RMSE_B(I)
         ENDDO
         BB(I)=B(I)/RMSE_B(I)
      ENDDO

C  Create correlation matrices.
      DO II=1,N
         DO JJ=1,II
            P(II,JJ)=0.
            DO I=1,N
               P(II,JJ)=P(II,JJ)+AA(I,II)*AA(I,JJ)
            ENDDO
         ENDDO
         Q(II)=0.
         DO I=1,N
            Q(II)=Q(II)+AA(I,II)*BB(I)
         ENDDO
      ENDDO
      IF(N.GE.2)THEN
         DO II=1,N-1
            DO JJ=II+1,N
               P(II,JJ)=P(JJ,II)
            ENDDO
         ENDDO
      ENDIF

C  Find the inverse of the correlation matrix.
      CALL MAT_INV(P,N,N,P_INV,N,CHECK_SINGULAR,SINGULAR_THRESHOLD,0,
     $INVERSE_SUCCESS)
      IF(.NOT.INVERSE_SUCCESS)THEN
         SUCCESS=.FALSE.
         RETURN
      ENDIF

C  Calculate the solution.
      DO II=1,N
         X(II)=0.
         DO JJ=1,N
            X(II)=X(II)+P_INV(II,JJ)*BB(JJ)
         ENDDO
      ENDDO

C  Calculate the rms errors of the solution.
      DO II=1,N
         RMSE_X(II)=P_INV(II,II)
      ENDDO

C  Done.
      SUCCESS=.TRUE.

      END SUBROUTINE LINEAR_SYSTEM_WITH_ERRORS
