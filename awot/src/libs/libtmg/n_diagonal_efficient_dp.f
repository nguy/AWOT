      SUBROUTINE N_DIAGONAL_EFFICIENT_DP(N_DIAGS,
     $AP_IN,AP_IN_DIM_1,N,B_IN,
     $SOLUTION_OUT,SUCCESS)

C  Thomas Matejka NOAA/NSSL 1 April 1998

C  This subroutine solves a set of simultaneous n-diagonal linear
C  equations.  The diagonal coefficients are stored efficiently.  Double
C  precision arithmetic is used.

C  Input:

C  N_DIAGS (integer) specifies the number of diagonals centered on the
C  principal diagonal with non-zero coefficients in the system of
C  equations.  It must be odd.

C  AP_IN (2d real array 1:AP_IN_DIM_1,-(N_DIAGS-1)/2:(N_DIAGS-1)/2).
C  AP_IN(I,J) specifies the coefficient for the (I+J)th variable in the
C  Ith equation.  AP_IN(1,0) must not be zero.

C  AP_IN_DIM_1 (integer) specifies the first dimension of AP_IN in the
C  calling program.

C  N (integer) specifies the number of equations, which is also the
C  number of variables.

C  B_IN (1d real array 1:N).  B_IN(I) specifies the right-hand side of
C  the Ith equation.

C  Output:

C  SOLUTION_OUT (1d real array 1:N).  SOLUTION_OUT(I) returns the
C  solution for the Ith variable.

C  SUCCESS (logical) returns .TRUE. if and only if a solution was found.

      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      LOGICAL::SUCCESS
      INTEGER::I,J,JJ,AP_IN_DIM_1,N,N_DIAGS,P,JMAX
      REAL,DIMENSION(N)::B_IN,SOLUTION_OUT
      REAL,DIMENSION(AP_IN_DIM_1,-(N_DIAGS-1)/2:(N_DIAGS-1)/2)::AP_IN
      DOUBLEPRECISION::S,SUM
      DOUBLEPRECISION,DIMENSION(N)::B,B_RR,SOLUTION
      DOUBLEPRECISION,DIMENSION(-(N_DIAGS-1)/2:(N_DIAGS-1)/2)::AP_HOLD
      DOUBLEPRECISION,DIMENSION(N,(N_DIAGS-1)/2)::AP_RR
      DOUBLEPRECISION,DIMENSION(N,-(N_DIAGS-1)/2:(N_DIAGS-1)/2)::AP

C  Check that the number of diagonals is odd.
      IF(MOD(N_DIAGS,2).EQ.0)THEN
         WRITE(TMMLIB_MESSAGE_UNIT,*)'N_DIAGONAL_EFFICIENT_DP:  THE ',
     $   'NUMBER OF DIAGONALS MUST BE ODD.'
         STOP
      ENDIF

C  Check that the first coefficient in the first equation is not zero.
      IF(AP_IN(1,0).EQ.0.)THEN
         SUCCESS=.FALSE.
         RETURN
      ENDIF

C  Calculate the number of non-zero diagonals on either side of the main
C  diagonal.
      P=(N_DIAGS-1)/2

C  Copy the coefficients and the right-hand side to double precision
C  arrays.
      DO I=1,N
         B(I)=DBLE(B_IN(I))
         DO J=-P,P
            AP(I,J)=DBLE(AP_IN(I,J))
         ENDDO
      ENDDO

C  Calculate the terms above the main diagonal and the right-hand side
C  of the augmented matrix row reduced to upper diagonal form.
      DO I=1,N
         DO J=-P,P
            AP_HOLD(J)=AP(I,J)
         ENDDO
         B_RR(I)=B(I)
         DO J=-P,-1
            S=AP_HOLD(J)
            AP_HOLD(J)=0.D0
            DO JJ=J+1,J+P
               AP_HOLD(JJ)=AP_HOLD(JJ)-S*AP_RR(I+J,JJ-J)
            ENDDO
            B_RR(I)=B_RR(I)-S*B_RR(I+J)
         ENDDO
         S=AP_HOLD(0)
         IF(S.EQ.0.D0)THEN
            SUCCESS=.FALSE.
            RETURN
         ENDIF
         AP_HOLD(0)=1.D0
         DO J=1,P
            AP_HOLD(J)=AP_HOLD(J)/S
         ENDDO
         B_RR(I)=B_RR(I)/S
         DO J=1,P
            AP_RR(I,J)=AP_HOLD(J)
         ENDDO
      ENDDO

C  Calculate the solution.
      DO I=N,1,-1
         SUM=0.D0
         JMAX=MIN0(P,N-I)
         IF(JMAX.GE.1)THEN
            DO J=1,JMAX
               SUM=SUM+AP_RR(I,J)*SOLUTION(I+J)
            ENDDO
         ENDIF
         SOLUTION(I)=B_RR(I)-SUM
      ENDDO

C  Copy the solution to a real array.
      DO I=1,N
         SOLUTION_OUT(I)=SNGL(SOLUTION(I))
      ENDDO

C  Done.
      SUCCESS=.TRUE.

      END SUBROUTINE N_DIAGONAL_EFFICIENT_DP
