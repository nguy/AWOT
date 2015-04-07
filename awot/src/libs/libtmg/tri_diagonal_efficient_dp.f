      SUBROUTINE TRI_DIAGONAL_EFFICIENT_DP(AP,AP_DIM_1,N,B,SOLUTION,
     $SUCCESS)

C  Thomas Matejka NOAA/NSSL 14 February 2002

C  This subroutine solves a set of simultaneous tri-diagonal linear
C  equations.  The diagonal coefficients are stored efficiently.  Double
C  precision arithmetic is used.

C  Input:

C  AP (2D doubleprecision array 1:AP_DIM_1,-1:1). AP(I,J) specifies the
C  coefficient for the (I+J)th variable in the Ith equation. AP(1,0)
C  must not be zero.

C  AP_DIM_1 (integer) specifies the first dimension of AP in the calling
C  program.

C  N (integer) specifies the number of equations, which is also the
C  number of variables.

C  B (1D doubleprecision array 1:N). B(I) specifies the right-hand side
C  of the Ith equation.

C  Output:

C  SOLUTION (1D doubleprecision array 1:N).  SOLUTION(I) returns the
C  solution for the Ith variable.

C  SUCCESS (logical) returns .TRUE. if and only if a solution was found.

      IMPLICIT NONE
      LOGICAL::SUCCESS
      INTEGER::I,AP_DIM_1,N
      DOUBLEPRECISION::TERM
      DOUBLEPRECISION,DIMENSION(N-1)::A_RR
      DOUBLEPRECISION,DIMENSION(N)::B,SOLUTION,B_RR
      DOUBLEPRECISION,DIMENSION(AP_DIM_1,-1:1)::AP

C  Initialize
      SUCCESS=.FALSE.

C  Check for form.
      IF(N.LT.1)THEN
         RETURN
      ENDIF
      IF(AP(1,-1).NE.0.D0)THEN
         RETURN
      ENDIF
      IF(AP(N,1).NE.0.D0)THEN
         RETURN
      ENDIF

C  Calculate the terms on the diagonal above the main diagonal and the
C  right hand side of the augmented matrix row reduced to upper diagonal
C  form.
      IF(AP(1,0).EQ.0.D0)THEN
         RETURN
      ENDIF
      A_RR(1)=AP(1,1)/AP(1,0)
      B_RR(1)=B(1)/AP(1,0)
      IF(N.GE.3)THEN
         DO I=2,N-1
            TERM=AP(I,0)-A_RR(I-1)*AP(I,-1)
            IF(TERM.EQ.0.D0)THEN
               RETURN
            ENDIF
            A_RR(I)=AP(I,1)/TERM
            B_RR(I)=(B(I)-B_RR(I-1)*AP(I,-1))/TERM
         ENDDO
      ENDIF
      IF(N.GE.2)THEN
         TERM=AP(N,0)-A_RR(N-1)*AP(N,-1)
         IF(TERM.EQ.0.D0)THEN
            RETURN
         ENDIF
         B_RR(N)=(B(N)-B_RR(N-1)*AP(N,-1))/TERM
      ENDIF

C  Calculate the solution.
      SOLUTION(N)=B_RR(N)
      IF(N.GE.2)THEN
         DO I=N-1,1,-1
            SOLUTION(I)=B_RR(I)-A_RR(I)*SOLUTION(I+1)
         ENDDO
      ENDIF
      SUCCESS=.TRUE.

      END SUBROUTINE TRI_DIAGONAL_EFFICIENT_DP
