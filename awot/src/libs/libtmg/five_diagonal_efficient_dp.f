      SUBROUTINE FIVE_DIAGONAL_EFFICIENT_DP(AP,AP_DIM_1,N,B,SOLUTION,
     $SUCCESS)

C  Thomas Matejka NOAA/NSSL 21 September 1994

C  This subroutine solves a set of simultaneous five-diagonal linear
C  equations.  The diagonal coefficients are stored efficiently.  Double
C  precision arithmetic is used.

C  Input:

C  AP is a two-dimensional double precision array.  The second dimension
C  of AP is indexed from -2 to 2.  AP(I,J) specifies the coefficient for
C  the (I+J)th variable in the Ith equation.  AP(1,0) must not be zero.

C  AP_DIM_1 is an integer variable that specifies the first dimension of
C  AP in the calling program.

C  N is an integer variable that specifies the number of equations,
C  which is also the number of variables.

C  B is a one-dimensional double precision array.  B(I) specifies the
C  right-hand side of the Ith equation.

C  Output:

C  SOLUTION is a one-dimensional double precision array.  SOLUTION(I)
C  returns the solution for the Ith variable.

C  SUCCESS is a logical variable that returns .TRUE. if and only if a
C  solution was found.

      IMPLICIT NONE
      LOGICAL SUCCESS
      INTEGER I,AP_DIM_1,N
      DOUBLEPRECISION TERM1,TERM2
      DOUBLEPRECISION A2_RR(N-2)
      DOUBLEPRECISION A1_RR(N-1)
      DOUBLEPRECISION B(N),SOLUTION(N),B_RR(N)
      DOUBLEPRECISION AP(AP_DIM_1,-2:2)

C  Calculate the terms on the diagonal above the main diagonal and the
C  right hand side of the augmented matrix row reduced to upper diagonal
C  form.
      IF(AP(1,0).EQ.0.D0)THEN
         SUCCESS=.FALSE.
         RETURN
      ENDIF
      A1_RR(1)=AP(1,1)/AP(1,0)
      A2_RR(1)=AP(1,2)/AP(1,0)
      B_RR(1)=B(1)/AP(1,0)
      TERM1=AP(2,0)-AP(2,-1)*A1_RR(1)
      A1_RR(2)=(AP(2,1)-AP(2,-1)*A2_RR(1))/TERM1
      A2_RR(2)=AP(2,2)/TERM1
      B_RR(2)=(B(2)-AP(2,-1)*B_RR(1))/TERM1
      DO I=3,N-2
         TERM1=AP(I,-1)-AP(I,-2)*A1_RR(I-2)
         TERM2=AP(I,0)-AP(I,-2)*A2_RR(I-2)-TERM1*A1_RR(I-1)
         IF(TERM2.EQ.0.D0)THEN
            SUCCESS=.FALSE.
            RETURN
         ENDIF
         A1_RR(I)=(AP(I,1)-TERM1*A2_RR(I-1))/TERM2
         A2_RR(I)=AP(I,2)/TERM2
         B_RR(I)=(B(I)-AP(I,-2)*B_RR(I-2)-TERM1*B_RR(I-1))/TERM2
      ENDDO
      TERM1=AP(N-1,-1)-AP(N-1,-2)*A1_RR(N-3)
      TERM2=AP(N-1,0)-AP(N-1,-2)*A2_RR(N-3)-TERM1*A1_RR(N-2)
      IF(TERM2.EQ.0.D0)THEN
         SUCCESS=.FALSE.
         RETURN
      ENDIF
      A1_RR(N-1)=(AP(N-1,1)-TERM1*A2_RR(N-2))/TERM2
      B_RR(N-1)=(B(N-1)-AP(N-1,-2)*B_RR(N-3)-TERM1*B_RR(N-2))/TERM2
      TERM1=AP(N,-1)-AP(N,-2)*A1_RR(N-2)
      TERM2=AP(N,0)-AP(N,-2)*A2_RR(N-2)-TERM1*A1_RR(N-1)
      IF(TERM2.EQ.0.D0)THEN
         SUCCESS=.FALSE.
         RETURN
      ENDIF
      B_RR(N)=(B(N)-AP(N,-2)*B_RR(N-2)-TERM1*B_RR(N-1))/TERM2

C  Calculate the solution.
      SOLUTION(N)=B_RR(N)
      SOLUTION(N-1)=B_RR(N-1)-A1_RR(N-1)*SOLUTION(I+1)
      DO I=N-2,1,-1
         SOLUTION(I)=B_RR(I)-A1_RR(I)*SOLUTION(I+1)-
     $   A2_RR(I)*SOLUTION(I+2)
      ENDDO

C  Done.
      SUCCESS=.TRUE.
      RETURN
      END
