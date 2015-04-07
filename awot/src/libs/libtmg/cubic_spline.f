      SUBROUTINE CUBIC_SPLINE(N,X,Q,GPP,SUCCESS)

C  Thomas Matejka NOAA/NSSL 29 April 1998
C=======================================================================
C  This subroutine calculates the second derivatives of a
C  one-dimensional natural cubic spline at the knots.  (For a natural
C  cubic spline, the second derivatives at the first and last knots are
C  assumed to be zero.)  These second derivatives and the values of the
C  cubic spline at the knots provide four coefficients per segment,
C  which is sufficient to define the cubic spline.

C  Input:

C  N (integer) specifies the number of knots.  N must be at least 2.

C  X (1d real array 1:N).  X(I) specifies the coordinate of the Ith
C  knot.  X(I) must be monotonically increasing for I=1,...,N.

C  Q (1d real array 1:N).  Q(I) specifies the value of the cubic spline
C  at X(I).

C  Output:

C  GPP (1d real array 1:N).  GPP(I) returns the value of the second
C  derivative of the cubic spline at X(I).

C  SUCCESS (logical) returns .TRUE. if and only if a solution was found.
C=======================================================================
      IMPLICIT NONE
      LOGICAL::SUCCESS
      INTEGER::N,I,J
      REAL,DIMENSION(1:N)::X
      REAL,DIMENSION(1:N)::Q
      REAL,DIMENSION(1:N)::GPP
      REAL,DIMENSION(1:N)::B
      REAL,DIMENSION(1:N,-1:1)::A

C  Check that there are at least two knots.
      IF(N-1.LT.1)THEN
         SUCCESS=.FALSE.
         RETURN
      ENDIF

C  Calculate the second derivative of the cubic spline at the knots.
      IF(N.GE.3)THEN
         DO I=1,N
            B(I)=0.
            DO J=-1,1
               A(I,J)=0.
            ENDDO
         ENDDO
         A(1,0)=1.
         DO I=2,N-1
            A(I,-1)=X(I)-X(I-1)
            A(I,0)=2.*(X(I+1)-X(I-1))
            A(I,1)=X(I+1)-X(I)
            B(I)=6.*
     $      ((Q(I+1)-Q(I))/(X(I+1)-X(I))-(Q(I)-Q(I-1))/(X(I)-X(I-1)))
         ENDDO
         A(N,0)=1.
         CALL N_DIAGONAL_EFFICIENT_DP(3,
     $   A,N,N,B,
     $   GPP,SUCCESS)
         IF(.NOT.SUCCESS)THEN
            RETURN
         ENDIF
      ELSE
         GPP(1)=0.
         GPP(N)=0.
      ENDIF

C  Done.
      SUCCESS=.TRUE.

      END SUBROUTINE CUBIC_SPLINE
