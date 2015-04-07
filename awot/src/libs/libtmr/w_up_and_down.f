      SUBROUTINE W_UP_AND_DOWN(Z1,TV1,Z2,TV2,Z3,RHO3,
     $N,Z,DIV,
     $I_BC,W_BC,
     $BADDATA,
     $W,
     $SUCCESS)

C  Thomas Matejka NOAA/NSSL 18 April 1997

C  This subroutine calculates vertical air velocity in a column using
C  W_UP and W_DOWN.

C  Input:

C  Z1 is a real variable that specifies an altitude (m).

C  TV1 is a real variable that specifies the virtual temperature (K) at
C  Z1.

C  Z2 is a real variable that specifies an altitude (m).  Z2 must not be
C  equal to Z1.

C  TV2 is a real variable that specifies the virtual temperature (K) at
C  Z2.

C  Z3 is a real variable that specifies an altitude (m).  Z3 may be
C  equal to Z1 or Z2.

C  RHO3 is a real variable that specifies the density (kg m**-3) at Z3.
C  Z1, TV1, Z2, TV2, Z3, and RHO3 define a density profile in a
C  constant-lapse-rate atmosphere.

C  N is an integer variable that specifies the number of data levels in
C  the column.

C  Z is a one-dimensional real array with elements indexed from 1 to N.
C  Z(I) specifies the altitude (m) of the Ith level, counting from the
C  bottom of the column.

C  DIV is a one-dimensional real array with elements indexed from 1 to
C  N.  DIV(I) specifies the horizontal divergence (s**-1) at Z(I).  If
C  it is missing, it should equal BADDATA.

C  I_BC is an integer variable that specifies that the boundary is at
C  Z(I_BC).

C  W_BC is a real variable that specifies the vertical air velocity
C  boundary condition (m s**-1) at Z(I_BC).  If it is missing, it should
C  equal BADDATA.

C  BADDATA is a real variable that indicates a missing value as
C  described.

C  Output:

C  W is a one-dimensional real array with elements indexed from 1 to N.
C  W(I) returns the vertical air velocity (m s**-1) at Z(I).  If it is
C  missing, it returns BADDATA.

C  SUCCESS is a logical variable that returns .TRUE. if and only if at
C  least one vertical air velocity above or below the boundary was
C  calculated.

      IMPLICIT NONE
      LOGICAL::SUCCESS,SUCCESS_DUM
      INTEGER::I,N,I_BC
      REAL::Z1,Z2,Z3,TV1,TV2,RHO3,W_BC,BADDATA
      REAL,DIMENSION(N)::DIV,Z,W,W_DUM

C  Initialize the vertical air velocity to missing.
      DO I=1,N
         W(I)=BADDATA
      ENDDO

C  Initialize the success to false.
      SUCCESS=.FALSE.

C  Upward integration.
      CALL W_UP(Z1,TV1,Z2,TV2,Z3,RHO3,
     $N,Z,DIV,
     $I_BC,W_BC,
     $BADDATA,
     $W_DUM,
     $SUCCESS_DUM)
      IF(SUCCESS_DUM)THEN
         SUCCESS=.TRUE.
         DO I=I_BC,N
            W(I)=W_DUM(I)
         ENDDO
      ENDIF
      
C  Downward integration.
      CALL W_DOWN(Z1,TV1,Z2,TV2,Z3,RHO3,
     $N,Z,DIV,
     $I_BC,W_BC,
     $BADDATA,
     $W_DUM,
     $SUCCESS_DUM)
      IF(SUCCESS_DUM)THEN
         SUCCESS=.TRUE.
         DO I=1,I_BC
            W(I)=W_DUM(I)
         ENDDO
      ENDIF

      END SUBROUTINE W_UP_AND_DOWN
