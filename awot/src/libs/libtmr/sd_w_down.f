      SUBROUTINE SD_W_DOWN(Z1,TV1,Z2,TV2,Z3,RHO3,
     $N,Z,SD_DIV,
     $I_BC,W_BC,SD_W_BC,
     $N_TRIALS,
     $BADDATA,
     $SD_W,SUCCESS)

C  Thomas Matejka NOAA/NSSL 29 November 2001

C  This subroutine calculates the rms error of vertical air velocity in
C  a column when W_DOWN is used.

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

C  SD_DIV is a one-dimensional real array with elements indexed from 1
C  to N.  SD_DIV(I) specifies the rms error of DIV(I) (s**-1).  If it is
C  missing, it should equal BADDATA.

C  I_BC is an integer variable that specifies that the boundary is at
C  Z(I_BC).

C  W_BC is a real variable that specifies the vertical air velocity
C  boundary condition (m s**-1) at Z(I_BC).  If it is missing, it should
C  equal BADDATA.

C  SD_W_BC is a real variable that specifies the rms error of W_BC (m
C  s**-1).  If it is missing, it should equal BADDATA.

C  N_TRIALS is an integer variable that specifies the number of Monte
C  Carlo trials to use for estimating rms errors.

C  BADDATA is a real variable that indicates a missing value as
C  described.

C  Output:

C  SD_W is a one-dimensional real array with elements indexed from 1 to
C  N.  SD_W(I) returns the rms error of vertical air velocity (m s**-1)
C  at Z(I).  If it is missing, it returns BADDATA.

C  SUCCESS is a logical variable that returns .TRUE. if and only if a
C  vertical air velocity could be calculated at at least one level.

      IMPLICIT NONE
      LOGICAL::SUCCESS
      INTEGER::I,N,I_BC,N_TRIALS,I_TRIAL
      INTEGER,DIMENSION(N)::COUNT_W
      REAL::Z1,Z2,Z3,TV1,TV2,RHO3,BADDATA,W_BC,SD_W_BC
      REAL,DIMENSION(N)::Z,SD_DIV,SD_W,DIV,W,SUM_SQ_ERR_W

C  Initialize the counters and total squared errors to zero.
      DO I=1,N
         COUNT_W(I)=0
         SUM_SQ_ERR_W(I)=0.
      ENDDO

C  Loop through the trials.
      DO I_TRIAL=1,N_TRIALS

C  Generate errors for horizontal divergence.
         DO I=1,N
            IF(SD_DIV(I).NE.BADDATA)THEN
               CALL RNGN(88,0.,SD_DIV(I),DIV(I))
            ELSE
               DIV(I)=BADDATA
            ENDIF
         ENDDO

C  Generate errors for the vertical air velocity boundary condition.
         IF(SD_W_BC.NE.BADDATA)THEN
            CALL RNGN(88,0.,SD_W_BC,W_BC)
         ELSE
            W_BC=BADDATA
         ENDIF

C  Perform a downward integration.
         CALL W_DOWN(Z1,TV1,Z2,TV2,Z3,RHO3,
     $   N,Z,DIV,
     $   I_BC,W_BC,
     $   BADDATA,
     $   W,
     $   SUCCESS)
         IF(.NOT.SUCCESS)THEN
            EXIT
         ENDIF

C  Accumulate the squared error.
         DO I=1,N
            IF(W(I).NE.BADDATA)THEN
               COUNT_W(I)=COUNT_W(I)+1
               SUM_SQ_ERR_W(I)=SUM_SQ_ERR_W(I)+W(I)**2
            ENDIF
         ENDDO

C  Return for the next trial.
      ENDDO

C  Calculate the rms errors.
      DO I=1,N
         IF(COUNT_W(I).NE.0)THEN
            SD_W(I)=SQRT(SUM_SQ_ERR_W(I)/FLOAT(COUNT_W(I)))
         ELSE
            SD_W(I)=BADDATA
         ENDIF
      ENDDO

      END SUBROUTINE SD_W_DOWN
