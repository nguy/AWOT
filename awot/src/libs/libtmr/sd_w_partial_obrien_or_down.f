      SUBROUTINE SD_W_PARTIAL_OBRIEN_OR_DOWN(Z1,TV1,Z2,TV2,Z3,RHO3,
     $N,Z,SD_DIV,
     $I_BCLO,W_BCLO,SD_W_BCLO,I_BCHI,W_BCHI,SD_W_BCHI,DELZ_MAX,
     $N_TRIALS,
     $BADDATA,
     $SD_DIV_OUT,SD_W,SUCCESS)

C  Thomas Matejka NOAA/NSSL 29 November 2001

C  This subroutine calculates the rms errors of adjusted horizontal
C  divergence and vertical air velocity in a column when
C  W_PARTIAL_OBRIEN_OR_DOWN is used.

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

C  I_BCLO is an integer variable that specifies that the lower boundary
C  is at Z(I_BCLO).

C  W_BCLO is a real variable that specifies the vertical air velocity
C  boundary condition (m s**-1) at Z(I_BCLO).  If it is missing, it
C  should equal BADDATA.

C  SD_W_BCLO is a real variable that specifies the rms error of W_BCLO
C  (m s**-1).  If it is missing, it should equal BADDATA.

C  I_BCHI is an integer variable that specifies that the upper boundary
C  is at Z(I_BCHI).

C  W_BCHI is a real variable that specifies the vertical air velocity
C  boundary condition (m s**-1) at Z(I_BCHI).  If it is missing, it
C  should equal BADDATA.

C  SD_W_BCHI is a real variable that specifies the rms error of W_BCHI
C  (m s**-1).  If it is missing, it should equal BADDATA.

C  DELZ_MAX is a real variable that specifies the depth (m) of the layer
C  of missing divergences at and above the lower boundary that must not
C  be equaled or exceeded for the partial O'Brien solution to be used.

C  N_TRIALS is an integer variable that specifies the number of Monte
C  Carlo trials to use for estimating rms errors.

C  BADDATA is a real variable that indicates a missing value as
C  described.

C  Output:

C  SD_DIV_OUT is a one-dimensional real array with elements indexed from
C  1 to N.  SD_DIV_OUT(I) returns the rms error of possibly adjusted
C  horizontal divergence (s**-1) at Z(I).  If it is missing, it returns
C  BADDATA.  SD_DIV_OUT may be the same as SD_DIV, in which case SD_DIV
C  is overwritten.

C  SD_W is a one-dimensional real array with elements indexed from 1 to
C  N.  SD_W(I) returns the rms error of vertical air velocity (m s**-1)
C  at Z(I).  If it is missing, it returns BADDATA.

C  SUCCESS is a logical variable that returns .TRUE. if and only if a
C  vertical air velocity could be calculated at at least one level.

      IMPLICIT NONE
      LOGICAL::LDUM,SUCCESS
      INTEGER::I,N,I_BCLO,I_BCHI,N_TRIALS,I_TRIAL
      INTEGER,DIMENSION(N)::COUNT_DIV,COUNT_W
      REAL::Z1,Z2,Z3,TV1,TV2,RHO3,DELZ_MAX,BADDATA,W_BCLO,W_BCHI,
     $SD_W_BCLO,SD_W_BCHI
      REAL,DIMENSION(N)::Z,SD_DIV,SD_DIV_OUT,SD_W,SD_DIV_TEMP,DIV,W,
     $SUM_SQ_ERR_DIV,SUM_SQ_ERR_W

C  Copy the input array in case it is overwritten.
      DO I=1,N
         SD_DIV_TEMP(I)=SD_DIV(I)
      ENDDO

C  Initialize the counters and total squared errors to zero.
      DO I=1,N
         COUNT_DIV(I)=0
         COUNT_W(I)=0
         SUM_SQ_ERR_DIV(I)=0.
         SUM_SQ_ERR_W(I)=0.
      ENDDO

C  Loop through the trials.
      DO I_TRIAL=1,N_TRIALS

C  Generate errors for horizontal divergence.
         DO I=1,N
            IF(SD_DIV_TEMP(I).NE.BADDATA)THEN
               CALL RNGN(88,0.,SD_DIV_TEMP(I),DIV(I))
            ELSE
               DIV(I)=BADDATA
            ENDIF
         ENDDO

C  Generate errors for the vertical air velocity boundary conditions.
         IF(SD_W_BCLO.NE.BADDATA)THEN
            CALL RNGN(88,0.,SD_W_BCLO,W_BCLO)
         ELSE
            W_BCLO=BADDATA
         ENDIF
         IF(SD_W_BCHI.NE.BADDATA)THEN
            CALL RNGN(88,0.,SD_W_BCHI,W_BCHI)
         ELSE
            W_BCHI=BADDATA
         ENDIF

C  Perform a partial O'Brien solution or downward integration.
         CALL W_PARTIAL_OBRIEN_OR_DOWN(Z1,TV1,Z2,TV2,Z3,RHO3,
     $   N,Z,DIV,SD_DIV_TEMP,
     $   I_BCLO,W_BCLO,I_BCHI,W_BCHI,DELZ_MAX,
     $   BADDATA,
     $   DIV,W,
     $   LDUM,LDUM,
     $   SUCCESS)
         IF(.NOT.SUCCESS)THEN
            EXIT
         ENDIF

C  Accumulate the squared error.
         DO I=1,N
            IF(DIV(I).NE.BADDATA)THEN
               COUNT_DIV(I)=COUNT_DIV(I)+1
               SUM_SQ_ERR_DIV(I)=SUM_SQ_ERR_DIV(I)+DIV(I)**2
            ENDIF
            IF(W(I).NE.BADDATA)THEN
               COUNT_W(I)=COUNT_W(I)+1
               SUM_SQ_ERR_W(I)=SUM_SQ_ERR_W(I)+W(I)**2
            ENDIF
         ENDDO

C  Return for the next trial.
      ENDDO

C  Calculate the rms errors.
      DO I=1,N
         IF(COUNT_DIV(I).NE.0)THEN
            SD_DIV_OUT(I)=SQRT(SUM_SQ_ERR_DIV(I)/FLOAT(COUNT_DIV(I)))
         ELSE
            SD_DIV_OUT(I)=BADDATA
         ENDIF
         IF(COUNT_W(I).NE.0)THEN
            SD_W(I)=SQRT(SUM_SQ_ERR_W(I)/FLOAT(COUNT_W(I)))
         ELSE
            SD_W(I)=BADDATA
         ENDIF
      ENDDO

      END SUBROUTINE SD_W_PARTIAL_OBRIEN_OR_DOWN
