      SUBROUTINE FIELD_STATISTICS_3D(A,MAXX,MAXY,MAXZ,IX_START,IX_END,
     $IY_START,IY_END,IZ_START,IZ_END,BADDATA,N_GOOD,N_BAD,MEAN,SD,
     $A_MIN,IX_A_MIN,IY_A_MIN,IZ_A_MIN,A_MAX,IX_A_MAX,IY_A_MAX,IZ_A_MAX,
     $IX_MIN,IX_MAX,IY_MIN,IY_MAX,IZ_MIN,IZ_MAX)

C  Thomas Matejka NOAA/NSSL 18 November 1996

C  This subroutine calculates statistics for a rectangular subset of a
C  three-dimensional field.

C  Input:

C  A is a three-dimensional real array.  A(I,J,K) specifies the value,
C  at the Ith grid point in the first dimension, the Jth grid point in
C  the second dimension, and the Kth grid point in the third dimension,
C  of the data field whose statistics are sought.  If it is missing, it
C  should equal BADDATA.

C  MAXX is an integer variable that specifies the first dimension of A
C  in the calling program.

C  MAXY is an integer variable that specifies the second dimension of A
C  in the calling program.

C  MAXZ is an integer variable that specifies the third dimension of A
C  in the calling program.

C  IX_START is an integer variable that specifies the first grid point
C  number in the first dimension of the rectangular subset of A for
C  which statistics will be compiled.

C  IX_END is an integer variable that specifies the last grid point
C  number in the first dimension of the rectangular subset of A for
C  which statistics will be compiled.

C  IY_START is an integer variable that specifies the first grid point
C  number in the second dimension of the rectangular subset of A for
C  which statistics will be compiled.

C  IY_END is an integer variable that specifies the last grid point
C  number in the second dimension of the rectangular subset of A for
C  which statistics will be compiled.

C  IZ_START is an integer variable that specifies the first grid point
C  number in the third dimension of the rectangular subset of A for
C  which statistics will be compiled.

C  IZ_END is an integer variable that specifies the last grid point
C  number in the third dimension of the rectangular subset of A for
C  which statistics will be compiled.

C  BADDATA is a real variable that indicates a missing value as
C  described.

C  Output:

C  N_GOOD is an integer variable that returns the number of points at
C  which data exist in the rectangular subset of A.

C  N_BAD is an integer variable that returns the number of points at
C  which data are missing in the rectangular subset of A.

C  MEAN is a real variable that returns the mean in the rectangular
C  subset of A.  If it cannot be calculated, it is returned as BADDATA.

C  SD is a real variable that returns the sample standard deviation in
C  the rectangular subset of A.  If it cannot be calculated, it is
C  returned as BADDATA.

C  A_MIN is a real variable that returns the minimum data value in the
C  rectangular subset of A.  If it does not exist, it is returned as
C  BADDATA.

C  IX_A_MIN is an integer variable that returns the grid point number in
C  the first dimension of the minimum data value in the rectangular
C  subset of A.  If it does not exist, it is returned as 0.

C  IY_A_MIN is an integer variable that returns the grid point number in
C  the second dimension of the minimum data value in the rectangular
C  subset of A.  If it does not exist, it is returned as 0.

C  IZ_A_MIN is an integer variable that returns the grid point number in
C  the third dimension of the minimum data value in the rectangular
C  subset of A.  If it does not exist, it is returned as 0.

C  A_MAX is a real variable that returns the maximum data value in the
C  rectangular subset of A.  If it does not exist, it is returned as
C  BADDATA.

C  IX_A_MAX is an integer variable that returns the grid point number in
C  the first dimension of the maximum data value in the rectangular
C  subset of A.  If it does not exist, it is returned as 0.

C  IY_A_MAX is an integer variable that returns the grid point number in
C  the second dimension of the maximum data value in the rectangular
C  subset of A.  If it does not exist, it is returned as 0.

C  IZ_A_MAX is an integer variable that returns the grid point number in
C  the third dimension of the maximum data value in the rectangular
C  subset of A.  If it does not exist, it is returned as 0.

C  IX_MIN is an integer variable that returns the smallest grid point
C  number in the first dimension at which data exist.

C  IX_MAX is an integer variable that returns the largest grid point
C  number in the first dimension at which data exist.

C  IY_MIN is an integer variable that returns the smallest grid point
C  number in the second dimension at which data exist.

C  IY_MAX is an integer variable that returns the largest grid point
C  number in the second dimension at which data exist.

C  IZ_MIN is an integer variable that returns the smallest grid point
C  number in the third dimension at which data exist.

C  IZ_MAX is an integer variable that returns the largest grid point
C  number in the third dimension at which data exist.

      IMPLICIT NONE
      INTEGER MAXX,MAXY,MAXZ,IX_START,IX_END,IY_START,IY_END,IZ_START,
     $IZ_END,IX,IY,IZ,N_BAD,N_GOOD,IX_A_MIN,IY_A_MIN,IZ_A_MIN,IX_A_MAX,
     $IY_A_MAX,IZ_A_MAX,IX_MIN,IX_MAX,IY_MIN,IY_MAX,IZ_MIN,IZ_MAX
      REAL BADDATA,SUM_A,SUM_A_SQ,MEAN,VAR_A,SD,A_MIN,A_MAX,S_AA
      REAL A(MAXX,MAXY,MAXZ)

      SUM_A=0.
      SUM_A_SQ=0.
      N_GOOD=0
      N_BAD=0
      DO IZ=IZ_START,IZ_END
         DO IY=IY_START,IY_END
            DO IX=IX_START,IX_END
               IF(A(IX,IY,IZ).NE.BADDATA)THEN
                  N_GOOD=N_GOOD+1
                  SUM_A=SUM_A+A(IX,IY,IZ)
                  SUM_A_SQ=SUM_A_SQ+A(IX,IY,IZ)**2
                  IF(N_GOOD.EQ.1)THEN
                     A_MIN=A(IX,IY,IZ)
                     IX_A_MIN=IX
                     IY_A_MIN=IY
                     IZ_A_MIN=IZ
                     A_MAX=A(IX,IY,IZ)
                     IX_A_MAX=IX
                     IY_A_MAX=IY
                     IZ_A_MAX=IZ
                     IX_MIN=IX
                     IX_MAX=IX
                     IY_MIN=IY
                     IY_MAX=IY
                     IZ_MIN=IZ
                     IZ_MAX=IZ
                  ELSE
                     IF(A(IX,IY,IZ).LT.A_MIN)THEN
                        A_MIN=A(IX,IY,IZ)
                        IX_A_MIN=IX
                        IY_A_MIN=IY
                        IZ_A_MIN=IZ
                     ENDIF
                     IF(A(IX,IY,IZ).GT.A_MAX)THEN
                        A_MAX=A(IX,IY,IZ)
                        IX_A_MAX=IX
                        IY_A_MAX=IY
                        IZ_A_MAX=IZ
                     ENDIF
                     IF(IX.LT.IX_MIN)THEN
                        IX_MIN=IX
                     ENDIF
                     IF(IX.GT.IX_MAX)THEN
                        IX_MAX=IX
                     ENDIF
                     IF(IY.LT.IY_MIN)THEN
                        IY_MIN=IY
                     ENDIF
                     IF(IY.GT.IY_MAX)THEN
                        IY_MAX=IY
                     ENDIF
                     IF(IZ.LT.IZ_MIN)THEN
                        IZ_MIN=IZ
                     ENDIF
                     IF(IZ.GT.IZ_MAX)THEN
                        IZ_MAX=IZ
                     ENDIF
                  ENDIF
               ELSE
                  N_BAD=N_BAD+1
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      IF(N_GOOD.GE.1)THEN
         MEAN=SUM_A/FLOAT(N_GOOD)
         IF(N_GOOD.GE.2)THEN
            S_AA=SUM_A_SQ-SUM_A**2/FLOAT(N_GOOD)
            IF(S_AA.LT.0.)THEN
               S_AA=0.
            ENDIF
            VAR_A=S_AA/FLOAT(N_GOOD-1)
            SD=SQRT(VAR_A)
         ELSE
            SD=BADDATA
         ENDIF
      ELSE
         MEAN=BADDATA
         SD=BADDATA
         A_MIN=BADDATA
         IX_A_MIN=0
         IY_A_MIN=0
         IZ_A_MIN=0
         A_MAX=BADDATA
         IX_A_MAX=0
         IY_A_MAX=0
         IZ_A_MAX=0
         IX_MIN=0
         IX_MAX=0
         IY_MIN=0
         IY_MAX=0
         IZ_MIN=0
         IZ_MAX=0
      ENDIF
      END
