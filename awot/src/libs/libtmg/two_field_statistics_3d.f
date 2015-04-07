      SUBROUTINE TWO_FIELD_STATISTICS_3D(A,B,MAXX,MAXY,MAXZ,IX_START,
     $IX_END,IY_START,IY_END,IZ_START,IZ_END,BADDATA,N_GOOD,N_BAD,
     $SQRT_COV,CORRELATION,A_MINUS_B_MIN,IX_MIN,IY_MIN,IZ_MIN,
     $A_MINUS_B_MAX,IX_MAX,IY_MAX,IZ_MAX,A_MINUS_B_MAX_MAG,IX_MAX_MAG,
     $IY_MAX_MAG,IZ_MAX_MAG)

C  Thomas Matejka NOAA/NSSL 24 April 1997

C  This subroutine calculates statistics for a rectangular subset of two
C  three-dimensional fields.

C  Input:

C  A is a three-dimensional real array.  A(I,J,K) specifies the value,
C  at the Ith grid point in the first dimension, the Jth grid point in
C  the second dimension, and the Kth grid point in the third dimension,
C  of the first data field whose difference statistics are sought.  If
C  it is missing, it should equal BADDATA.

C  B is a three-dimensional real array.  B(I,J,K) specifies the value,
C  at the Ith grid point in the first dimension, the Jth grid point in
C  the second dimension, and the Kth grid point in the third dimension,
C  of the second data field whose difference statistics are sought.  If
C  it is missing, it should equal BADDATA.

C  MAXX is an integer variable that specifies the first dimension of A
C  and B in the calling program.

C  MAXY is an integer variable that specifies the second dimension of A
C  and B in the calling program.

C  MAXZ is an integer variable that specifies the third dimension of A
C  and B in the calling program.

C  IX_START is an integer variable that specifies the first grid point
C  number in the first dimension of the rectangular subset of A and B
C  for which statistics will be compiled.

C  IX_END is an integer variable that specifies the last grid point
C  number in the first dimension of the rectangular subset of A and B
C  for which statistics will be compiled.

C  IY_START is an integer variable that specifies the first grid point
C  number in the second dimension of the rectangular subset of A and B
C  for which statistics will be compiled.

C  IY_END is an integer variable that specifies the last grid point
C  number in the second dimension of the rectangular subset of A and B
C  for which statistics will be compiled.

C  IZ_START is an integer variable that specifies the first grid point
C  number in the third dimension of the rectangular subset of A and B
C  for which statistics will be compiled.

C  IZ_END is an integer variable that specifies the last grid point
C  number in the third dimension of the rectangular subset of A and B
C  for which statistics will be compiled.

C  BADDATA is a real variable that indicates a missing value as
C  described.

C  Output:

C  N_GOOD is an integer variable that returns the number of points at
C  which data exist in both of the rectangular subsets of A and B.

C  N_BAD is an integer variable that returns the number of points at
C  which data are missing in either of the rectangular subsets of A or
C  B.

C  SQRT_COV is a real variable that returns the square root of the
C  sample covariance of the two fields.

C  CORRELATION is a real variable that returns the correlation between
C  the two fields.

C  A_MINUS_B_MIN is a real variable that returns the minimum value of A
C  minus B in the rectangular subsets of A and B.  If it does not exist,
C  it is returned as BADDATA.

C  IX_MIN is an integer variable that returns the grid point number in
C  the first dimension of the minimum value of A minus B in the
C  rectangular subsets of A and B.  If it does not exist, it is returned
C  as 0.

C  IY_MIN is an integer variable that returns the grid point number in
C  the second dimension of the minimum value of A minus B in the
C  rectangular subsets of A and B.  If it does not exist, it is returned
C  as 0.

C  IZ_MIN is an integer variable that returns the grid point number in
C  the third dimension of the minimum value of A minus B in the
C  rectangular subsets of A and B.  If it does not exist, it is returned
C  as 0.

C  A_MINUS_B_MAX is a real variable that returns the maximum value of A
C  minus B in the rectangular subsets of A and B.  If it does not exist,
C  it is returned as BADDATA.

C  IX_MAX is an integer variable that returns the grid point number in
C  the first dimension of the maximum value of A minus B in the
C  rectangular subsets of A and B.  If it does not exist, it is returned
C  as 0.

C  IY_MAX is an integer variable that returns the grid point number in
C  the second dimension of the maximum value of A minus B in the
C  rectangular subsets of A and B.  If it does not exist, it is returned
C  as 0.

C  IZ_MAX is an integer variable that returns the grid point number in
C  the third dimension of the maximum value of A minus B in the
C  rectangular subsets of A and B.  If it does not exist, it is returned
C  as 0.

C  A_MINUS_B_MAX_MAG is a real variable that returns the magnitude of
C  the maximum difference between A and B in the rectangular subsets of
C  A and B.  If it does not exist, it is returned as BADDATA.

C  IX_MAX_MAG is an integer variable that returns the grid point number
C  in the first dimension of the maximum difference between A and B in
C  the rectangular subsets of A and B.  If it does not exist, it is
C  returned as 0.

C  IY_MAX_MAG is an integer variable that returns the grid point number
C  in the second dimension of the maximum difference between A and B in
C  the rectangular subsets of A and B.  If it does not exist, it is
C  returned as 0.

C  IZ_MAX_MAG is an integer variable that returns the grid point number
C  in the third dimension of the maximum difference between A and B in
C  the rectangular subsets of A and B.  If it does not exist, it is
C  returned as 0.

      IMPLICIT NONE
      INTEGER MAXX,MAXY,MAXZ,IX_START,IX_END,IY_START,IY_END,IZ_START,
     $IZ_END,IX,IY,IZ,N_BAD,N_GOOD,IX_MIN,IY_MIN,IZ_MIN,IX_MAX,IY_MAX,
     $IZ_MAX,IX_MAX_MAG,IY_MAX_MAG,IZ_MAX_MAG,N_GOOD_TEMP
      REAL BADDATA,SQRT_COV,CORRELATION,A_MINUS_B_MIN,A_MINUS_B_MAX,
     $A_MINUS_B_MAX_MAG,SUM_A,SUM_B,SUM_AA,SUM_BB,SUM_AB,COV_AB,S_AB,
     $S_AA,S_BB
      REAL A(MAXX,MAXY,MAXZ),B(MAXX,MAXY,MAXZ)

      SUM_A=0.
      SUM_B=0.
      SUM_AB=0.
      SUM_AA=0.
      SUM_BB=0.
      N_GOOD_TEMP=0
      N_BAD=0
      DO IZ=IZ_START,IZ_END
         DO IY=IY_START,IY_END
            DO IX=IX_START,IX_END
               IF(A(IX,IY,IZ).NE.BADDATA.AND.
     $         B(IX,IY,IZ).NE.BADDATA)THEN
                  N_GOOD_TEMP=N_GOOD_TEMP+1
                  SUM_A=SUM_A+A(IX,IY,IZ)
                  SUM_B=SUM_B+B(IX,IY,IZ)
                  SUM_AB=SUM_AB+A(IX,IY,IZ)*B(IX,IY,IZ)
                  SUM_AA=SUM_AA+A(IX,IY,IZ)**2
                  SUM_BB=SUM_BB+B(IX,IY,IZ)**2
                  IF(N_GOOD_TEMP.EQ.1)THEN
                     A_MINUS_B_MIN=A(IX,IY,IZ)-B(IX,IY,IZ)
                     IX_MIN=IX
                     IY_MIN=IY
                     IZ_MIN=IZ
                     A_MINUS_B_MAX=A(IX,IY,IZ)-B(IX,IY,IZ)
                     IX_MAX=IX
                     IY_MAX=IY
                     IZ_MAX=IZ
                     A_MINUS_B_MAX_MAG=ABS(A(IX,IY,IZ)-B(IX,IY,IZ))
                     IX_MAX_MAG=IX
                     IY_MAX_MAG=IY
                     IZ_MAX_MAG=IZ
                  ELSE
                     IF(A(IX,IY,IZ)-B(IX,IY,IZ).LT.A_MINUS_B_MIN)THEN
                        A_MINUS_B_MIN=A(IX,IY,IZ)-B(IX,IY,IZ)
                        IX_MIN=IX
                        IY_MIN=IY
                        IZ_MIN=IZ
                     ENDIF
                     IF(A(IX,IY,IZ)-B(IX,IY,IZ).GT.A_MINUS_B_MAX)THEN
                        A_MINUS_B_MAX=A(IX,IY,IZ)-B(IX,IY,IZ)
                        IX_MAX=IX
                        IY_MAX=IY
                        IZ_MAX=IZ
                     ENDIF
                     IF(ABS(A(IX,IY,IZ)-B(IX,IY,IZ)).GT.
     $               A_MINUS_B_MAX_MAG)THEN
                        A_MINUS_B_MAX_MAG=ABS(A(IX,IY,IZ)-B(IX,IY,IZ))
                        IX_MAX_MAG=IX
                        IY_MAX_MAG=IY
                        IZ_MAX_MAG=IZ
                     ENDIF
                  ENDIF
               ELSE
                  N_BAD=N_BAD+1
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      IF(N_GOOD_TEMP.GE.2)THEN
         S_AB=SUM_AB-SUM_A*SUM_B/FLOAT(N_GOOD_TEMP)
         IF(S_AB.LT.0.)THEN
            S_AB=0.
         ENDIF
         S_AA=SUM_AA-SUM_A**2/FLOAT(N_GOOD_TEMP)
         IF(S_AA.LT.0.)THEN
            S_AA=0.
         ENDIF
         S_BB=SUM_BB-SUM_B**2/FLOAT(N_GOOD_TEMP)
         IF(S_BB.LT.0.)THEN
            S_BB=0.
         ENDIF
         COV_AB=S_AB/FLOAT(N_GOOD_TEMP-1)
         SQRT_COV=SQRT(COV_AB)
         IF(S_AA.NE.0..AND.
     $   S_BB.NE.0.)THEN
            CORRELATION=S_AB/SQRT(S_AA*S_BB)
         ELSE
            CORRELATION=BADDATA
         ENDIF
      ELSE
         SQRT_COV=BADDATA
         CORRELATION=BADDATA
         A_MINUS_B_MIN=BADDATA
         IX_MIN=0
         IY_MIN=0
         IZ_MIN=0
         A_MINUS_B_MAX=BADDATA
         IX_MAX=0
         IY_MAX=0
         IZ_MAX=0
         A_MINUS_B_MAX_MAG=BADDATA
         IX_MAX_MAG=0
         IY_MAX_MAG=0
         IZ_MAX_MAG=0
      ENDIF
      N_GOOD=N_GOOD_TEMP
      END
