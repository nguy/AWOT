      SUBROUTINE ABSMX_3D(A,MAXX,MAXY,MAXZ,IX_START,IX_END,IY_START,
     $IY_END,IZ_START,IZ_END,BADDATA,CUMULATIVE,A_MAX_MAG_IN,A_MAX_MAG)

C  Thomas Matejka NOAA/NSSL 11 March 1994

C  This subroutine finds the maximum absolute value in a rectangular
C  subset of a three-dimensional data field.

C  Input:

C  A is a three-dimensional real array.  A(I,J,K) specifies the value,
C  at the Ith grid point in the first dimension, the Jth grid point in
C  the second dimension, and the Kth grid point in the third dimension,
C  of the data field whose maximum absolute value is sought.  If it is
C  missing, it should equal BADDATA.

C  MAXX is an integer variable that specifies the first dimension of A
C  in the calling program.

C  MAXY is an integer variable that specifies the second dimension of A
C  in the calling program.

C  MAXZ is an integer variable that specifies the third dimension of A
C  in the calling program.

C  IX_START is an integer variable that specifies the first grid point
C  number in the first dimension of the rectangular subset of A in which
C  a maximum absolute value will be sought.

C  IX_END is an integer variable that specifies the last grid point
C  number in the first dimension of the rectangular subset of A in which
C  a maximum absolute value will be sought.

C  IY_START is an integer variable that specifies the first grid point
C  number in the second dimension of the rectangular subset of A in which
C  a maximum absolute value will be sought.

C  IY_END is an integer variable that specifies the last grid point
C  number in the second dimension of the rectangular subset of A in which
C  a maximum absolute value will be sought.

C  IZ_START is an integer variable that specifies the first grid point
C  number in the third dimension of the rectangular subset of A in which
C  a maximum absolute value will be sought.

C  IZ_END is an integer variable that specifies the last grid point
C  number in the third dimension of the rectangular subset of A in which
C  a maximum absolute value will be sought.

C  BADDATA is a real variable that indicates a missing value as
C  described.

C  CUMULATIVE is a logical variable.  If CUMULATIVE is .TRUE., then a
C  cumulative maximum absolute value can be determined by specifying an
C  initial maximum absolute value.

C  A_MAX_MAG_IN is a real variable that is relevant only when CUMULATIVE
C  is .TRUE..  A_MAX_MAG_IN is assumed initially to be the maximum
C  absolute value.

C  Output:

C  A_MAX_MAG is a real variable that returns the maximum absolute value
C  in the rectangular subset of A.  If it does not exist, it is returned
C  as BADDATA.

      IMPLICIT NONE
      LOGICAL CUMULATIVE,FOUND
      INTEGER MAXX,MAXY,MAXZ,IX_START,IX_END,IY_START,IY_END,IZ_START,
     $IZ_END,IX,IY,IZ
      REAL A_MAX_MAG,BADDATA,A_MAX_MAG_IN
      REAL A(MAXX,MAXY,MAXZ)

C  Initialize with the cumulative value if requested.
      IF(CUMULATIVE)THEN
         A_MAX_MAG=A_MAX_MAG_IN
         FOUND=.TRUE.
      ELSE
         FOUND=.FALSE.
      ENDIF

C  Loop through the points and find the maximum absolute value.
      DO 1 IZ=IZ_START,IZ_END
         DO 2 IY=IY_START,IY_END
            DO 3 IX=IX_START,IX_END
               IF(A(IX,IY,IZ).NE.BADDATA)THEN
                  IF(FOUND)THEN
                     A_MAX_MAG=AMAX1(A_MAX_MAG,ABS(A(IX,IY,IZ)))
                  ELSE     
                     A_MAX_MAG=ABS(A(IX,IY,IZ))
                     FOUND=.TRUE.
                  ENDIF
               ENDIF
3           CONTINUE
2        CONTINUE
1     CONTINUE

C  Check whether a maximum absolute value was found.
      IF(.NOT.FOUND)THEN
         A_MAX_MAG=BADDATA
      ENDIF

C  Done.
      RETURN
      END
