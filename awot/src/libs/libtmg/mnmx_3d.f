      SUBROUTINE MNMX_3D(A,MAXX,MAXY,MAXZ,IX_START,IX_END,IY_START,
     $IY_END,IZ_START,IZ_END,BADDATA,CUMULATIVE,A_MIN_IN,A_MAX_IN,A_MIN,
     $A_MAX)

C  Thomas Matejka NOAA/NSSL 11 March 1994

C  This subroutine finds the minimum and maximum values in a rectangular
C  subset of a three-dimensional data field.

C  Input:

C  A is a three-dimensional real array.  A(I,J,K) specifies the value,
C  at the Ith grid point in the first dimension, the Jth grid point in
C  the second dimension, and the Kth grid point in the third dimension,
C  of the data field whose minimum and maximum values are sought.  If it
C  is missing, it should equal BADDATA.

C  MAXX is an integer variable that specifies the first dimension of A
C  in the calling program.

C  MAXY is an integer variable that specifies the second dimension of A
C  in the calling program.

C  MAXZ is an integer variable that specifies the third dimension of A
C  in the calling program.

C  IX_START is an integer variable that specifies the first grid point
C  number in the first dimension of the rectangular subset of A in which
C  the minimum and maximum values will be found.

C  IX_END is an integer variable that specifies the last grid point
C  number in the first dimension of the rectangular subset of A in which
C  the minimum and maximum values will be found.

C  IY_START is an integer variable that specifies the first grid point
C  number in the second dimension of the rectangular subset of A in
C  which the minimum and maximum values will be found.

C  IY_END is an integer variable that specifies the last grid point
C  number in the second dimension of the rectangular subset of A in
C  which the minimum and maximum values will be found.

C  IZ_START is an integer variable that specifies the first grid point
C  number in the third dimension of the rectangular subset of A in
C  which the minimum and maximum values will be found.

C  IZ_END is an integer variable that specifies the last grid point
C  number in the third dimension of the rectangular subset of A in
C  which the minimum and maximum values will be found.

C  BADDATA is a real variable that indicates a missing value as
C  described.

C  CUMULATIVE is a logical variable.  If CUMULATIVE is .TRUE., then
C  cumulative minimum and maximum values can be determined by specifying
C  initial minimum and maximum values.

C  A_MIN_IN is a real variable that is relevant only when CUMULATIVE is
C  .TRUE..  A_MIN_IN is assumed initially to be the minimum value.

C  A_MAX_IN is a real variable that is relevant only when CUMULATIVE is
C  .TRUE..  A_MAX_IN is assumed initially to be the maximum value.

C  Output:

C  A_MIN is a real variable that returns the minimum value in the
C  rectangular subset of A.  If it does not exist, it is returned as
C  BADDATA.

C  A_MAX is a real variable that returns the maximum value in the
C  rectangular subset of A.  If it does not exist, it is returned as
C  BADDATA.

      IMPLICIT NONE
      LOGICAL CUMULATIVE,FOUND
      INTEGER MAXX,MAXY,MAXZ,IX_START,IX_END,IY_START,IY_END,IZ_START,
     $IZ_END,IX,IY,IZ
      REAL A_MIN,A_MAX,BADDATA,A_MIN_IN,A_MAX_IN
      REAL A(MAXX,MAXY,MAXZ)

C  Initialize with the cumulative values if requested.
      IF(CUMULATIVE)THEN
         A_MIN=A_MIN_IN
         A_MAX=A_MAX_IN
         FOUND=.TRUE.
      ELSE
         FOUND=.FALSE.
      ENDIF

C  Loop through the points and find the minimum and maximum values.
      DO 1 IZ=IZ_START,IZ_END
         DO 2 IY=IY_START,IY_END
            DO 3 IX=IX_START,IX_END
               IF(A(IX,IY,IZ).NE.BADDATA)THEN
                  IF(FOUND)THEN
                     A_MIN=AMIN1(A_MIN,A(IX,IY,IZ))
                     A_MAX=AMAX1(A_MAX,A(IX,IY,IZ))
                  ELSE     
                     A_MIN=A(IX,IY,IZ)
                     A_MAX=A(IX,IY,IZ)
                     FOUND=.TRUE.
                  ENDIF
               ENDIF
3           CONTINUE
2        CONTINUE
1     CONTINUE

C  Check whether minimum and maximum values were found.
      IF(.NOT.FOUND)THEN
         A_MIN=BADDATA
         A_MAX=BADDATA
      ENDIF

C  Done.
      RETURN
      END
