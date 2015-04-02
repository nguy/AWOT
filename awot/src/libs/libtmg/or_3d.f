      SUBROUTINE OR_3D(A_IN,SD_A_IN,B_IN,SD_B_IN,MAXX_IN,MAXY_IN,
     $MAXZ_IN,NX,NY,NZ,BADDATA,MAXX_OUT,MAXY_OUT,MAXZ_OUT,A_OUT,
     $SD_A_OUT)

C  Thomas Matejka NOAA/NSSL 6 April 1994

C  This subroutine creates a three-dimensional data field and its
C  standard deviation that consists at each point of the more reliable
C  value at that point in a two input fields.

C  Input:

C  A_IN is a three-dimensional real array.  A_IN(I,J,K) specifies the
C  value, at the Ith grid point in the first dimension, the Jth grid
C  point in the second dimension, and the Kth grid point in the third
C  dimension, of the first data field to be considered for copying to
C  the output field.  If it is missing, it should equal BADDATA.

C  SD_A_IN is a three-dimensional real array.  SD_A_IN(I,J,K) specifies
C  the standard deviation of A_IN(I,J,K).  If it is missing, it should
C  equal BADDATA.

C  B_IN is a three-dimensional real array.  B_IN(I,J,K) specifies the
C  value, at the Ith grid point in the first dimension, the Jth grid
C  point in the second dimension, and the Kth grid point in the third
C  dimension, of the second data field to be considered for copying to
C  the output field.  If it is missing, it should equal BADDATA.

C  SD_B_IN is a three-dimensional real array.  SD_B_IN(I,J,K) specifies
C  the standard deviation of B_IN(I,J,K).  If it is missing, it should
C  equal BADDATA.

C  MAXX_IN is an integer variable that specifies the first dimension of
C  A_IN, SD_A_IN, B_IN, and SD_B_IN in the calling program.

C  MAXY_IN is an integer variable that specifies the second dimension of
C  A_IN, SD_A_IN, B_IN, and SD_B_IN in the calling program.

C  MAXZ_IN is an integer variable that specifies the third dimension of
C  A_IN, SD_A_IN, B_IN, and SD_B_IN in the calling program.

C  NX is an integer variable that specifies the number of grid points in
C  the first dimension for A_IN, SD_A_IN, B_IN, SD_B_IN, A_OUT, and
C  SD_A_OUT.

C  NY is an integer variable that specifies the number of grid points in
C  the second dimension for A_IN, SD_A_IN, B_IN, SD_B_IN, A_OUT, and
C  SD_A_OUT.

C  NZ is an integer variable that specifies the number of grid points in
C  the third dimension for A_IN, SD_A_IN, B_IN, SD_B_IN, A_OUT, and
C  SD_A_OUT.

C  BADDATA is a real variable that indicates a missing value as
C  described.

C  MAXX_OUT is an integer variable that specifies the first dimension of
C  A_OUT and SD_A_OUT in the calling program.

C  MAXY_OUT is an integer variable that specifies the second dimension
C  of A_OUT and SD_A_OUT in the calling program.

C  MAXZ_OUT is an integer variable that specifies the third dimension of
C  A_OUT and SD_A_OUT and the calling program.

C  Output:

C  A_OUT is a three-dimensional real array.  A_OUT(I,J,K) returns the
C  value, at the Ith grid point in the first dimension, the Jth grid
C  point in the second dimension, and the Kth grid point in the third
C  dimension, of the more reliable of A_IN(I,J,K) and B_IN(I,J,K).  If
C  it is missing, it is returned as BADDATA.  A_OUT may be identical to
C  A_IN or B_IN, in which case A_IN or B_IN is overwritten.

C  SD_A_OUT is a three-dimensional real array.  SD_A_OUT(I,J,K) returns
C  the standard deviation of A_OUT(I,J,K).  If it is missing, it is
C  returned as BADDATA.  SD_A_OUT may be identical to SD_A_IN or
C  SD_B_IN, in which case SD_A_IN or SD_B_IN is overwritten.

      IMPLICIT NONE
      INTEGER MAXX_IN,MAXY_IN,MAXZ_IN,MAXX_OUT,MAXY_OUT,MAXZ_OUT,NX,NY,
     $NZ,IX,IY,IZ
      REAL BADDATA
      REAL A_IN(MAXX_IN,MAXY_IN,MAXZ_IN),
     $SD_A_IN(MAXX_IN,MAXY_IN,MAXZ_IN),B_IN(MAXX_IN,MAXY_IN,MAXZ_IN),
     $SD_B_IN(MAXX_IN,MAXY_IN,MAXZ_IN)
      REAL A_OUT(MAXX_OUT,MAXY_OUT,MAXZ_OUT),
     $SD_A_OUT(MAXX_OUT,MAXY_OUT,MAXZ_OUT)

C  Loop through all the points.
      DO 1 IZ=1,NZ
         DO 2 IY=1,NY
            DO 3 IX=1,NX

C  Find the more reliable value of A_IN and B_IN.
               IF(A_IN(IX,IY,IZ).NE.BADDATA.AND.
     $         SD_A_IN(IX,IY,IZ).NE.BADDATA)THEN
                  IF(B_IN(IX,IY,IZ).NE.BADDATA.AND.
     $            SD_B_IN(IX,IY,IZ).NE.BADDATA)THEN
                     IF(SD_A_IN(IX,IY,IZ).LE.SD_B_IN(IX,IY,IZ))THEN
                        A_OUT(IX,IY,IZ)=A_IN(IX,IY,IZ)
                        SD_A_OUT(IX,IY,IZ)=SD_A_IN(IX,IY,IZ)
                     ELSE
                        A_OUT(IX,IY,IZ)=B_IN(IX,IY,IZ)
                        SD_A_OUT(IX,IY,IZ)=SD_B_IN(IX,IY,IZ)
                     ENDIF
                  ELSE
                     A_OUT(IX,IY,IZ)=A_IN(IX,IY,IZ)
                     SD_A_OUT(IX,IY,IZ)=SD_A_IN(IX,IY,IZ)
                  ENDIF
               ELSE
                  A_OUT(IX,IY,IZ)=B_IN(IX,IY,IZ)
                  SD_A_OUT(IX,IY,IZ)=SD_B_IN(IX,IY,IZ)
               ENDIF
3           CONTINUE
2        CONTINUE
1     CONTINUE

C  Done.
      RETURN
      END
