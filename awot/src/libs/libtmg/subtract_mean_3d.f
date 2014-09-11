      SUBROUTINE SUBTRACT_MEAN_3D(A_IN,MAXX_IN,MAXY_IN,MAXZ_IN,NX,NY,NZ,
     $BADDATA,MAXX_OUT,MAXY_OUT,MAXZ_OUT,A_OUT)

C  Thomas Matejka NOAA/NSSL 16 May 1996

C  This subroutine creates a three-dimensional data field that consists
C  of a data field from which its mean has been subtracted.

C  Input:

C  A_IN is a three-dimensional real array.  A_IN(I,J,K) specifies the
C  value, at the Ith grid point in the first dimension, the Jth grid
C  point in the second dimension, and the Kth grid point in the third
C  dimension, of the data field whose mean is to be subtracted.  If it
C  is missing, it should equal BADDATA.

C  MAXX_IN is an integer variable that specifies the first dimension of
C  A_IN in the calling program.

C  MAXY_IN is an integer variable that specifies the second dimension of
C  A_IN in the calling program.

C  MAXZ_IN is an integer variable that specifies the third dimension of
C  A_IN in the calling program.

C  NX is an integer variable that specifies the number of grid points in
C  the first dimension for A_IN and A_OUT.

C  NY is an integer variable that specifies the number of grid points in
C  the second dimension for A_IN and A_OUT.

C  NZ is an integer variable that specifies the number of grid points in
C  the third dimension for A_IN and A_OUT.

C  BADDATA is a real variable that indicates a missing value as
C  described.

C  MAXX_OUT is an integer variable that specifies the first dimension of
C  A_OUT in the calling program.

C  MAXY_OUT is an integer variable that specifies the second dimension
C  of A_OUT in the calling program.

C  MAXZ_OUT is an integer variable that specifies the third dimension of
C  A_OUT and the calling program.

C  Output:

C  A_OUT is a three-dimensional real array.  A_OUT(I,J,K) returns the
C  value, at the Ith grid point in the first dimension, the Jth grid
C  point in the second dimension, and the Kth grid point in the third
C  dimension, of A_IN(I,J,K) less the mean of A_IN.  If it is missing,
C  it is returned as BADDATA.  A_OUT may be identical to A_IN, in which
C  case A_IN is overwritten.

      IMPLICIT NONE
      INTEGER::MAXX_IN,MAXY_IN,MAXZ_IN,
     $MAXX_OUT,MAXY_OUT,MAXZ_OUT,
     $NX,NY,NZ,IX,IY,IZ,
     $COUNT
      REAL BADDATA,
     $SUM,MEAN
      REAL A_IN(MAXX_IN,MAXY_IN,MAXZ_IN)
      REAL A_OUT(MAXX_OUT,MAXY_OUT,MAXZ_OUT)

C  Loop through all the points and compute the mean.
      SUM=0.
      COUNT=0
      DO IZ=1,NZ
         DO IY=1,NY
            DO IX=1,NX
               IF(A_IN(IX,IY,IZ).NE.BADDATA)THEN
                  SUM=SUM+A_IN(IX,IY,IZ)
                  COUNT=COUNT+1
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      IF(COUNT.NE.0)THEN
         MEAN=SUM/FLOAT(COUNT)
      ENDIF

C  Subtract the mean.
      DO IZ=1,NZ
         DO IY=1,NY
            DO IX=1,NX
               IF(A_IN(IX,IY,IZ).NE.BADDATA)THEN
                  A_OUT(IX,IY,IZ)=A_IN(IX,IY,IZ)-MEAN
               ELSE
                  A_OUT(IX,IY,IZ)=BADDATA
               ENDIF
            ENDDO
         ENDDO
      ENDDO

C  Done.
      END
