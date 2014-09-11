      SUBROUTINE THRESH_ON_FIELD_3D(A_IN,B,MAXX_IN,MAXY_IN,MAXZ_IN,NX,
     $NY,NZ,BADDATA,THRESH_VALUE_LOW,THRESH_VALUE_HIGH,
     $REPLACE_IF_MISSING,REPLACEMENT_VALUE_MISSING,
     $REPLACEMENT_VALUE_LOW,REPLACEMENT_VALUE_HIGH,MAXX_OUT,MAXY_OUT,
     $MAXZ_OUT,A_OUT)

C  Thomas Matejka NOAA/NSSL 1 April 1994

C  This subroutine deletes or replaces data in a three-dimensional data
C  field based on the data in another three-dimensional data field.  A
C  datum in the first field is set equal to a designated value if the
C  datum at the same point in the second field is missing.  A datum in
C  the first field is set equal to a designated value if the datum at
C  the same point in the second field is larger than a designated value.
C  A datum in the first field is set equal to a designated value if the
C  datum at the same point in the second field is smaller than a
C  designated value.  Replacement of data in the first field can be
C  performed or prevented if data in the first field are missing, as
C  desired.

C  Input:

C  A_IN is a three-dimensional real array.  A_IN(I,J,K) specifies the
C  value, at the Ith grid point in the first dimension, the Jth grid
C  point in the second dimension, and the Kth grid point in the third
C  dimension, of the data field to be thresholded.  If it is missing, it
C  should equal BADDATA.

C  B is a three-dimensional real array.  B(I,J,K) specifies the value,
C  at the Ith grid point in the first dimension, the Jth grid point in
C  the second dimension, and the Kth grid point in the third dimension,
C  of the data field governing the thresholding.  If it is missing, it
C  should equal BADDATA.  B may be identical to A, in which case A is
C  thresholded on itself.

C  MAXX_IN is an integer variable that specifies the first dimension of
C  A_IN and B in the calling program.

C  MAXY_IN is an integer variable that specifies the second dimension of
C  A_IN and B in the calling program.

C  MAXZ_IN is an integer variable that specifies the third dimension of
C  A_IN and B in the calling program.

C  NX is an integer variable that specifies the number of grid points in
C  the first dimension for A_IN, B, and A_OUT.

C  NY is an integer variable that specifies the number of grid points in
C  the second dimension for A_IN, B, and A_OUT.

C  NZ is an integer variable that specifies the number of grid points in
C  the third dimension for A_IN, B, and A_OUT.

C  BADDATA is a real variable that indicates a missing value as
C  described.

C  THRESH_VALUE_LOW is a real variable that specifies the minimum value
C  of data in field B for which data in field A are to be retained.  If
C  it equals BADDATA, then thresholding on THRESH_VALUE_LOW is not
C  performed.

C  THRESH_VALUE_HIGH is a real variable that specifies the maximum value
C  of data in field B for which data in field A are to be retained.  If
C  it equals BADDATA, then thresholding on THRESH_VALUE_LOW is not
C  performed.

C  REPLACE_IF_MISSING is a logical variable.  If REPLACE_IF_MISSING is
C  .FALSE., missing data in A_IN remain missing regardless of values in
C  B.  If REPLACE_IF_MISSING is .TRUE., missing data in A_IN are
C  replaced by the specified values according to values in B.

C  REPLACEMENT_VALUE_MISSING is a real variable that specifies the value
C  to be assigned to a point in field A when the point in field B is
C  missing.  It may equal BADDATA.

C  REPLACEMENT_VALUE_LOW is a real variable that specifies the value to
C  be assigned to a point in field A when the point in field B is less
C  than THRESH_VALUE_LOW.  It is relevant only when THRESH_VALUE_LOW is
C  not equal to BADDATA.  It may equal BADDATA.

C  REPLACEMENT_VALUE_HIGH is a real variable that specifies the value to
C  be assigned to a point in field A when the point in field B is
C  greater than THRESH_VALUE_HIGH.  It is relevant only when
C  THRESH_VALUE_HIGH is not equal to BADDATA.  It may equal BADDATA.

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
C  dimension, of the thresholded value of A(I,J,K).  If it is missing,
C  it is returned as BADDATA.  A_OUT may be identical to A, in which
C  case A is overwritten.

      IMPLICIT NONE
      LOGICAL REPLACE_IF_MISSING
      INTEGER MAXX_IN,MAXY_IN,MAXZ_IN,MAXX_OUT,MAXY_OUT,MAXZ_OUT,NX,NY,
     $NZ,IX,IY,IZ
      REAL BADDATA,THRESH_VALUE_LOW,THRESH_VALUE_HIGH,
     $REPLACEMENT_VALUE_MISSING,REPLACEMENT_VALUE_LOW,
     $REPLACEMENT_VALUE_HIGH
      REAL A_IN(MAXX_IN,MAXY_IN,MAXZ_IN),B(MAXX_IN,MAXY_IN,MAXZ_IN)
      REAL A_OUT(MAXX_OUT,MAXY_OUT,MAXZ_OUT)

C  Loop through all the points.
      DO 1 IZ=1,NZ
         DO 2 IY=1,NY
            DO 3 IX=1,NX

C  Check whether A is missing and is not to be replaced.
               IF(.NOT.REPLACE_IF_MISSING.AND.
     $         A_IN(IX,IY,IZ).EQ.BADDATA)THEN
                  A_OUT(IX,IY,IZ)=BADDATA

C  Threshold A when B is missing.
               ELSEIF(B(IX,IY,IZ).EQ.BADDATA)THEN
                  A_OUT(IX,IY,IZ)=REPLACEMENT_VALUE_MISSING

C  Threshold A on low values of B.
               ELSEIF(THRESH_VALUE_LOW.NE.BADDATA.AND.
     $         B(IX,IY,IZ).LT.THRESH_VALUE_LOW)THEN
                  A_OUT(IX,IY,IZ)=REPLACEMENT_VALUE_LOW

C  Threshold A on high values of B.
               ELSEIF(THRESH_VALUE_HIGH.NE.BADDATA.AND.
     $         B(IX,IY,IZ).GT.THRESH_VALUE_HIGH)THEN
                  A_OUT(IX,IY,IZ)=REPLACEMENT_VALUE_HIGH

C  Retain A.
               ELSE
                  A_OUT(IX,IY,IZ)=A_IN(IX,IY,IZ)
               ENDIF
3           CONTINUE
2        CONTINUE
1     CONTINUE

C  Done.
      RETURN
      END
