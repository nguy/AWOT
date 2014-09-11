      SUBROUTINE DDY_3D(A,MAXX_IN,MAXY_IN,MAXZ_IN,NX,NY,NZ,DELY,BADDATA,
     $MAXX_OUT,MAXY_OUT,MAXZ_OUT,DDY_A)

C  Thomas Matejka NOAA/NSSL 24 January 1994

C  This subroutine calculates the partial derivative with respect to the
C  second dimension of a three-dimensional data field.  A centered
C  difference over two grid spaces is used where possible.  At data
C  edges, a two-point difference over one grid space is used where
C  possible.

C  Input:

C  A is a three-dimensional real array.  A(I,J,K) specifies the value,
C  at the Ith grid point in the first dimension, the Jth grid point in
C  the second dimension, and the Kth grid point in the third dimension,
C  of the data field whose derivative is sought.  If it is missing, it
C  should equal BADDATA.

C  MAXX_IN is an integer variable that specifies the first dimension of
C  A in the calling program.

C  MAXY_IN is an integer variable that specifies the second dimension of
C  A in the calling program.

C  MAXZ_IN is an integer variable that specifies the third dimension of
C  A in the calling program.

C  NX is an integer variable that specifies the number of grid points in
C  the first dimension for A and DDY_A.

C  NY is an integer variable that specifies the number of grid points in
C  the second dimension for A and DDY_A.

C  NZ is an integer variable that specifies the number of grid points in
C  the third dimension for A and DDY_A.

C  DELY is a real variable that specifies the grid spacing in the second
C  dimension for A and DDY_A.

C  BADDATA is a real variable that specifies the value used to indicate
C  a missing value as described.

C  MAXX_OUT is an integer variable that specifies the first dimension of
C  DDY_A in the calling program.

C  MAXY_OUT is an integer variable that specifies the second dimension
C  of DDY_A in the calling program.

C  MAXZ_OUT is an integer variable that specifies the third dimension of
C  DDY_A in the calling program.

C  Output:

C  DDY_A is a three-dimensional real array.  DDY_A(I,J,K) returns the
C  value, at the Ith grid point in the first dimension, the Jth grid
C  point in the second dimension, and the Kth grid point in the third
C  dimension, of the partial derivative of A with respect to the second
C  dimension.  If it is missing, it is returned as BADDATA.

      IMPLICIT NONE
      INTEGER MAXX_IN,MAXY_IN,MAXZ_IN,MAXX_OUT,MAXY_OUT,MAXZ_OUT,NX,NY,
     $NZ,IX,IY,IZ
      REAL DELY,BADDATA
      REAL A(MAXX_IN,MAXY_IN,MAXZ_IN)
      REAL DDY_A(MAXX_OUT,MAXY_OUT,MAXZ_OUT)

      DO 1 IZ=1,NZ
         DO 2 IY=1,NY
            DO 3 IX=1,NX
               IF(IY.GT.1.AND.
     $         A(IX,IY-1,IZ).NE.BADDATA)THEN
                  IF(IY.LT.NY.AND.
     $            A(IX,IY+1,IZ).NE.BADDATA)THEN
                     DDY_A(IX,IY,IZ)=(A(IX,IY+1,IZ)-A(IX,IY-1,IZ))/
     $               (2.*DELY)
                  ELSEIF(A(IX,IY,IZ).NE.BADDATA)THEN
                     DDY_A(IX,IY,IZ)=(A(IX,IY,IZ)-A(IX,IY-1,IZ))/DELY
                  ELSE
                     DDY_A(IX,IY,IZ)=BADDATA
                  ENDIF
               ELSEIF(IY.LT.NY.AND.
     $         A(IX,IY+1,IZ).NE.BADDATA.AND.
     $         A(IX,IY,IZ).NE.BADDATA)THEN
                  DDY_A(IX,IY,IZ)=(A(IX,IY+1,IZ)-A(IX,IY,IZ))/DELY
               ELSE
                  DDY_A(IX,IY,IZ)=BADDATA
               ENDIF
3           CONTINUE
2        CONTINUE
1     CONTINUE
      RETURN
      END              
