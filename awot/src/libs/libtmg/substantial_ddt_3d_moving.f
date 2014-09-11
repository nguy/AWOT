      SUBROUTINE SUBSTANTIAL_DDT_3D_MOVING(A1,A2,A3,U2,V2,W2,MAXX_IN,
     $MAXY_IN,MAXZ_IN,XMIN1,YMIN1,ZMIN1,DELX1,DELY1,DELZ1,NX1,NY1,NZ1,
     $XMIN2,YMIN2,ZMIN2,DELX2,DELY2,DELZ2,NX2,NY2,NZ2,XMIN3,YMIN3,ZMIN3,
     $DELX3,DELY3,DELZ3,NX3,NY3,NZ3,T1,T2,T3,BADDATA,USYS,VSYS,WSYS,
     $MAXX_OUT,MAXY_OUT,MAXZ_OUT,SUBSTANTIAL_DDT_A2,DDT_A2,ADVX_A2,
     $ADVY_A2,ADVZ_A2)

C  Thomas Matejka NOAA/NSSL 22 April 1994

C  This subroutine calculates the substantial derivative and the terms
C  that compose it of a three-dimensional data field input at three
C  times.  The data spacing in time does not have to be uniform.  The
C  substantial derivative may be calculated with respect to a moving
C  frame of reference.

C  Input:

C  A1 is a three-dimensional real array.  A1(I,J,K) specifies the value,
C  at the Ith grid point in the first dimension, the Jth grid point in
C  the second dimension, the Kth grid point in the third dimension, and
C  the first time, of the data field whose substantial derivative is
C  sought.  If it is missing, it should equal BADDATA.

C  A2 is a three-dimensional real array.  A2(I,J,K) specifies the value,
C  at the Ith grid point in the first dimension, the Jth grid point in
C  the second dimension, the Kth grid point in the third dimension, and
C  the second time, of the data field whose derivative is sought.  If it
C  is missing, it should equal BADDATA.

C  A3 is a three-dimensional real array.  A3(I,J,K) specifies the value,
C  at the Ith grid point in the first dimension, the Jth grid point in
C  the second dimension, the Kth grid point in the third dimension, and
C  the third time, of the data field whose derivative is sought.  If it
C  is missing, it should equal BADDATA.

C  U2 is a three-dimensional real array.  U2(I,J,K) specifies the value,
C  at the Ith grid point in the first dimension, the Jth grid point in
C  the second dimension, the Kth grid point in the third dimension, and
C  the second time, of the velocity component in the first dimension.
C  If it is missing, it should equal BADDATA.

C  V2 is a three-dimensional real array.  V2(I,J,K) specifies the value,
C  at the Ith grid point in the first dimension, the Jth grid point in
C  the second dimension, the Kth grid point in the third dimension, and
C  the second time, of the velocity component in the second dimension.
C  If it is missing, it should equal BADDATA.

C  W2 is a three-dimensional real array.  W2(I,J,K) specifies the value,
C  at the Ith grid point in the first dimension, the Jth grid point in
C  the second dimension, the Kth grid point in the third dimension, and
C  the second time, of the velocity component in the third dimension.
C  If it is missing, it should equal BADDATA.

C  MAXX_IN is an integer variable that specifies the first dimension of
C  A1, A2, A3, U2, V2, and W2 in the calling program.

C  MAXY_IN is an integer variable that specifies the second dimension of
C  A1, A2, A3, U2, V2, and W2 in the calling program.

C  MAXZ_IN is an integer variable that specifies the third dimension of
C  A1, A2, A3, U2, V2, and W2 in the calling program.

C  XMIN1 is a real variable that specifies the minimum coordinate in the
C  first dimension for A1.

C  YMIN1 is a real variable that specifies the minimum coordinate in the
C  second dimension for A1.

C  ZMIN1 is a real variable that specifies the minimum coordinate in the
C  third dimension for A1.

C  DELX1 is a real variable that specifies the grid spacing in the first
C  dimension for A1.

C  DELY1 is a real variable that specifies the grid spacing in the
C  second dimension for A1.

C  DELZ1 is a real variable that specifies the grid spacing in the third
C  dimension for A1.

C  NX1 is an integer variable that specifies the number of grid points
C  in the first dimension for A1.

C  NY1 is an integer variable that specifies the number of grid points
C  in the second dimension for A1.

C  NZ1 is an integer variable that specifies the number of grid points
C  in the third dimension for A1.

C  XMIN2 is a real variable that specifies the minimum coordinate in the
C  first dimension for A2, U2, V2, W2, SUBSTANTIAL_DDT_A2, DDT_A2,
C  ADVX_A2, ADVY_A2, and ADVZ_A2.

C  YMIN2 is a real variable that specifies the minimum coordinate in the
C  second dimension for A2, U2, V2, W2, SUBSTANTIAL_DDT_A2, DDT_A2,
C  ADVX_A2, ADVY_A2, and ADVZ_A2.

C  ZMIN2 is a real variable that specifies the minimum coordinate in the
C  third dimension for A2, U2, V2, W2, SUBSTANTIAL_DDT_A2, DDT_A2,
C  ADVX_A2, ADVY_A2, and ADVZ_A2.

C  DELX2 is a real variable that specifies the grid spacing in the first
C  dimension for A2, U2, V2, W2, SUBSTANTIAL_DDT_A2, DDT_A2, ADVX_A2,
C  ADVY_A2, and ADVZ_A2.

C  DELY2 is a real variable that specifies the grid spacing in the
C  second dimension for A2, U2, V2, W2, SUBSTANTIAL_DDT_A2, DDT_A2,
C  ADVX_A2, ADVY_A2, and ADVZ_A2.

C  DELZ2 is a real variable that specifies the grid spacing in the third
C  dimension for A2, U2, V2, W2, SUBSTANTIAL_DDT_A2, DDT_A2, ADVX_A2,
C  ADVY_A2, and ADVZ_A2.

C  NX2 is an integer variable that specifies the number of grid points
C  in the first dimension for A2, U2, V2, W2, SUBSTANTIAL_DDT_A2,
C  DDT_A2, ADVX_A2, ADVY_A2, and ADVZ_A2.

C  NY2 is an integer variable that specifies the number of grid points
C  in the second dimension for A2, U2, V2, W2, SUBSTANTIAL_DDT_A2,
C  DDT_A2, ADVX_A2, ADVY_A2, and ADVZ_A2.

C  NZ2 is an integer variable that specifies the number of grid points
C  in the third dimension for A2, U2, V2, W2, SUBSTANTIAL_DDT_A2,
C  DDT_A2, ADVX_A2, ADVY_A2, and ADVZ_A2.

C  XMIN3 is a real variable that specifies the minimum coordinate in the
C  first dimension for A3.

C  YMIN3 is a real variable that specifies the minimum coordinate in the
C  second dimension for A3.

C  ZMIN3 is a real variable that specifies the minimum coordinate in the
C  third dimension for A3.

C  DELX3 is a real variable that specifies the grid spacing in the first
C  dimension for A3.

C  DELY3 is a real variable that specifies the grid spacing in the
C  second dimension for A3.

C  DELZ3 is a real variable that specifies the grid spacing in the third
C  dimension for A3.

C  NX3 is an integer variable that specifies the number of grid points
C  in the first dimension for A3.

C  NY3 is an integer variable that specifies the number of grid points
C  in the second dimension for A3.

C  NZ3 is an integer variable that specifies the number of grid points
C  in the third dimension for A3.

C  T1 is a real variable that specifies the time for A1.

C  T2 is a real variable that specifies the time for A2, U2, V2, W2, and
C  SUBSTANTIAL_DDT_A2, DDT_A2, ADVX_A2, ADVY_A2, and ADVZ_A2.

C  T3 is a real variable that specifies the time for A3.

C  BADDATA is a real variable that specifies the value used to indicate
C  a missing value as described.

C  USYS is a real variable that specifies the component of velocity of
C  the frame of reference in the first dimension.

C  VSYS is a real variable that specifies the component of velocity of
C  the frame of reference in the second dimension.

C  WSYS is a real variable that specifies the component of velocity of
C  the frame of reference in the third dimension.

C  MAXX_OUT is an integer variable that specifies the first dimension of
C  SUBSTANTIAL_DDT_A2, DDT_A2, ADVX_A2, ADVY_A2, and ADVZ_A2 in the
C  calling program.

C  MAXY_OUT is an integer variable that specifies the second dimension
C  of SUBSTANTIAL_DDT_A2, DDT_A2, ADVX_A2, ADVY_A2, and ADVZ_A2 in the
C  calling program.

C  MAXZ_OUT is an integer variable that specifies the third dimension of
C  SUBSTANTIAL_DDT_A2, DDT_A2, ADVX_A2, ADVY_A2, and ADVZ_A2 in the
C  calling program.

C  Output:

C  SUBSTANTIAL_DDT_A2 is a three-dimensional real array.
C  SUBSTANTIAL_DDT_A2(I,J,K) returns the value, at the Ith grid point in
C  the first dimension, the Jth grid point in the second dimension, the
C  Kth grid point in the third dimension, and the second time, of the
C  substantial derivative of the field A1, A2, A3.  If it is missing, it
C  is returned as BADDATA.

C  DDT_A2 is a three-dimensional real array.  DDT_A2(I,J,K) returns the
C  value, at the Ith grid point in the first dimension, the Jth grid
C  point in the second dimension, the Kth grid point in the third
C  dimension, and the second time, of the time derivative component of
C  SUBSTANTIAL_DDT_A2(I,J,K).  If SUBSTANTIAL_DDT_A2(I,J,K) is missing,
C  it is returned as BADDATA.

C  ADVX_A2 is a three-dimensional real array.  ADVX_A2(I,J,K) returns
C  the value, at the Ith grid point in the first dimension, the Jth grid
C  point in the second dimension, the Kth grid point in the third
C  dimension, and the second time, of the first-dimensional advective
C  component of SUBSTANTIAL_DDT_A2(I,J,K).  If SUBSTANTIAL_DDT_A2(I,J,K)
C  is missing, it is returned as BADDATA.

C  ADVY_A2 is a three-dimensional real array.  ADVY_A2(I,J,K) returns
C  the value, at the Ith grid point in the first dimension, the Jth grid
C  point in the second dimension, the Kth grid point in the third
C  dimension, and the second time, of the second-dimensional advective
C  component of SUBSTANTIAL_DDT_A2(I,J,K).  If SUBSTANTIAL_DDT_A2(I,J,K)
C  is missing, it is returned as BADDATA.

C  ADVZ_A2 is a three-dimensional real array.  ADVZ_A2(I,J,K) returns
C  the value, at the Ith grid point in the first dimension, the Jth grid
C  point in the second dimension, the Kth grid point in the third
C  dimension, and the second time, of the third-dimensional advective
C  component of SUBSTANTIAL_DDT_A2(I,J,K).  If SUBSTANTIAL_DDT_A2(I,J,K)
C  is missing, it is returned as BADDATA.

      IMPLICIT NONE
      INTEGER MAXX_IN,MAXY_IN,MAXZ_IN,MAXX_OUT,MAXY_OUT,MAXZ_OUT,NX1,
     $NY1,NZ1,NX2,NY2,NZ2,NX3,NY3,NZ3,IX,IY,IZ
      REAL BADDATA,T1,T2,T3,XMIN1,YMIN1,ZMIN1,XMIN2,YMIN2,ZMIN2,XMIN3,
     $YMIN3,ZMIN3,DELX1,DELY1,DELZ1,DELX2,DELY2,DELZ2,DELX3,DELY3,DELZ3,
     $USYS,VSYS,WSYS
      REAL A1(MAXX_IN,MAXY_IN,MAXZ_IN),A2(MAXX_IN,MAXY_IN,MAXZ_IN),
     $A3(MAXX_IN,MAXY_IN,MAXZ_IN),U2(MAXX_IN,MAXY_IN,MAXZ_IN),
     $V2(MAXX_IN,MAXY_IN,MAXZ_IN),W2(MAXX_IN,MAXY_IN,MAXZ_IN)
      REAL SUBSTANTIAL_DDT_A2(MAXX_OUT,MAXY_OUT,MAXZ_OUT),
     $DDT_A2(MAXX_OUT,MAXY_OUT,MAXZ_OUT),
     $ADVX_A2(MAXX_OUT,MAXY_OUT,MAXZ_OUT),
     $ADVY_A2(MAXX_OUT,MAXY_OUT,MAXZ_OUT),
     $ADVZ_A2(MAXX_OUT,MAXY_OUT,MAXZ_OUT)
      REAL DDX_A2(NX2,NY2,NZ2),DDY_A2(NX2,NY2,NZ2),DDZ_A2(NX2,NY2,NZ2)

C  Calculate the partial derivatives with respect to all four dimensions
C  at the second time of the data field A1, A2, and A3.
      CALL DDX_3D(A2,MAXX_IN,MAXY_IN,MAXZ_IN,NX2,NY2,NZ2,DELX2,BADDATA,
     $NX2,NY2,NZ2,DDX_A2)
      CALL DDY_3D(A2,MAXX_IN,MAXY_IN,MAXZ_IN,NX2,NY2,NZ2,DELY2,BADDATA,
     $NX2,NY2,NZ2,DDY_A2)
      CALL DDZ_3D(A2,MAXX_IN,MAXY_IN,MAXZ_IN,NX2,NY2,NZ2,DELZ2,BADDATA,
     $NX2,NY2,NZ2,DDZ_A2)
      CALL DDT_3D_MOVING(A1,A2,A3,MAXX_IN,MAXY_IN,MAXZ_IN,XMIN1,YMIN1,
     $ZMIN1,DELX1,DELY1,DELZ1,NX1,NY1,NZ1,XMIN2,YMIN2,ZMIN2,DELX2,DELY2,
     $DELZ2,NX2,NY2,NZ2,XMIN3,YMIN3,ZMIN3,DELX3,DELY3,DELZ3,NX3,NY3,NZ3,
     $T1,T2,T3,BADDATA,USYS,VSYS,WSYS,MAXX_OUT,MAXY_OUT,MAXZ_OUT,DDT_A2)

C  Loop through the points in the volume and calculate the advection
C  terms and the substantial derivative.
      DO 1 IZ=1,NZ2
         DO 2 IY=1,NY2
            DO 3 IX=1,NX2
               IF(U2(IX,IY,IZ).NE.BADDATA.AND.
     $         V2(IX,IY,IZ).NE.BADDATA.AND.
     $         W2(IX,IY,IZ).NE.BADDATA.AND.
     $         DDT_A2(IX,IY,IZ).NE.BADDATA.AND.
     $         DDX_A2(IX,IY,IZ).NE.BADDATA.AND.
     $         DDY_A2(IX,IY,IZ).NE.BADDATA.AND.
     $         DDZ_A2(IX,IY,IZ).NE.BADDATA)THEN
                  ADVX_A2(IX,IY,IZ)=(U2(IX,IY,IZ)-USYS)*DDX_A2(IX,IY,IZ)
                  ADVY_A2(IX,IY,IZ)=(V2(IX,IY,IZ)-VSYS)*DDY_A2(IX,IY,IZ)
                  ADVZ_A2(IX,IY,IZ)=(W2(IX,IY,IZ)-WSYS)*DDZ_A2(IX,IY,IZ)
                  SUBSTANTIAL_DDT_A2(IX,IY,IZ)=DDT_A2(IX,IY,IZ)+
     $            ADVX_A2(IX,IY,IZ)+ADVY_A2(IX,IY,IZ)+ADVZ_A2(IX,IY,IZ)
               ELSE
                  DDT_A2(IX,IY,IZ)=BADDATA
                  ADVX_A2(IX,IY,IZ)=BADDATA
                  ADVY_A2(IX,IY,IZ)=BADDATA
                  ADVZ_A2(IX,IY,IZ)=BADDATA
                  SUBSTANTIAL_DDT_A2(IX,IY,IZ)=BADDATA
               ENDIF
3           CONTINUE
2        CONTINUE
1     CONTINUE
      RETURN
      END
