      SUBROUTINE D2DT2_3D_MOVING(A1,A2,A3,MAXX_IN,MAXY_IN,MAXZ_IN,XMIN1,
     $YMIN1,ZMIN1,DELX1,DELY1,DELZ1,NX1,NY1,NZ1,XMIN2,YMIN2,ZMIN2,DELX2,
     $DELY2,DELZ2,NX2,NY2,NZ2,XMIN3,YMIN3,ZMIN3,DELX3,DELY3,DELZ3,NX3,
     $NY3,NZ3,T1,T2,T3,BADDATA,USYS,VSYS,WSYS,MAXX_OUT,MAXY_OUT,
     $MAXZ_OUT,D2DT2_A2)

C  Thomas Matejka NOAA/NSSL 22 August 1994

C  This subroutine calculates the partial second derivative with respect
C  to time of a three-dimensional data field input at three times.  The
C  data spacing in time does not have to be uniform.  A three-point
C  second-order linear regression is used.  Time derivatives may be
C  calculated with respect to a moving frame of reference.

C  Input:

C  A1 is a three-dimensional real array.  A1(I,J,K) specifies the value,
C  at the Ith grid point in the first dimension, the Jth grid point in
C  the second dimension, the Kth grid point in the third dimension, and
C  the first time of the data field whose second derivative is sought.
C  If it is missing, it should equal BADDATA.

C  A2 is a three-dimensional real array.  A2(I,J,K) specifies the value,
C  at the Ith grid point in the first dimension, the Jth grid point in
C  the second dimension, the Kth grid point in the third dimension, and
C  the second time of the data field whose second derivative is sought.
C  If it is missing, it should equal BADDATA.

C  A3 is a three-dimensional real array.  A3(I,J,K) specifies the value,
C  at the Ith grid point in the first dimension, the Jth grid point in
C  the second dimension, the Kth grid point in the third dimension, and
C  the third time of the data field whose second derivative is sought.
C  If it is missing, it should equal BADDATA.

C  MAXX_IN is an integer variable that specifies the first dimension of
C  A1, A2, and A3 in the calling program.

C  MAXY_IN is an integer variable that specifies the second dimension of
C  A1, A2, and A3 in the calling program.

C  MAXZ_IN is an integer variable that specifies the third dimension of
C  A1, A2, and A3 in the calling program.

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
C  first dimension for A2 and D2DT2_A2.

C  YMIN2 is a real variable that specifies the minimum coordinate in the
C  second dimension for A2 and D2DT2_A2.

C  ZMIN2 is a real variable that specifies the minimum coordinate in the
C  third dimension for A2 and D2DT2_A2.

C  DELX2 is a real variable that specifies the grid spacing in the first
C  dimension for A2 and D2DT2_A2.

C  DELY2 is a real variable that specifies the grid spacing in the
C  second dimension for A2 and D2DT2_A2.

C  DELZ2 is a real variable that specifies the grid spacing in the third
C  dimension for A2 and D2DT2_A2.

C  NX2 is an integer variable that specifies the number of grid points
C  in the first dimension for A2 and D2DT2_A2.

C  NY2 is an integer variable that specifies the number of grid points
C  in the second dimension for A2 and D2DT2_A2.

C  NZ2 is an integer variable that specifies the number of grid points
C  in the third dimension for A2 and D2DT2_A2.

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

C  T2 is a real variable that specifies the time for A2 and D2DT2_A2.

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
C  D2DT2_A2 in the calling program.

C  MAXY_OUT is an integer variable that specifies the second dimension
C  of D2DT2_A2 in the calling program.

C  MAXZ_OUT is an integer variable that specifies the third dimension of
C  D2DT2_A2 in the calling program.

C  Output:

C  D2DT2_A2 is a three-dimensional real array.  D2DT2_A2(I,J,K) returns
C  the value, at the Ith grid point in the first dimension, the Jth grid
C  point in the second dimension, the Kth grid point in the third
C  dimension, and the second time, of the partial second derivative of
C  A2 with respect to time in the moving frame of reference.  If it is
C  missing, it is returned as BADDATA.

      IMPLICIT NONE
      REAL LINEAR_INTERPOLATION_3D
      INTEGER MAXX_IN,MAXY_IN,MAXZ_IN,MAXX_OUT,MAXY_OUT,MAXZ_OUT,NX1,
     $NY1,NZ1,NX2,NY2,NZ2,NX3,NY3,NZ3,IX,IY,IZ
      REAL BADDATA,T1,T2,T3,F1,F2,F3,USYS,VSYS,WSYS,XMIN1,YMIN1,ZMIN1,
     $DELX1,DELY1,DELZ1,XMIN2,YMIN2,ZMIN2,DELX2,DELY2,DELZ2,XMIN3,YMIN3,
     $ZMIN3,DELX3,DELY3,DELZ3,X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3
      REAL A1(MAXX_IN,MAXY_IN,MAXZ_IN),A2(MAXX_IN,MAXY_IN,MAXZ_IN),
     $A3(MAXX_IN,MAXY_IN,MAXZ_IN)
      REAL D2DT2_A2(MAXX_OUT,MAXY_OUT,MAXZ_OUT)

      DO 1 IZ=1,NZ2
         Z2=ZMIN2+FLOAT(IZ-1)*DELZ2
         Z1=Z2-WSYS*(T2-T1)
         Z3=Z2-WSYS*(T2-T3)
         DO 2 IY=1,NY2
            Y2=YMIN2+FLOAT(IY-1)*DELY2
            Y1=Y2-VSYS*(T2-T1)
            Y3=Y2-VSYS*(T2-T3)
            DO 3 IX=1,NX2
               X2=XMIN2+FLOAT(IX-1)*DELX2
               X1=X2-USYS*(T2-T1)
               X3=X2-USYS*(T2-T3)
               F1=LINEAR_INTERPOLATION_3D(A1,MAXX_IN,MAXY_IN,MAXZ_IN,
     $         XMIN1,YMIN1,ZMIN1,DELX1,DELY1,DELZ1,NX1,NY1,NZ1,BADDATA,
     $         X1,Y1,Z1,.FALSE.,8,.FALSE.,.FALSE.,.FALSE.)
               F2=A2(IX,IY,IZ)
               F3=LINEAR_INTERPOLATION_3D(A3,MAXX_IN,MAXY_IN,MAXZ_IN,
     $         XMIN3,YMIN3,ZMIN3,DELX3,DELY3,DELZ3,NX3,NY3,NZ3,BADDATA,
     $         X3,Y3,Z3,.FALSE.,8,.FALSE.,.FALSE.,.FALSE.)
               IF(F1.NE.BADDATA.AND.
     $         F3.NE.BADDATA.AND.
     $         F2.NE.BADDATA)THEN
                  D2DT2_A2(IX,IY,IZ)=2.*(F1*(T3-T2)+F2*(T1-T3)+
     $            F3*(T2-T1))/((T3-T2)*(T2-T1)*(T3-T1))
               ELSE
                  D2DT2_A2(IX,IY,IZ)=BADDATA
               ENDIF
3           CONTINUE
2        CONTINUE
1     CONTINUE
      RETURN
      END
