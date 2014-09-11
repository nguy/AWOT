      SUBROUTINE D2DT2_3D(A1,A2,A3,MAXX_IN,MAXY_IN,MAXZ_IN,NX,NY,NZ,T1,
     $T2,T3,BADDATA,MAXX_OUT,MAXY_OUT,MAXZ_OUT,D2DT2_A2)

C  Thomas Matejka NOAA/NSSL 22 August 1994

C  This subroutine calculates the partial second derivative with respect
C  to a fourth dimension of a three-dimensional data field input at
C  three positions in the fourth dimension.  The data spacing in the
C  fourth dimension does not have to be uniform.  A three-point
C  second-order linear regression is used.  Note that, if the first
C  three dimensions are spatial and the fourth dimension is time, the
C  results will be the second time derivative with respect to a
C  stationary frame of reference.  To get the second time derivative
C  with respect to a moving frame of reference (for example, the most
C  stationary frame of reference), data fields A1 and A3 should first be
C  shifted to time T2 at the velocity of the frame of reference.

C  Input:

C  A1 is a three-dimensional real array.  A1(I,J,K) specifies the value,
C  at the Ith grid point in the first dimension, the Jth grid point in
C  the second dimension, the Kth grid point in the third dimension, and
C  the first grid point in the fourth dimension of the data field whose
C  second derivative is sought.  If it is missing, it should equal
C  BADDATA.

C  A2 is a three-dimensional real array.  A2(I,J,K) specifies the value,
C  at the Ith grid point in the first dimension, the Jth grid point in
C  the second dimension, the Kth grid point in the third dimension, and
C  the second grid point in the fourth dimension of the data field whose
C  second derivative is sought.  If it is missing, it should equal
C  BADDATA.

C  A3 is a three-dimensional real array.  A3(I,J,K) specifies the value,
C  at the Ith grid point in the first dimension, the Jth grid point in
C  the second dimension, the Kth grid point in the third dimension, and
C  the third grid point in the fourth dimension of the data field whose
C  second derivative is sought.  If it is missing, it should equal
C  BADDATA.

C  MAXX_IN is an integer variable that specifies the first dimension of
C  A1, A2, and A3 in the calling program.

C  MAXY_IN is an integer variable that specifies the second dimension of
C  A1, A2, and A3 in the calling program.

C  MAXZ_IN is an integer variable that specifies the third dimension of
C  A1, A2, and A3 in the calling program.

C  NX is an integer variable that specifies the number of grid points in
C  the first dimension for A1, A2, A3, and D2DT2_A2.

C  NY is an integer variable that specifies the number of grid points in
C  the second dimension for A1, A2, A3, and D2DT2_A2.

C  NZ is an integer variable that specifies the number of grid points in
C  the third dimension for A1, A2, A3, and D2DT2_A2.

C  T1 is a real variable that specifies the value of the coordinate in
C  the fourth dimension for A1.

C  T2 is a real variable that specifies the value of the coordinate in
C  the fourth dimension for A2 and D2DT2_A2.

C  T3 is a real variable that specifies the value of the coordinate in
C  the fourth dimension for A3.

C  BADDATA is a real variable that specifies the value used to indicate
C  a missing value as described.

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
C  dimension, and the second grid point in the fourth dimension, of the
C  partial second derivative of A2 with respect to the fourth dimension.
C  If it is missing, it is returned as BADDATA.

      IMPLICIT NONE
      INTEGER MAXX_IN,MAXY_IN,MAXZ_IN,MAXX_OUT,MAXY_OUT,MAXZ_OUT,NX,NY,
     $NZ,IX,IY,IZ
      REAL BADDATA,T1,T2,T3,F1,F2,F3
      REAL A1(MAXX_IN,MAXY_IN,MAXZ_IN),A2(MAXX_IN,MAXY_IN,MAXZ_IN),
     $A3(MAXX_IN,MAXY_IN,MAXZ_IN)
      REAL D2DT2_A2(MAXX_OUT,MAXY_OUT,MAXZ_OUT)

      DO 1 IZ=1,NZ
         DO 2 IY=1,NY
            DO 3 IX=1,NX
               F1=A1(IX,IY,IZ)
               F2=A2(IX,IY,IZ)
               F3=A3(IX,IY,IZ)
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
