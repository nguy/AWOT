      SUBROUTINE DDT_3D(A1,A2,A3,MAXX_IN,MAXY_IN,MAXZ_IN,NX,NY,NZ,T1,T2,
     $T3,BADDATA,MAXX_OUT,MAXY_OUT,MAXZ_OUT,DDT_A2)

C  Thomas Matejka NOAA/NSSL 22 August 1995

C  This subroutine calculates the partial derivative with respect to a
C  fourth dimension of a three-dimensional data field input at three
C  positions in the fourth dimension.  The data spacing in the fourth
C  dimension does not have to be uniform.  A three-point first-order
C  linear regression is used where possible.  Otherwise, a two-point
C  difference over two data intervals or one data interval in the fourth
C  dimension is used where possible.  Note that, if the first three
C  dimensions are spatial and the fourth dimension is time, the results
C  will be the time derivative with respect to a stationary frame of
C  reference.  To get the time derivative with respect to a moving frame
C  of reference (for example, the most stationary frame of reference),
C  data fields A1 and A3 should first be shifted to time T2 at the
C  velocity of the frame of reference.

C  Input:

C  A1 is a three-dimensional real array.  A1(I,J,K) specifies the value,
C  at the Ith grid point in the first dimension, the Jth grid point in
C  the second dimension, the Kth grid point in the third dimension, and
C  the first grid point in the fourth dimension of the data field whose
C  derivative is sought.  If it is missing, it should equal BADDATA.

C  A2 is a three-dimensional real array.  A2(I,J,K) specifies the value,
C  at the Ith grid point in the first dimension, the Jth grid point in
C  the second dimension, the Kth grid point in the third dimension, and
C  the second grid point in the fourth dimension of the data field whose
C  derivative is sought.  If it is missing, it should equal BADDATA.

C  A3 is a three-dimensional real array.  A3(I,J,K) specifies the value,
C  at the Ith grid point in the first dimension, the Jth grid point in
C  the second dimension, the Kth grid point in the third dimension, and
C  the third grid point in the fourth dimension of the data field whose
C  derivative is sought.  If it is missing, it should equal BADDATA.

C  MAXX_IN is an integer variable that specifies the first dimension of
C  A1, A2, and A3 in the calling program.

C  MAXY_IN is an integer variable that specifies the second dimension of
C  A1, A2, and A3 in the calling program.

C  MAXZ_IN is an integer variable that specifies the third dimension of
C  A1, A2, and A3 in the calling program.

C  NX is an integer variable that specifies the number of grid points in
C  the first dimension for A1, A2, A3, and DDT_A2.

C  NY is an integer variable that specifies the number of grid points in
C  the second dimension for A1, A2, A3, and DDT_A2.

C  NZ is an integer variable that specifies the number of grid points in
C  the third dimension for A1, A2, A3, and DDT_A2.

C  T1 is a real variable that specifies the value of the coordinate in
C  the fourth dimension for A1.

C  T2 is a real variable that specifies the value of the coordinate in
C  the fourth dimension for A2 and DDT_A2.

C  T3 is a real variable that specifies the value of the coordinate in
C  the fourth dimension for A3.

C  BADDATA is a real variable that specifies the value used to indicate
C  a missing value as described.

C  MAXX_OUT is an integer variable that specifies the first dimension of
C  DDT_A2 in the calling program.

C  MAXY_OUT is an integer variable that specifies the second dimension
C  of DDT_A2 in the calling program.

C  MAXZ_OUT is an integer variable that specifies the third dimension of
C  DDT_A2 in the calling program.

C  Output:

C  DDT_A2 is a three-dimensional real array.  DDT_A2(I,J,K) returns the
C  value, at the Ith grid point in the first dimension, the Jth grid
C  point in the second dimension, the Kth grid point in the third
C  dimension, and the second grid point in the fourth dimension, of the
C  partial derivative of A2 with respect to the fourth dimension.  If it
C  is missing, it is returned as BADDATA.

      IMPLICIT NONE
      INTEGER MAXX_IN,MAXY_IN,MAXZ_IN,MAXX_OUT,MAXY_OUT,MAXZ_OUT,NX,NY,
     $NZ,IX,IY,IZ
      REAL BADDATA,T1,T2,T3,F1,F2,F3
      REAL A1(MAXX_IN,MAXY_IN,MAXZ_IN),A2(MAXX_IN,MAXY_IN,MAXZ_IN),
     $A3(MAXX_IN,MAXY_IN,MAXZ_IN)
      REAL DDT_A2(MAXX_OUT,MAXY_OUT,MAXZ_OUT)

      DO 1 IZ=1,NZ
         DO 2 IY=1,NY
            DO 3 IX=1,NX
               F1=A1(IX,IY,IZ)
               F2=A2(IX,IY,IZ)
               F3=A3(IX,IY,IZ)
               IF(F1.NE.BADDATA)THEN
                  IF(F3.NE.BADDATA)THEN
                     IF(F2.NE.BADDATA)THEN
                        DDT_A2(IX,IY,IZ)=(3.*(T1*F1+T2*F2+T3*F3)-
     $                  (T1+T2+T3)*(F1+F2+F3))/(3.*(T1**2+T2**2+T3**2)-
     $                  (T1+T2+T3)**2)
                     ELSE
                        DDT_A2(IX,IY,IZ)=(F3-F1)/(T3-T1)
                     ENDIF
                  ELSEIF(F2.NE.BADDATA)THEN
                     DDT_A2(IX,IY,IZ)=(F2-F1)/(T2-T1)
                  ELSE
                     DDT_A2(IX,IY,IZ)=BADDATA
                  ENDIF
               ELSEIF(F2.NE.BADDATA.AND.
     $         F3.NE.BADDATA)THEN
                  DDT_A2(IX,IY,IZ)=(F3-F2)/(T3-T2)
               ELSE
                  DDT_A2(IX,IY,IZ)=BADDATA
               ENDIF
3           CONTINUE
2        CONTINUE
1     CONTINUE
      RETURN
      END
