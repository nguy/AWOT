      SUBROUTINE LINEAR_INTERPOLATION_4D(A,MAXX,MAXY,MAXZ,MAXT,
     $XMIN,YMIN,ZMIN,DELX,DELY,DELZ,NX,NY,NZ,T_DATA,NT,
     $USYS,VSYS,WSYS,
     $BADDATA,
     $X,Y,Z,T,
     $DB,
     $STRICT_CONTROL,
     $EXTEND_X,EXTEND_Y,EXTEND_Z,EXTEND_T,
     $VALUE,
     $IT_FACE,T_FACE,X_DATATIME,Y_DATATIME,Z_DATATIME)

C  Thomas Matejka NOAA/NSSL 11 August 1995

C  This subroutine interpolates a value VALUE from the four-dimensional
C  data field A for the point (X,Y,Z,T).

C  MAXX, MAXY, MAXZ, and MAXT are the first, second, third, and fourth
C  dimensions of A in the calling program.

C  The fourth dimension is interpreted as time.  The interpolation is
C  performed relative to a frame of reference (ideally the most
C  stationary frame of reference of the data field) that is moving at a
C  speed of USYS in the first dimension, VSYS in the second dimension,
C  and WSYS in the third dimension.  If the data field is translating
C  but not evolving, the interpolation will yield identical results to a
C  three-dimensional interpolation of the data field at one time.  If
C  the data field is evolving, the interpolation in time will be of the
C  evolution alone without interference from the translation.

C  The data field A is defined at NT times, that is, at NT grid points
C  in the fourth dimension.  The data field A is defined at NX(L),
C  NY(L), and NZ(L) grid points in the first, second, and third
C  dimensions at the Lth time.  XMIN(L), YMIN(L), and ZMIN(L) are the
C  minimum grid coordinates in the first, second, and third dimensions
C  at the Lth time.  The grid points are spaced DELX(L), DELY(L), and
C  DELZ(L) in the first, second, and third dimensions at the Lth time.
C  The grid may be irregularly spaced in time.  T_DATA(L) specifies the
C  Lth data time in ascending order.

C  When DB is .FALSE., the interpolation is performed directly on the
C  data in A.  When DB is .TRUE., the data in A are assumed to be in
C  decibels, and the interpolation is performed on the linear, not the
C  decibel, values and then converted back to decibels.

C  BADDATA is the value used to indicate missing data in A.  When an
C  interpolated value cannot be obtained, BADDATA is returned.

C  STRICT_CONTROL should be an integer from 1 to 16 that specifies how
C  many of the sixteen surrounding values must be present for the
C  interpolation to proceed.

C  EXTEND_X controls how the interpolation is handled when X lies
C  outside of the domain of A.  If EXTEND_X is .TRUE., then the value of
C  A at the edge of the domain is used.  If EXTEND_X is .FALSE., then
C  BADDATA is returned

C  EXTEND_Y controls how the interpolation is handled when Y lies
C  outside of the domain of A.  If EXTEND_Y is .TRUE., then the value of
C  A at the edge of the domain is used.  If EXTEND_Y is .FALSE., then
C  BADDATA is returned

C  EXTEND_Z controls how the interpolation is handled when Z lies
C  outside of the domain of A.  If EXTEND_Z is .TRUE., then the value of
C  A at the edge of the domain is used.  If EXTEND_Z is .FALSE., then
C  BADDATA is returned

C  EXTEND_T controls how the interpolation is handled when T lies
C  outside of the domain of A.  If EXTEND_T is .TRUE., then the value of
C  A at the edge of the domain is used.  If EXTEND_T is .FALSE., then
C  BADDATA is returned

C  In IT_FACE(1) and IT_FACE(2) are returned the earlier and later time
C  grid point numbers that bracket (X,Y,Z,T).  In T_FACE(1) and
C  T_FACE(2) are returned the earlier and later data times that bracket
C  (X,Y,Z,T).

C  In X_DATATIME(1), Y_DATATIME(1), and Z_DATATIME(1) are returned the
C  first three coordinates of (X,Y,Z,T) projected in the moving frame of
C  reference onto the earlier data time.  In X_DATATIME(2),
C  Y_DATATIME(2), and Z_DATATIME(2) are returned the first three
C  coordinates of (X,Y,Z,T) projected in the moving frame of reference
C  onto the later data time.

      IMPLICIT NONE
      REAL LINEAR_INTERPOLATION
      LOGICAL DB,EXTEND_X,EXTEND_Y,EXTEND_Z,EXTEND_T
      INTEGER MAXX,MAXY,MAXZ,MAXT,NT,I,J,K,L,COUNT,STRICT_CONTROL
      INTEGER IT_FACE(2)
      INTEGER NX(NT),NY(NT),NZ(NT)
      INTEGER IX_FACE(2,2),IY_FACE(2,2),IZ_FACE(2,2)
      REAL USYS,VSYS,WSYS,X,Y,Z,T,BADDATA,VALUE
      REAL CORNER_DATA_1D(2),T_FACE(2),X_DATATIME(2),Y_DATATIME(2),
     $Z_DATATIME(2)
      REAL T_DATA(NT),XMIN(NT),YMIN(NT),ZMIN(NT),DELX(NT),DELY(NT),
     $DELZ(NT)
      REAL CORNER_DATA_2D(2,2),X_FACE(2,2),Y_FACE(2,2),Z_FACE(2,2)
      REAL CORNER_DATA_3D(2,2,2)
      REAL CORNER_DATA_4D(2,2,2,2)
      REAL A(MAXX,MAXY,MAXZ,MAXT)

C  Find the time grid point numbers that bracket T.
      IF(NT.GT.1)THEN
         IF(T.LE.T_DATA(1))THEN
            IT_FACE(1)=1
            IT_FACE(2)=1
         ELSEIF(T.GE.T_DATA(NT))THEN
            IT_FACE(1)=NT
            IT_FACE(2)=NT
         ELSE
            DO L=1,NT-1
               IF(T_DATA(L).EQ.T)THEN
                  IT_FACE(1)=L
                  IT_FACE(2)=L
                  EXIT
               ELSEIF(T_DATA(L).LT.T.AND.
     $         T_DATA(L+1).GT.T)THEN
                  IT_FACE(1)=L
                  IT_FACE(2)=L+1
                  EXIT
               ENDIF
            ENDDO
         ENDIF
      ELSE
         IT_FACE(1)=1
         IT_FACE(2)=1
      ENDIF

C  Find the time coordinates that bracket (X,Y,Z,T).
      DO L=1,2
         T_FACE(L)=T_DATA(IT_FACE(L))
      ENDDO

C  Calculate the coordinates that are obtained when the position
C  (X,Y,Z,T) is projected in the moving frame of reference onto the
C  bracketing data times.
      DO L=1,2
         X_DATATIME(L)=X+USYS*(T_FACE(L)-T)
         Y_DATATIME(L)=Y+VSYS*(T_FACE(L)-T)
         Z_DATATIME(L)=Z+WSYS*(T_FACE(L)-T)
      ENDDO

C  Find the grid point numbers and coordinates that bracket the
C  projected positions at each of the bracketing data times.
      DO L=1,2
         CALL FIND_GRIDBOX_1D(X_DATATIME(L),XMIN(IT_FACE(L)),
     $   DELX(IT_FACE(L)),NX(IT_FACE(L)),IX_FACE(1,L),X_FACE(1,L))
         CALL FIND_GRIDBOX_1D(Y_DATATIME(L),YMIN(IT_FACE(L)),
     $   DELY(IT_FACE(L)),NY(IT_FACE(L)),IY_FACE(1,L),Y_FACE(1,L))
         CALL FIND_GRIDBOX_1D(Z_DATATIME(L),ZMIN(IT_FACE(L)),
     $   DELZ(IT_FACE(L)),NZ(IT_FACE(L)),IZ_FACE(1,L),Z_FACE(1,L))
      ENDDO

C  Find the field values at the sixteen surrounding points.
      DO L=1,2
         DO K=1,2
            DO J=1,2
               DO I=1,2
                  CORNER_DATA_4D(I,J,K,L)=A(IX_FACE(I,L),IY_FACE(J,L),
     $            IZ_FACE(K,L),IT_FACE(L))
               ENDDO
            ENDDO
         ENDDO
      ENDDO

C  Count the number of good values.
      COUNT=0
      DO L=1,2
         DO K=1,2
            DO J=1,2
               DO I=1,2
                  IF(CORNER_DATA_4D(I,J,K,L).NE.BADDATA)THEN
                     COUNT=COUNT+1
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO

C  If there are enough good values, proceed with the interpolation.
      IF(COUNT.GE.STRICT_CONTROL)THEN

C  Interpolate in the first dimension along the eight edges with
C  constant coordinates in the second, third, and fourth dimensions.
         DO L=1,2
            DO K=1,2
               DO J=1,2
                  CORNER_DATA_3D(J,K,L)=
     $            LINEAR_INTERPOLATION(X_FACE(1,L),
     $            CORNER_DATA_4D(1,J,K,L),X_FACE(2,L),
     $            CORNER_DATA_4D(2,J,K,L),X_DATATIME(L),BADDATA,DB,
     $            EXTEND_X)
               ENDDO
            ENDDO
         ENDDO

C  Interpolate in the second dimension along the four resulting edges
C  with constant coordinates in the third and fourth dimensions.
         DO L=1,2
            DO K=1,2
               CORNER_DATA_2D(K,L)=LINEAR_INTERPOLATION(Y_FACE(1,L),
     $         CORNER_DATA_3D(1,K,L),Y_FACE(2,L),CORNER_DATA_3D(2,K,L),
     $         Y_DATATIME(L),BADDATA,DB,EXTEND_Y)
            ENDDO
         ENDDO

C  Interpolate in the third dimension along the two resulting edges with
C  a contant coordinate in the fourth dimension.
         DO L=1,2
            CORNER_DATA_1D(L)=LINEAR_INTERPOLATION(Z_FACE(1,L),
     $      CORNER_DATA_2D(1,L),Z_FACE(2,L),CORNER_DATA_2D(2,L),
     $      Z_DATATIME(L),BADDATA,DB,EXTEND_Z)
         ENDDO

C  Interpolate in the fourth dimension along the one resulting edge.
         VALUE=LINEAR_INTERPOLATION(T_FACE(1),CORNER_DATA_1D(1),
     $   T_FACE(2),CORNER_DATA_1D(2),T,BADDATA,DB,EXTEND_T)
      ELSE

C  There are not enough good values.
         VALUE=BADDATA
      ENDIF

C  Done.
      RETURN
      END
