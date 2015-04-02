      REAL FUNCTION CUBIC_SPLINE_INTERP(N,X,Q,X0,
     $BADDATA,
     $EXTRAP_MODE)

C  Thomas Matejka NOAA/NSSL 29 April 1998
C=======================================================================
C  This real function performs natural cubic spline interpolation on a
C  set of discrete data in one dimension and returns the interpolated
C  value.  If an interpolated value cannot be calculated, it returns
C  BADDATA.

C  Input:

C  N (integer) is the number of data to serve as knots for the cubic
C  spline.

C  X (1d real array 1:N).  X(I) specifies the coordinate of the Ith
C  datum.  X(I) must monotonically increase for I=1,...,N.

C  Q (1d real array 1:N).  Q(I) specifies the value of the cubic spline
C  at X(I).  It should never be missing.

C  X0 (real) specifies the coordinate to interpolate to.

C  BADDATA (real) indicates missing values as described.

C  EXTRAP_MODE (character string) controls what the function returns
C  when X lies outside of [X(1),X(N)] or [X(N),X(1)].  If it is EXTEND,
C  the nearest data value is returned.  If it is NOTHING, BADDATA is
C  returned.
C=======================================================================
      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      CHARACTER(LEN=*)::EXTRAP_MODE
      LOGICAL::SUCCESS
      INTEGER::N,I
      REAL::X0,DX,BADDATA
      REAL,DIMENSION(1:N)::X,Q,GPP

      IF(EXTRAP_MODE.NE.'EXTEND'.AND.
     $EXTRAP_MODE.NE.'NOTHING')THEN
         WRITE(TMMLIB_MESSAGE_UNIT,*)'CUBIC_SPLINE_INTERP:  ILLEGAL ',
     $   'EXTRAP_MODE.'
         STOP
      ENDIF
      IF(X0.LT.X(1))THEN
         IF(EXTRAP_MODE.EQ.'EXTEND')THEN
            CUBIC_SPLINE_INTERP=Q(1)
         ELSE
            CUBIC_SPLINE_INTERP=BADDATA
         ENDIF
      ELSEIF(X0.GT.X(N))THEN
         IF(EXTRAP_MODE.EQ.'EXTEND')THEN
            CUBIC_SPLINE_INTERP=Q(N)
         ELSE
            CUBIC_SPLINE_INTERP=BADDATA
         ENDIF
      ELSEIF(N.EQ.1)THEN
         IF(X0.EQ.X(1))THEN
            CUBIC_SPLINE_INTERP=Q(1)  
         ELSE
            CUBIC_SPLINE_INTERP=BADDATA
         ENDIF
      ELSE
         CALL CUBIC_SPLINE(N,X,Q,GPP,SUCCESS)
         IF(.NOT.SUCCESS)THEN
            CUBIC_SPLINE_INTERP=BADDATA
            RETURN
         ENDIF
         DO I=2,N
            IF(X0.LE.X(I))THEN
               DX=X(I)-X(I-1)
               CUBIC_SPLINE_INTERP=
     $         GPP(I-1)*((X(I)-X0)**3/DX-DX*(X(I)-X0))/6.+
     $         GPP(I)*((X0-X(I-1))**3/DX-DX*(X0-X(I-1)))/6.+
     $         Q(I-1)*(X(I)-X0)/DX+
     $         Q(I)*(X0-X(I-1))/DX
               EXIT
            ENDIF
         ENDDO
      ENDIF

      END FUNCTION CUBIC_SPLINE_INTERP
