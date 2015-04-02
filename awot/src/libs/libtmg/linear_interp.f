      REAL FUNCTION LINEAR_INTERP(X1,Q1,
     $X2,Q2,
     $X,
     $BADDATA,
     $EXTRAP_MODE,INTERP_ONE_GOOD)

C  Thomas Matejka NOAA/NSSL 23 June 1998
C=======================================================================
C  This real function performs linear interpolation on two discrete data
C  in one dimension and returns the linearly interpolated value.  If an
C  interpolated value cannot be calculated, it returns BADDATA.

C  Input:

C  X1 (real) specifies a coordinate value.

C  Q1 (real) specifies the data value at X1.  If it is missing, it
C  should equal BADDATA.

C  X2 (real) specifies a coordinate value.

C  Q2 (real) specifies the data value at X2.  If it is missing, it
C  should equal BADDATA.

C  X (real) specifies the coordinate value to interpolate to.

C  BADDATA (real) indicates missing values as described.

C  EXTRAP_MODE (character string) controls what the function returns
C  when X lies outside of [X1,X2] or [X2,X1].  If it is LINEAR,
C  extrapolation is performed with the linear trend used for
C  interpolation.  If it is EXTEND, the nearer data value is returned.
C  If it is NOTHING, BADDATA is returned.

C  INTERP_ONE_GOOD (logical) controls what the function returns when X
C  lies inside of [X1,X2] or [X2,X1] but when only one of Q1 and Q2
C  exists.  If it is .TRUE., the one good data value is returned.  If it
C  is .FALSE., BADDATA is returned.
C=======================================================================
      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      CHARACTER(LEN=*)::EXTRAP_MODE
      LOGICAL::INTERP_ONE_GOOD
      REAL::X1,X2,X,Q1,Q2,BADDATA

      IF(EXTRAP_MODE.NE.'LINEAR'.AND.
     $EXTRAP_MODE.NE.'EXTEND'.AND.
     $EXTRAP_MODE.NE.'NOTHING')THEN
         WRITE(TMMLIB_MESSAGE_UNIT,*)'LINEAR_INTERP:  ILLEGAL ',
     $   'EXTRAP_MODE.'
         STOP
      ENDIF
      IF(X1.EQ.X2)THEN
         IF(Q1.EQ.Q2)THEN
            IF(X.EQ.X1.OR.
     $      EXTRAP_MODE.EQ.'LINEAR'.OR.
     $      EXTRAP_MODE.EQ.'EXTEND')THEN
               LINEAR_INTERP=Q1
            ELSE
               LINEAR_INTERP=BADDATA
            ENDIF
         ELSE
            LINEAR_INTERP=BADDATA
         ENDIF
      ELSE
         IF(Q1.NE.BADDATA.AND.
     $   Q2.NE.BADDATA)THEN
            IF((X.GE.X1.AND.X.LE.X2).OR.
     $      (X.LE.X1.AND.X.GE.X2).OR.
     $      EXTRAP_MODE.EQ.'LINEAR')THEN
               LINEAR_INTERP=Q1+(X-X1)*(Q2-Q1)/(X2-X1)
            ELSEIF(EXTRAP_MODE.EQ.'EXTEND')THEN
               IF(ABS(X-X1).LT.ABS(X-X2))THEN
                  LINEAR_INTERP=Q1
               ELSE
                  LINEAR_INTERP=Q2
               ENDIF
            ELSE
               LINEAR_INTERP=BADDATA
            ENDIF
         ELSEIF(Q1.NE.BADDATA)THEN
            IF((X.GE.X1.AND.X.LE.X2).OR.
     $      (X.LE.X1.AND.X.GE.X2))THEN
               IF(INTERP_ONE_GOOD)THEN
                  LINEAR_INTERP=Q1
               ELSE
                  LINEAR_INTERP=BADDATA
               ENDIF
            ELSEIF(EXTRAP_MODE.EQ.'EXTEND')THEN
               IF(ABS(X-X1).LT.ABS(X-X2))THEN
                  LINEAR_INTERP=Q1
               ELSE
                  LINEAR_INTERP=BADDATA
               ENDIF
            ELSE
               LINEAR_INTERP=BADDATA
            ENDIF
         ELSEIF(Q2.NE.BADDATA)THEN
            IF((X.GE.X1.AND.X.LE.X2).OR.
     $      (X.LE.X1.AND.X.GE.X2))THEN
               IF(INTERP_ONE_GOOD)THEN
                  LINEAR_INTERP=Q2
               ELSE
                  LINEAR_INTERP=BADDATA
               ENDIF
            ELSEIF(EXTRAP_MODE.EQ.'EXTEND')THEN
               IF(ABS(X-X1).LT.ABS(X-X2))THEN
                  LINEAR_INTERP=BADDATA
               ELSE
                  LINEAR_INTERP=Q2
               ENDIF
            ELSE
               LINEAR_INTERP=BADDATA
            ENDIF
         ELSE
            LINEAR_INTERP=BADDATA
         ENDIF            
      ENDIF

      END FUNCTION LINEAR_INTERP
