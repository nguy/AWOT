      REAL FUNCTION LINEAR_INTERPOLATION(X1,Q1,X2,Q2,X,BADDATA,DB,
     $EXTEND)

C  Thomas Matejka NOAA/NSSL 26 October 1993

C  This function returns a linearly interpolated value at X, given
C  values Q1 and Q2 at X1 and X2.  X may lie inside or outside of
C  [X1,X2].

C  When DB is .FALSE., the interpolation is performed directly on Q1 and
C  Q2.  When DB is .TRUE., Q1 and Q2 are assumed to be in decibels, and
C  the interpolation is performed on the linear, not the decibel, values
C  and then converted back to decibels.

C  When Q1 or Q2 equals BADDATA, it is treated as missing.  Q1 is
C  returned if Q2 is missing, and Q2 is returned if Q1 is missing.

C  X1 equalling X2 and not equalling X is treated as a signal that X
C  lies outside of the data domain whose edge is at X1.  EXTEND controls
C  how the interpolation is handled if X1 equals X2, Q1 equals Q2, and X
C  does not equal X1.  If, in this situation, EXTEND is .FALSE., then
C  BADDATA is returned.  If EXTEND is .TRUE., then Q1 is returned,
C  thereby extrapolating the edge value outside of the domain.  If X1
C  equals X2 and Q1 does not equal Q2, then BADDATA is always returned.

      IMPLICIT NONE
      REAL X_TO_DBX,DBX_TO_X
      LOGICAL DB,EXTEND
      REAL X1,X2,X,Q1,Q2,Q1_LINEAR,BADDATA

      IF(Q1.NE.BADDATA.AND.Q2.NE.BADDATA)THEN
         IF(X1.NE.X2)THEN
            IF(.NOT.DB)THEN
               LINEAR_INTERPOLATION=Q1+(X-X1)*(Q2-Q1)/(X2-X1)
            ELSE
               Q1_LINEAR=DBX_TO_X(Q1)
               LINEAR_INTERPOLATION=X_TO_DBX(Q1_LINEAR+
     $         (X-X1)*(DBX_TO_X(Q2)-Q1_LINEAR)/(X2-X1))
            ENDIF
         ELSE
            IF(Q1.EQ.Q2)THEN
               IF(EXTEND)THEN
                  LINEAR_INTERPOLATION=Q1
               ELSE
                  IF(X.EQ.X1)THEN
                     LINEAR_INTERPOLATION=Q1
                  ELSE
                     LINEAR_INTERPOLATION=BADDATA
                  ENDIF
               ENDIF
            ELSE
               LINEAR_INTERPOLATION=BADDATA
            ENDIF
         ENDIF
      ELSE
         IF(X1.NE.X2)THEN
            IF(Q1.NE.BADDATA)THEN
               LINEAR_INTERPOLATION=Q1
            ELSE
               LINEAR_INTERPOLATION=Q2
            ENDIF
         ELSE
            LINEAR_INTERPOLATION=BADDATA
         ENDIF
      ENDIF
      RETURN
      END
