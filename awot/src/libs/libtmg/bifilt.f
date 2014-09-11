      SUBROUTINE BIFILT(A_IN,ISTART,IEND,N,BADDATA,A_OUT)

C  Thomas Matejka NOAA/NSSL 31 October 1994

C  This subroutine runs an N-point binomial filter on the
C  one-dimensional array A_IN from element ISTART to IEND.  N must be
C  positive and odd.  Values equal to BADDATA will be considered
C  missing.  The filtered vector is returned in the one-dimensional
C  array A_OUT.

      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      INTEGER ISTART,IEND,N,I,J,ITERM,II
      REAL X,Y,SUM1,SUM2,BADDATA
      REAL A_IN(IEND),A_OUT(IEND)
      REAL C(N)

C  Check whether the number of points in the filter is legal.
      IF(N.LT.1.OR.
     $MOD(N,2).EQ.0)THEN
         WRITE(TMMLIB_MESSAGE_UNIT,*)'BIFILT:  N MUST BE POSITIVE AND ',
     $   'ODD.'
         STOP
      ENDIF

C  Calculate the filter weights.
      C(1)=1.
      IF(N.GT.1)THEN
         DO I=2,N
            X=1.
            DO J=N-I+1,N-1
               X=X*FLOAT(J)
            ENDDO
            Y=1.
            DO J=1,I-1
               Y=Y*FLOAT(J)
            ENDDO
            C(I)=X/Y
         ENDDO
      ENDIF

C  Filter the data.
      DO I=ISTART,IEND
         ITERM=I-N/2-1
         SUM1=0.
         SUM2=0.
         DO J=1,N
            II=ITERM+J
            IF(II.GE.ISTART.AND.
     $      II.LE.IEND)THEN
               IF(A_IN(II).NE.BADDATA)THEN
                  SUM1=SUM1+C(J)*A_IN(II)
                  SUM2=SUM2+C(J)
               ENDIF
            ENDIF
         ENDDO
         IF(SUM2.NE.0.)THEN
            A_OUT(I)=SUM1/SUM2
         ELSE
            A_OUT(I)=BADDATA
         ENDIF
      ENDDO
      RETURN
      END
