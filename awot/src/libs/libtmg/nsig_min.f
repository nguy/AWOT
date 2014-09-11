      INTEGER FUNCTION NSIG_MIN(A,N)

C  Thomas Matejka NOAA/NSSL 13 January 1995

C  This function returns the minimum number of significant digits that
C  are required to distinguish the numbers in the one-dimensional array
C  A.

C  Duplicate elements of A do not affect the result.

      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      REAL ROUND_SIGNIF
      LOGICAL AGAIN
      INTEGER N,I,M
      REAL B_PREV,B
      REAL A(N)
      REAL A_SORTED(N)

      DO I=1,N
         A_SORTED(I)=A(I)
      ENDDO
      CALL SORT_F(A_SORTED,N,.TRUE.,M)

      NSIG_MIN=0
      AGAIN=.TRUE.
      DOWHILE(AGAIN)
         AGAIN=.FALSE.
         NSIG_MIN=NSIG_MIN+1
         IF(NSIG_MIN.GT.MAX_SIGNIF)THEN
            WRITE(TMMLIB_MESSAGE_UNIT,*)'NSIG_MIN:  EXCEEDED MAXIMUM ',
     $      'NUMBER OF SIGNIFICANT DIGITS FOR A REAL NUMBER.'
            STOP
         ENDIF
         B_PREV=ROUND_SIGNIF(A_SORTED(1),NSIG_MIN)
         DO I=2,M
            B=ROUND_SIGNIF(A_SORTED(I),NSIG_MIN)
            IF(B_PREV.EQ.B)THEN
               AGAIN=.TRUE.
               EXIT
            ENDIF
            B_PREV=B
         ENDDO
      ENDDO
      RETURN
      END
