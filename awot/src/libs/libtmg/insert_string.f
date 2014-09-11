      SUBROUTINE INSERT_STRING(STRING,STRING_TO_INSERT,I)

C  Thomas Matejka NOAA/NSSL 16 October 1996

C  This subroutine inserts the character string STRING_TO_INSERT into
C  the character string STRING starting at character number I.  If I
C  occurs before the last of the characters in STRING, the rest of
C  STRING follows the insertion.  If I occurs after all the characters
C  in STRING, the intervening blanks are retained.  All trailing blanks
C  in STRING_TO_INSERT are inserted.  Trailing blanks in STRING that no
C  longer fit are discarded.  If STRING is not long enough, the
C  subroutine stops.

      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      INTEGER,EXTERNAL::S_L
      CHARACTER(LEN=*)::STRING,STRING_TO_INSERT
      INTEGER::I,
     $K,L,M

      K=S_L(STRING)
      L=LEN(STRING)
      M=LEN(STRING_TO_INSERT)
      IF(I.LE.K)THEN
         IF(K+M.GT.L)THEN
            WRITE(7,*)'INSERT_STRING:  STRING IS TOO SHORT.'
            STOP
         ENDIF
         STRING(I+M:K+M)=STRING(I:K)
      ELSEIF(I-1+M.GT.L)THEN
         WRITE(7,*)'INSERT_STRING:  STRING IS TOO SHORT.'
         STOP
      ENDIF
      STRING(I:I-1+M)=STRING_TO_INSERT(1:M)
      END
