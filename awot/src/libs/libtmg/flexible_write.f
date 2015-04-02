      SUBROUTINE FLEXIBLE_WRITE(SUBSTRING,I_START,I_STOP,
     $LU,N_INDENT,TEXT,NEW_LINE)

C  Thomas Matejka NOAA/NSSL 5 April 2002

C  This subroutine writes specified substrings to a specified output
C  device.

C  Input:

C  SUBSTRING (1D character string array 1:I_STOP). SUBSTRING(I)
C  specifies the Ith substring.

C  I_START (integer) specifies that SUBSTRING(I_START) will be the first
C  substring written.

C  I_STOP (integer) specifies that SUBSTRING(I_STOP) will be the last
C  substring written.

C  LU (integer) specifies the unit number of the device to write to.

C  N_INDENT (integer) specifies the number of blanks to write at the
C  beginning of each output record.

C  TEXT (character string) specifies a text to write before the
C  substrings in each output record.

C  NEW_LINE (logical). .FALSE. indicates that the substrings will be
C  written, separated by a blank, in the same output record after
C  N_INDENT blanks and TEXT. .TRUE. indicates that each substring will
C  be written in a separate output record after N_INDENT blanks and
C  TEXT.

      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      INTEGER,EXTERNAL::S_L
      CHARACTER(LEN=MAX_STRING+1)::STRING
      CHARACTER(LEN=*)::TEXT
      CHARACTER(LEN=*),DIMENSION(1:I_STOP)::SUBSTRING
      LOGICAL::NEW_LINE
      INTEGER::I_START,I_STOP,I,LU,N_INDENT

      IF(I_STOP.LT.1.OR.
     $I_STOP.LT.I_START)THEN
         WRITE(TMMLIB_MESSAGE_UNIT,*)'FLEXIBLE_WRITE: ILLEGAL I_STOP.'
         STOP
      ENDIF

      IF(NEW_LINE)THEN
         IF(S_L(TEXT).EQ.0)THEN
            DO I=I_START,I_STOP
               STRING=''
               CALL APPEND_STRING(N_INDENT,SUBSTRING(I),STRING)
               WRITE(LU,"(A)")STRING(1:S_L(STRING))
            ENDDO
         ELSE
            DO I=I_START,I_STOP
               STRING=''
               CALL APPEND_STRING(N_INDENT,TEXT,STRING)
               CALL APPEND_STRING(1,SUBSTRING(I),STRING)
               WRITE(LU,"(A)")STRING(1:S_L(STRING))
            ENDDO
         ENDIF
      ELSE
         IF(S_L(TEXT).EQ.0)THEN
            STRING=''
            DO I=I_START,I_STOP
               IF(I.EQ.I_START)THEN
                  CALL APPEND_STRING(N_INDENT,SUBSTRING(I),STRING)
               ELSE
                  CALL APPEND_STRING(1,SUBSTRING(I),STRING)
               ENDIF
            ENDDO
            WRITE(LU,"(A)")STRING(1:S_L(STRING))
         ELSE
            STRING=''
            CALL APPEND_STRING(N_INDENT,TEXT,STRING)
            DO I=I_START,I_STOP
               CALL APPEND_STRING(1,SUBSTRING(I),STRING)
            ENDDO
            WRITE(LU,"(A)")STRING(1:S_L(STRING))
         ENDIF
      ENDIF

      END SUBROUTINE FLEXIBLE_WRITE
