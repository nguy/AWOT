      SUBROUTINE FLEXIBLE_READ(LU,LITERAL_CHAR,MAX_SUBSTRINGS,
     $N_SUBSTRINGS,SUBSTRING)

C  Thomas Matejka NOAA/NSSL 2 April 2002

C  This subroutine reads a record from a specified input device and
C  parses the record into substrings. The record is interpreted as a
C  string of characters. Substrings should be separated by a comma
C  and/or one or more blanks. For a substring to contain commas or
C  blanks, it must be preceded and followed by the specified special
C  character. The substrings are returned individually in the elements
C  of an array.

C  Input:

C  LU_READ (integer) specifies the unit number of the device to read
C  from.

C  LITERAL_CHARACTER (character string) specifies the character that
C  precedes and follows a substring that contains literal commas or
C  blanks.

C  MAX_SUBSTRINGS (integer) spcifies maximum number of substrings in the
C  input record.

C  Output:

C  N_SUBSTRINGS (integer) returns the number of substrings in the input
C  record.

C  SUBSTRING (1D character string array 1:MAX_SUBSTRINGS). SUBSTRING(I)
C  returns the Ith substring.

      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      INTEGER,EXTERNAL::S_L
      CHARACTER(LEN=1)::LITERAL_CHAR
      CHARACTER(LEN=MAX_STRING+1)::STRING
      CHARACTER(LEN=*),DIMENSION(MAX_SUBSTRINGS)::SUBSTRING
      LOGICAL::LITERAL
      INTEGER::LU,MAX_SUBSTRINGS,N_SUBSTRINGS,IEND,I,ISTART_SUBSTRING,
     $IEND_SUBSTRING

C  Initialize the output character strings as null.
      DO N_SUBSTRINGS=1,MAX_SUBSTRINGS
         SUBSTRING(N_SUBSTRINGS)=''
      ENDDO

C  Read the string.
      READ(LU,"(A)")STRING
      IEND=S_L(STRING)
      IF(IEND.GT.MAX_STRING)THEN
         WRITE(7,*)'FLEXIBLE_READ:  MEMORY EXCEEDED.  INCREASE ',
     $   'MAX_STRING.'
         STOP
      ENDIF

C  Check that literal strings are properly defined.
      IF(IEND.GE.1)THEN
         LITERAL=.FALSE.
         DO I=1,IEND
            IF(STRING(I:I).EQ.LITERAL_CHAR)THEN
               IF(LITERAL)THEN
                  LITERAL=.FALSE.
               ELSE
                  LITERAL=.TRUE.
               ENDIF
            ENDIF
         ENDDO
         IF(LITERAL)THEN
            WRITE(7,*)'FLEXIBLE_READ:  NO END TO LITERAL STRING.'
            STOP
         ENDIF

C  Make the separations between substrings only commas. First change
C  multiple blanks to single blanks. Then remove blanks adjacent to
C  commas. Then replace blanks with commas.
         I=0
         LITERAL=.FALSE.
         DO
            I=I+1
            IF(STRING(I:I).EQ.LITERAL_CHAR)THEN
               IF(LITERAL)THEN
                  LITERAL=.FALSE.
               ELSE
                  LITERAL=.TRUE.
               ENDIF
            ELSE
               IF(.NOT.LITERAL)THEN
                  IF(I.GE.IEND)THEN
                     EXIT
                  ENDIF
                  IF(STRING(I:I+1).EQ.'  ')THEN
                     STRING(I:IEND-1)=STRING(I+1:IEND)
                     IEND=IEND-1
                     I=I-1
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
         I=0
         LITERAL=.FALSE.
         DO
            I=I+1
            IF(STRING(I:I).EQ.LITERAL_CHAR)THEN
               IF(LITERAL)THEN
                  LITERAL=.FALSE.
               ELSE
                  LITERAL=.TRUE.
               ENDIF
            ELSE
               IF(.NOT.LITERAL)THEN
                  IF(I.GE.IEND)THEN
                     EXIT
                  ENDIF
                  IF(STRING(I:I+1).EQ.' ,')THEN
                     STRING(I:IEND-1)=STRING(I+1:IEND)
                     IEND=IEND-1
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
         I=0
         LITERAL=.FALSE.
         DO
            I=I+1
            IF(STRING(I:I).EQ.LITERAL_CHAR)THEN
               IF(LITERAL)THEN
                  LITERAL=.FALSE.
               ELSE
                  LITERAL=.TRUE.
               ENDIF
            ELSE
               IF(.NOT.LITERAL)THEN
                  IF(I.GE.IEND)THEN
                     EXIT
                  ENDIF
                  IF(STRING(I:I+1).EQ.', ')THEN
                     IF(I.LE.IEND-2)THEN
                        STRING(I+1:IEND-1)=STRING(I+2:IEND)
                     ENDIF
                     IEND=IEND-1
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
         LITERAL=.FALSE.
         DO I=1,IEND
            IF(STRING(I:I).EQ.LITERAL_CHAR)THEN
               IF(LITERAL)THEN
                  LITERAL=.FALSE.
               ELSE
                  LITERAL=.TRUE.
               ENDIF
            ELSE
               IF(.NOT.LITERAL)THEN
                  IF(STRING(I:I).EQ.' ')THEN
                     STRING(I:I)=','
                  ENDIF
               ENDIF
            ENDIF
         ENDDO

C  Parse the string into substrings.
         N_SUBSTRINGS=0
         ISTART_SUBSTRING=1
         IEND_SUBSTRING=0
         LITERAL=.FALSE.
         DO I=1,IEND
            IF(STRING(I:I).EQ.LITERAL_CHAR)THEN
               IF(LITERAL)THEN
                  IEND_SUBSTRING=I
                  LITERAL=.FALSE.
               ELSE
                  ISTART_SUBSTRING=I
                  LITERAL=.TRUE.
               ENDIF
            ENDIF
            IF(.NOT.LITERAL)THEN
               IF(STRING(I:I).EQ.',')THEN
                  N_SUBSTRINGS=N_SUBSTRINGS+1
                  IF(N_SUBSTRINGS.GT.MAX_SUBSTRINGS)THEN
                     WRITE(7,*)'FLEXIBLE_READ:  MAX_SUBSTRINGS IS TOO ',
     $               'SMALL.'
                     STOP
                  ENDIF
                  IF(STRING(ISTART_SUBSTRING:ISTART_SUBSTRING).EQ.
     $            LITERAL_CHAR)THEN
                     IF(IEND_SUBSTRING-1.GE.ISTART_SUBSTRING+1)THEN
                        SUBSTRING(N_SUBSTRINGS)=
     $                  STRING(ISTART_SUBSTRING+1:IEND_SUBSTRING-1)
                     ENDIF
                  ELSE
                     IF(IEND_SUBSTRING.GE.ISTART_SUBSTRING)THEN
                        SUBSTRING(N_SUBSTRINGS)=
     $                  STRING(ISTART_SUBSTRING:IEND_SUBSTRING)
                     ENDIF
                  ENDIF
                  IF(I.EQ.IEND)THEN
                     N_SUBSTRINGS=N_SUBSTRINGS+1
                     IF(N_SUBSTRINGS.GT.MAX_SUBSTRINGS)THEN
                        WRITE(7,*)'FLEXIBLE_READ:  MAX_SUBSTRINGS IS ',
     $                 'TOO SMALL.'
                       STOP
                     ENDIF
                  ELSE
                     ISTART_SUBSTRING=I+1
                     IEND_SUBSTRING=I
                  ENDIF
               ELSE
                  IEND_SUBSTRING=I
                  IF(I.EQ.IEND)THEN
                     N_SUBSTRINGS=N_SUBSTRINGS+1
                     IF(N_SUBSTRINGS.GT.MAX_SUBSTRINGS)THEN
                        WRITE(7,*)'FLEXIBLE_READ:  MAX_SUBSTRINGS IS ',
     $                  'TOO SMALL.'
                        STOP
                     ENDIF
                     IF(STRING(ISTART_SUBSTRING:ISTART_SUBSTRING).EQ.
     $               LITERAL_CHAR)THEN
                        IF(IEND_SUBSTRING-1.GE.ISTART_SUBSTRING+1)THEN
                           SUBSTRING(N_SUBSTRINGS)=
     $                     STRING(ISTART_SUBSTRING+1:IEND_SUBSTRING-1)
                        ENDIF
                     ELSE
                        IF(IEND_SUBSTRING.GE.ISTART_SUBSTRING)THEN
                           SUBSTRING(N_SUBSTRINGS)=
     $                     STRING(ISTART_SUBSTRING:IEND_SUBSTRING)
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
      ELSE
         N_SUBSTRINGS=0
      ENDIF

      END SUBROUTINE FLEXIBLE_READ
