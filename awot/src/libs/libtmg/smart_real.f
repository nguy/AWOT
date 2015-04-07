      SUBROUTINE SMART_REAL(A,E,STRING_OUT)

C  Thomas Matejka NOAA/NSSL 22 May 2000

      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      INTEGER,EXTERNAL::S_L
      CHARACTER(LEN=MAX_STRING)::STRING
      CHARACTER(LEN=*)::STRING_OUT
      LOGICAL::SCIENTIFIC
      INTEGER::E,I,IE,L,I1,J,Y,IEND,IDP,K
      REAL::A,B

C  Initialize.
      STRING=''
      STRING_OUT=''
      L=LEN(STRING_OUT)

C  Zero.
      IF(A.EQ.0.)THEN
         STRING_OUT='0'

C  Nonzero.
      ELSE
         B=ABS(A)
         WRITE(STRING,*)B
         IF(S_L(STRING).GE.MAX_STRING)THEN
            WRITE(TMMLIB_MESSAGE_UNIT,*)'SMART_REAL: MEMORY EXCEEDED. ',
     $      'INCREASE MAX_STRING.'
            STOP
         ENDIF
         CALL LEFT_JUSTIFY(STRING)
         IEND=S_L(STRING)
         DO I=1,IEND
            IF(STRING(I:I).EQ.'.')THEN
               IDP=I
               EXIT
            ENDIF
         ENDDO
         DO I=1,IEND
            IF(STRING(I:I).NE.'0'.AND.
     $      STRING(I:I).NE.'.')THEN
               I1=I
               EXIT
            ENDIF
         ENDDO
         CALL FINDCHAR(STRING,'E',1,IE,SCIENTIFIC)
         IF(.NOT.SCIENTIFIC)THEN
            CALL FINDCHAR(STRING,'e',1,IE,SCIENTIFIC)
         ENDIF

C  The unformatted write produced scientific notation. Determine its
C  standard scientific notation.
         IF(SCIENTIFIC)THEN
            READ(STRING(IE+1:IEND),*)Y
            IF(IDP.GT.I1)THEN
               Y=Y+IDP-I1
            ELSE
               Y=Y+IDP-I1+1
            ENDIF
            J=0
            DO I=I1,IE-1
               IF(STRING(I:I).NE.'.')THEN
                  J=J+1
                  STRING(J:J)=STRING(I:I)
               ENDIF
            ENDDO
            IEND=J

C  The unformatted write produced decimal notation. Determine its
C  standard scientific notation.
         ELSE
            IF(IDP.GT.I1)THEN
               Y=IDP-I1
            ELSE
               Y=IDP-I1+1
            ENDIF
            J=0
            DO I=I1,IEND
               IF(STRING(I:I).NE.'.')THEN
                  J=J+1
                  STRING(J:J)=STRING(I:I)
               ENDIF
            ENDDO
            IEND=J
         ENDIF
         
C  Truncate zeroes.
         DO I=IEND,1,-1
            IF(STRING(I:I).NE.'0')THEN
               IEND=I
               STRING(IEND+1:MAX_STRING)=''
               EXIT
            ENDIF
         ENDDO

C  Negative number.
         J=0
         IF(A.LT.0.)THEN
            J=J+1
            IF(J.GT.L)THEN
               WRITE(TMMLIB_MESSAGE_UNIT,*)'SMART_REAL: STRING_OUT IS ',
     $         'TOO SHORT.'
               STOP
            ENDIF
            STRING_OUT(J:J)='-'
         ENDIF

C  Large or small nubmer.
         IF(B.GE.10.**E.OR.
     $   B.LT.10.**(-E))THEN
            J=J+1
            IF(J.GT.L)THEN
               WRITE(TMMLIB_MESSAGE_UNIT,*)'SMART_REAL: STRING_OUT IS ',
     $         'TOO SHORT.'
               STOP
            ENDIF
            STRING_OUT(J:J)='.'
            J=J+1
            IF(J+IEND-1.GT.L)THEN
               WRITE(TMMLIB_MESSAGE_UNIT,*)'SMART_REAL: STRING_OUT IS ',
     $         'TOO SHORT.'
               STOP
            ENDIF
            STRING_OUT(J:J+IEND-1)=STRING(1:IEND)
            J=J+IEND
            IF(J.GT.L)THEN
               WRITE(TMMLIB_MESSAGE_UNIT,*)'SMART_REAL: STRING_OUT IS ',
     $         'TOO SHORT.'
               STOP
            ENDIF
            STRING_OUT(J:J)='E'
            J=J+1
            IF(J.GT.L)THEN
               WRITE(TMMLIB_MESSAGE_UNIT,*)'SMART_REAL: STRING_OUT IS ',
     $         'TOO SHORT.'
               STOP
            ENDIF
            IF(Y.GE.0)THEN
               STRING_OUT(J:J)='+'
            ELSE
               STRING_OUT(J:J)='-'
            ENDIF
            J=J+1
            IF(J+1.GT.L)THEN
               WRITE(TMMLIB_MESSAGE_UNIT,*)'SMART_REAL: STRING_OUT IS ',
     $         'TOO SHORT.'
               STOP
            ENDIF
            IF(IABS(Y).LT.10)THEN
               STRING_OUT(J:J)='0'
               J=J+1
               WRITE(STRING_OUT(J:J),"(I1)")IABS(Y)
            ELSE
               WRITE(STRING_OUT(J:J+1),"(I2)")IABS(Y)
            ENDIF

C  Medium number.
         ELSE
            IF(Y.GT.0)THEN
               DO I=1,Y
                  J=J+1
                  IF(J.GT.L)THEN
                     WRITE(TMMLIB_MESSAGE_UNIT,*)'SMART_REAL: ',
     $               'STRING_OUT IS TOO SHORT.'
                     STOP
                  ENDIF
                  IF(I.LE.IEND)THEN
                     STRING_OUT(J:J)=STRING(I:I)
                  ELSE
                     STRING_OUT(J:J)='0'
                  ENDIF
                  K=I
               ENDDO
               IF(K.LT.IEND)THEN
                  J=J+1
                  IF(J.GT.L)THEN
                     WRITE(TMMLIB_MESSAGE_UNIT,*)'SMART_REAL: ',
     $               'STRING_OUT IS TOO SHORT.'
                     STOP
                  ENDIF
                  STRING_OUT(J:J)='.'
                  DO I=K+1,IEND
                     J=J+1
                     IF(J.GT.L)THEN
                        WRITE(TMMLIB_MESSAGE_UNIT,*)'SMART_REAL: ',
     $                  'STRING_OUT IS TOO SHORT.'
                        STOP
                     ENDIF
                     STRING_OUT(J:J)=STRING(I:I)
                  ENDDO
               ENDIF
            ELSEIF(Y.LE.0)THEN
               J=J+1
               IF(J+1.GT.L)THEN
                  WRITE(TMMLIB_MESSAGE_UNIT,*)'SMART_REAL: STRING_OUT ',
     $            'IS TOO SHORT.'
                  STOP
               ENDIF
               STRING_OUT(J:J+1)='0.'
               J=J+1
               IF(Y.LT.0)THEN
                  DO I=1,-Y
                     J=J+1
                     IF(J.GT.L)THEN
                        WRITE(TMMLIB_MESSAGE_UNIT,*)'SMART_REAL: ',
     $                  'STRING_OUT IS TOO SHORT.'
                        STOP
                     ENDIF
                     STRING_OUT(J:J)='0'
                  ENDDO
               ENDIF
               DO I=1,IEND
                  J=J+1
                  IF(J.GT.L)THEN
                     WRITE(TMMLIB_MESSAGE_UNIT,*)'SMART_REAL: ',
     $               'STRING_OUT IS TOO SHORT.'
                     STOP
                  ENDIF
                  STRING_OUT(J:J)=STRING(I:I)
               ENDDO
            ENDIF            
         ENDIF
      ENDIF

      END SUBROUTINE SMART_REAL
