      SUBROUTINE EXTEND_DOWN(N_IN,Z_IN,A_IN,
     $Z_EXTEND,
     $BADDATA,
     $N_OUT,Z_OUT,A_OUT)

C  Thomas Matejka NOAA/NSSL 29 January 2002

C  This subroutine extends data down to a specified altitude if the data
C  at that altitude are missing or if the level is new.

C  Input:

C  N_IN is an integer variable that specifies the number of levels in
C  the input column.

C  Z_IN is a one-dimensional real array with elements indexed from 1 to
C  N_IN.  Z_IN(I) specifies the altitude of the Ith level of the input
C  column, counting from the bottom of the column.

C  A_IN is a one-dimensional real array with elements indexed from 1 to
C  N_IN.  A_IN(I) specifies the datum at Z_IN(I) in the input column..

C  Z_EXTEND is a real variable that specifies the altitude to which to
C  extend data.

C  BADDATA is a real variable that indicates a missing datum.

C  Output:

C  N_OUT is an integer variable that returns the number of levels in the
C  output column.  N_OUT may be the same as N_IN, in which case N_IN is
C  overwritten.

C  Z_OUT is a one-dimensional real array with elements indexed from 1 to
C  N_OUT.  Z_OUT(I) returns the altitude of the Ith level of the output
C  column, counting from the bottom of the column.  Z_OUT may be the
C  same as Z_IN, in which case Z_IN is overwritten.

C  A_OUT is a one-dimensional real array with elements indexed from 1 to
C  N_OUT.  A_OUT(I) returns the datum at Z_OUT(I) in the output column.
C  A_OUT may be the same as A_IN, in which case A_IN is overwritten.

      IMPLICIT NONE
      INCLUDE 'tmrlib.inc'
      INTEGER::N_IN,N_OUT,IZ,JZ,IZ_EXTEND
      REAL::Z_EXTEND,BADDATA
      REAL,DIMENSION(N_IN)::Z_IN,A_IN
      REAL,DIMENSION(N_IN+1)::Z_OUT,A_OUT

C  Insert a new level if necessary.
      CALL INSERT_LEVEL(N_IN,Z_IN,A_IN,
     $Z_EXTEND,BADDATA,
     $BADDATA,
     $N_OUT,Z_OUT,A_OUT)

C  Find the index of the level to extend to.
      DO IZ=1,N_OUT
         IF(Z_OUT(IZ).EQ.Z_EXTEND)THEN
            IZ_EXTEND=IZ
            EXIT
         ENDIF
      ENDDO

C  If the datum at the level to extend to is missing, search upward for
C  the next datum.
      IF(A_OUT(IZ_EXTEND).EQ.BADDATA)THEN
         IF(IZ_EXTEND.LT.N_OUT)THEN
            DO IZ=IZ_EXTEND+1,N_OUT
               IF(A_OUT(IZ).NE.BADDATA)THEN
                  DO JZ=IZ_EXTEND,IZ-1
                     A_OUT(JZ)=A_OUT(IZ)
                  ENDDO
                  RETURN
               ENDIF
            ENDDO
         ENDIF
      ENDIF

      END SUBROUTINE EXTEND_DOWN
