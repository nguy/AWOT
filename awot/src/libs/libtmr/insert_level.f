      SUBROUTINE INSERT_LEVEL(N_IN,Z_IN,A_IN,
     $Z_INSERT,A_INSERT,
     $BADDATA,
     $N_OUT,Z_OUT,A_OUT)

C  Thomas Matejka NOAA/NSSL 17 April 1997

C  This subroutine inserts a datum at a specified altitude into a
C  column.  If the altitude at which the datum is to be inserted already
C  exists and if the datum in that level is not missing, the column is
C  unchanged.  If the altitude at which the datum is to be inserted
C  already exists and if the datum in that level is missing, the
C  specified value is inserted in the existing level.  If the altitude
C  at which the datum is to be inserted does not exist, the specified
C  datum is inserted into a new level and the data in the column are
C  shifted.

C  Input:

C  N_IN is an integer variable that specifies the number of levels in
C  the input column.

C  Z_IN is a one-dimensional real array with elements indexed from 1 to
C  N_IN.  Z_IN(I) specifies the altitude of the Ith level of the input
C  column, counting from the bottom of the column.

C  A_IN is a one-dimensional real array with elements indexed from 1 to
C  N_IN.  A_IN(I) specifies the datum at Z_IN(I) in the input column..

C  Z_INSERT is a real variable that specifies the altitude at which a
C  data value is to be inserted.

C  A_INSERT is a real variable that specifies the datum to insert.

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
      INTEGER::N_IN,N_OUT,IZ,JZ,N_IN_TEMP
      REAL::BADDATA,Z_INSERT,A_INSERT
      REAL,DIMENSION(N_IN)::Z_IN,A_IN,Z_IN_TEMP,A_IN_TEMP
      REAL,DIMENSION(N_IN+1)::Z_OUT,A_OUT

C  Copy the input variables in case they are overwritten.
      N_IN_TEMP=N_IN
      DO IZ=1,N_IN
         Z_IN_TEMP(IZ)=Z_IN(IZ)
         A_IN_TEMP(IZ)=A_IN(IZ)
      ENDDO

C  Check whether the level of insertion already exists.
      DO IZ=1,N_IN_TEMP
         IF(Z_IN_TEMP(IZ).EQ.Z_INSERT)THEN
            N_OUT=N_IN_TEMP
            DO JZ=1,N_IN_TEMP
               Z_OUT(JZ)=Z_IN_TEMP(JZ)
               A_OUT(JZ)=A_IN_TEMP(JZ)
            ENDDO
            IF(A_OUT(IZ).EQ.BADDATA)THEN
               A_OUT(IZ)=A_INSERT
            ENDIF
            RETURN
         ENDIF
      ENDDO

C  The level of insertion is new.
      N_OUT=N_IN_TEMP+1

C  Check whether the level of insertion is below the top of the column.
      DO IZ=1,N_IN_TEMP
         IF(Z_IN_TEMP(IZ).LT.Z_INSERT)THEN
            Z_OUT(IZ)=Z_IN_TEMP(IZ)
            A_OUT(IZ)=A_IN_TEMP(IZ)
         ELSE
            Z_OUT(IZ)=Z_INSERT
            A_OUT(IZ)=A_INSERT
            DO JZ=IZ+1,N_OUT
               Z_OUT(JZ)=Z_IN_TEMP(JZ-1)
               A_OUT(JZ)=A_IN_TEMP(JZ-1)
            ENDDO
            RETURN
         ENDIF
      ENDDO

C  The level of insertion is above the top of the column.
      Z_OUT(N_OUT)=Z_INSERT
      A_OUT(N_OUT)=A_INSERT

      END SUBROUTINE INSERT_LEVEL
