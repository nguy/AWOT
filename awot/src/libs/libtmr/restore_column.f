      SUBROUTINE RESTORE_COLUMN(N_EXPANDED,Z_EXPANDED,A_EXPANDED,
     $N_ORIGINAL,Z_ORIGINAL,A_ORIGINAL)

C  Thomas Matejka NOAA/NSSL 17 April 1997

C  This subroutine removes inserted levels in a column.

C  Input:

C  N_EXPANDED is an integer variable that specifies the number of levels
C  in the expanded column.

C  Z_EXPANDED is a one-dimensional real array with elements indexed from
C  1 to N_EXPANDED.  Z_EXPANDED(I) specifies the altitude of the Ith
C  level of the expanded column, counting from the bottom of the column.

C  A_EXPANDED is a one-dimensional real array with elements indexed from
C  1 to N_EXPANDED.  A_EXPANDED(I) specifies the datum at Z_EXPANDED(I)
C  in the expanded column..

C  N_ORIGINAL is an integer variable that specifies the number of levels
C  in the original column.

C  Z_ORIGINAL is a one-dimensional real array with elements indexed from
C  1 to N_ORIGINAL.  Z_ORIGINAL(I) specifies the altitude of the Ith
C  level of the original column, counting from the bottom of the column.

C  Output:

C  A_ORIGINAL is a one-dimensional real array with elements indexed from
C  1 to N_ORIGINAL.  A_ORIGINAL(I) returns the datum at Z_ORIGINAL(I) in
C  the original column.  A_ORIGINAL may be the same as A_EXPANDED, in
C  which case A_EXPANDED is overwritten.

      IMPLICIT NONE
      INCLUDE 'tmrlib.inc'
      INTEGER::N_EXPANDED,N_ORIGINAL,IZ,JZ
      REAL,DIMENSION(N_EXPANDED)::Z_EXPANDED,A_EXPANDED,A_EXPANDED_TEMP
      REAL,DIMENSION(N_ORIGINAL)::Z_ORIGINAL,A_ORIGINAL

C  Copy the expanded array in case it is overwritten.
      DO IZ=1,N_EXPANDED
         A_EXPANDED_TEMP(IZ)=A_EXPANDED(IZ)
      ENDDO

C  Compare the altitudes of the original and expanded data arrays to
C  determine which levels of the expanded data array should be retained.
      DO IZ=1,N_ORIGINAL
         DO JZ=1,N_EXPANDED
            IF(Z_ORIGINAL(IZ).EQ.Z_EXPANDED(JZ))THEN
               A_ORIGINAL(IZ)=A_EXPANDED_TEMP(JZ)
            ENDIF
         ENDDO
      ENDDO

      END SUBROUTINE RESTORE_COLUMN
