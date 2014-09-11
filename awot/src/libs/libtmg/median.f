      REAL FUNCTION MEDIAN(A,N)

C  Thomas Matejka NOAA/NSSL 13 January 1995

C  This function returns the median of the N elements of array A.

      IMPLICIT NONE
      INTEGER I,N,NDUM
      REAL A(N)
      REAL B(N)

      DO I=1,N
         B(I)=A(I)
      ENDDO
      CALL SORT_F(B,N,.FALSE.,NDUM)
      IF(MOD(N,2).EQ.0)THEN
         MEDIAN=(B(N/2)+B(N/2+1))/2.
      ELSE
         MEDIAN=B(N/2+1)
      ENDIF
      RETURN
      END
