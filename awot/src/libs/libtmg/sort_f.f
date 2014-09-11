      SUBROUTINE SORT_F(A,N,ELIM_DUPLICATES,M)

C  Thomas Matejka NOAA/NSSL 3 July 1995

C  This subroutine sorts the N elements of array A in ascending order.

C  If ELIM_DUPLICATES is .TRUE., then duplicate values of A are
C  eliminated.

C  The number of elements in the sorted list is returned as M.

      IMPLICIT NONE
      LOGICAL ELIM_DUPLICATES
      INTEGER J,K,K_MIN,N,M
      REAL STORE,A_MIN
      REAL A(N)

C  Sort the list.
      IF(N.GT.1)THEN
         DO J=1,N-1
            A_MIN=A(J)
            K_MIN=J
            DO K=J+1,N
               IF(A(K).LT.A_MIN)THEN
                  A_MIN=A(K)
                  K_MIN=K
               ENDIF
            ENDDO
            IF(K_MIN.NE.J)THEN
               STORE=A(J)
               A(J)=A(K_MIN)
               A(K_MIN)=STORE
            ENDIF
         ENDDO

C  Eliminate duplicates.
         IF(ELIM_DUPLICATES)THEN
            M=N
            J=1
            DOWHILE(J.LT.M)
               IF(A(J).EQ.A(J+1))THEN
                  IF(J+1.LT.M)THEN
                     DO K=J+1,M-1
                        A(K)=A(K+1)
                     ENDDO
                  ENDIF
                  M=M-1
               ELSE
                  J=J+1
               ENDIF
            ENDDO
         ELSE
            M=N
         ENDIF
      ELSE
         M=N
      ENDIF

C  Done.
      RETURN
      END

