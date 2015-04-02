      Subroutine Qsort (A, B, N)

C SORTS N A-B     RECORDS ASCENDING USING THE QUICKEST SORT ALGORITHM
C KNOWN TO MANKIND.
C A - ARRAY OF LENGTH N WHICH IS THE SORT KEY
C B - ARRAY OF LENGTH N WHICH IS PART OF THE SORT RECORD
C N - NUMBER OF RECORDS TO BE SORTED
C      MAXIMUM N IS 2**(DIMENSIONS OF LV AND IUV). IF THOSE ARRAYS ARE
C DIMENSIONED 18, THE MAXIMUM NUMBER OF SORT RECORDS IS 2**18, OR 262144
C VERSION 10 NOV 1976
C
      DIMENSION A(N), B(N)
      DIMENSION LV(18), IUV(18)
C                                      INITIALIZE PUSHDOWN LIST
      LV(1) = 1
      IUV(1) = N
      IP = 1
C                                      PARTITION NEXT SEGMENT
 10   IF (IP .LT. 1) RETURN
 20   IF (IUV(IP) - LV(IP) .GE. 1) GO TO 30
      IP = IP - 1
      GO TO 10
 30   LP = LV(IP) - 1
      IUP = IUV(IP)
C                                      CHOOSE BOUND
      A1 = A(IUP)
      B1 = B(IUP)
C                                      MOVE LOWER POINTER
 40   IF (IUP - LP .LT. 2) GO TO 70
      LP = LP + 1
      IF (A(LP) .LE. A1) GO TO 40
      A(IUP) = A(LP)
      B(IUP) = B(LP)
C                                      MOVE UPPER POINTER
 50   IF (IUP - LP .LT. 2) GO TO 60
      IUP = IUP - 1
      IF (A(IUP) .GE. A1) GO TO 50
      A(LP) = A(IUP)
      B(LP) = B(IUP)
      GO TO 40
C                                      FINISH UP
 60   IUP = IUP - 1
 70   A(IUP) = A1
      B(IUP) = B1
      IF (IUP - LV(IP) .LT. IUV(IP) - IUP) GO TO 80
      LV(IP + 1) = IUP + 1
      IUV(IP + 1) = IUV(IP)
      IUV(IP) = IUP - 1
      GO TO 90
 80   LV(IP + 1) = LV(IP)
      IUV(IP + 1) = IUP - 1
      LV(IP) = IUP + 1
 90   IP = IP + 1
      GO TO 20
      END
