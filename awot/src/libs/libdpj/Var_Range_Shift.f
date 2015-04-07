C  "<08-Oct-1993 21:54:37><UTC>"
**************************************************************************
      Subroutine  Var_range_shift(Iarray,Ind,I1,I2,Incre)
c
c  This subroutine will do shifting and possible filtering for array
c     Iarray (Integer*2).
c  Input
c     Ind: number, after incrementing, for next slot in Iarray to
c         shift to.
c     The slots indexed by I1..I2, with increment Incre, will be shifted.
c     Incre should be >=1.
c     For example
c       On entry:
c         Iarray(1..16) = 2,4,6,8,  10,12,14,0,  -, -, -, -,  -, -, -, -,
c         Ind =8, I1=4, I2=8, Incre =2
c       On exit:
c         Iarray(1..16) = 2,4,6,8,  10,12,14,0,  8,12, 0, -,  -, -, -, -,
c         Ind =11, I1=4, I2=8, Incre =2
c
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  
      Integer*2  Iarray(*)
 
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      Do II = I1,I2,Incre
          Ind =Ind+1
          Iarray(Ind) = Iarray(II)
      End Do
      Return
      End ! subroutine Var_range_Shift ends
