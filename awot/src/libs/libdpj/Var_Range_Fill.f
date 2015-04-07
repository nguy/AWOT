C  "<08-Oct-1993 21:54:36><UTC>"
**************************************************************************
      Subroutine  Var_Range_Fill(Iarray,Ind,I1,I2,Ifill)
c
c  This subroutine will do filling and shifting for integer*2 Iarray.
c  Input
c    I1 should be < I2.
c
c     Ind: number, after incrementing, for next slot in Iarray.
c     The slots indexed by I1..I2 will be shifted, but before each shift
c     there will be  Ifill points filled between slot and previous.
c     When filling is done between two points, the point with the lower
c     index will be used.
c     If Ifill is less than 0:
c         The process goes in reverse (this would be used when an array is
c         being expanded, say Iarray(5..8) to Iarray(5..12), and care has
c         be taken so that expansion does not write over good data).  In this
c         case, Ind, after decrementing, will be next slot in Iarray.
c         The slots are indexed by I1..I2 will be shifted, starting with
c         I2 and going down.  AFTER each shift there will be Abs(Ifill)
c         points between point and next lower Point.
c           For example
c             On entry:
c               Iarray(1..16) = 2,4,6,8,  10,12,14,0,  -, -, -, -,  -, -, -, -
c               Ind =13, I1=5, I2=8, Ifill =-1
c             On exit:
c               Iarray(1..16) = 2,4,6,8,  8,10,10,12, 12,14, 0, 0,  -, -, -, -
c               Ind =5, I1=5, I2=8, Ifill = -1
c
c     If there is filling between two points and one or both are
c     zero (missing data), then all filled points will be zero.
c     For example
c       On entry:
c         Iarray(1..16) = 2,4,6,8,  10,12,14,0,  -, -, -, -,  -, -, -, -,
c         Ind =8, I1=5, I2=8, Ifill =1
c       On exit:
c         Iarray(1..16) = 2,4,6,8,  10,12,14,0,  8,10,10,12, 12,14, 0, 0,
c         Ind =16, I1=5, I2=8, Ifill =1
c
c     Note:  Iarray (I1-1,..,I2) needs to be defined on entry.
c         If Ifill >=0, then the targets will be:
c             Iarray (Ind+1,..,Ind+(I2-I1+1)*(1+Ifill) )
c             and on exit Ind = Ind+(I2-I1+1)*(1+Ifill).
c             The target slots will be filled in increasing order.
c         If Ifill <0, then the targets will be:
c             Iarray (Ind+(I2-I1+1)*(-1 +Ifill),..,Ind-1) and
c             on exit Ind = Ind+(I2-I1+1)*(-1+Ifill).
c             The target slots will be filled in decreasing order.
c
* *  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  
      Integer*2  Iarray(*)
  
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      If (Ifill.ge.0) Then
  
c     Fill up slots with increasing index
 
          Ichange = 1
          Istart = I1
          Iend = I2
          Ifillp = Ifill +1
  
      Else
  
c     Fill up slots with decreasing index
  
          Ichange = -1   ! This will change direction of filling
          Istart = I2
          Iend = I1
          Ifillp = Ifill -1
      End If
      Do II = Istart,Iend,Ichange
          If (Ichange .eq. -1) Then
  
c  Switch point before the filled points have been put in.
  
              Ind =Ind+Ichange
              Iarray(Ind) = Iarray(II)
          End If
          JJ= 0
          Do While (Ichange*JJ.lt.Ichange*Ifill)
  
c  Get filled points
  
              JJ=JJ+IChange
              Ind=Ind+Ichange
              If ((Iarray(II-1).eq.0).or.(Iarray(II).eq.0)) Then
                  Iarray(Ind) =0
              Else
                  Iarray(Ind)= Iarray(II-1)  ! fill with lower-indexed element
              End If
  
          End Do
          If (Ichange .eq. 1) Then
  
c  Switch point after the filled points have been put in.
  
              Ind =Ind+Ichange
              Iarray(Ind) = Iarray(II)
          End If
      End Do
      Return
      End !  subroutine Var_Range_Fill ends
