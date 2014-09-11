C  "<12-May-1999 14:30:21><MDT>"
*************************************************************************       
      Integer*4 Function FirstNonNum(String)
                                                                                
c  Find first non-number in character string.               
c   Ex.     FirstNonNum('B012d4') will equal 1
c   Ex.     FirstNonNum(' 912d4') will equal 1
c   Ex.     FirstNonNum('22p2') will equal 3   
c   Ex.     FirstNonNum('12345F12') will equal 6  
c   Ex.     FirstNonNum('01234567890') will equal 0,
c      i.e., zero will be returned if only numeric characters are
c      present in string (although one might have logically
c      guessed the string_length + 1 would be returned).
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *         
      Implicit none                                                             
      Character String*(*),Char*1
      Integer*4 LenString
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *         
      LenString = Len(String)   
      FirstNonNum = 1
      Do While (FirstNonNum .le. LenString)
	  Char = String(FirstNonNum:FirstNonNum)
          If ( ('0' .eq. Char) .or.
     >         ('1' .eq. Char) .or.
     >         ('2' .eq. Char) .or.
     >         ('3' .eq. Char) .or.
     >         ('4' .eq. Char) .or.
     >         ('5' .eq. Char) .or.
     >         ('6' .eq. Char) .or.
     >         ('7' .eq. Char) .or.
     >         ('8' .eq. Char) .or.
     >         ('9' .eq. Char)) then
              FirstNonNum = FirstNonNum +1
	  else 
	      Return
	  end if
      End Do                                                                    
      FirstNonNum = 0    ! non-numeric not found!!
      Return                                                                    
      End ! Integer*4 Function FirstNonNum ends
