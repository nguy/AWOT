C  "<25-Jul-1995 20:33:06><UTC>"
*************************************************************************       
      Integer*4 Function FirstAnyChar(String,AnyChars)
                                                                                
c  Find first occurence of any character from AnyChars in character string.
c   Ex.     FirstAnyChar('001234','32') will equal 4
c   Ex.     FirstAnyChar('0012','32') will equal 4   
c   Ex.     FirstAnyChar('223456','32') will equal 1  
c   Ex.     FirstAnyChar('001234','2') will equal 4
c   Ex.     FirstAnyChar('0012','2') will equal 4   
c   Ex.     FirstAnyChar('223456','2') will equal 1  
c   Ex.     FirstAnyChar('123','4') will equal 0,
c      i.e., zero will be returned if no character is present in string
c      (although one might have logically guessed the string_length + 1
c      would be returned).
c   Ex.     FirstAnyChar('123','48') will equal 0.
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *         
      Implicit none                                                             
      Character String*(*),AnyChars*(*)
      Integer*4 LenString, LenAnyChars,II
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *         
      LenString = Len(String)   
      LenAnyChars = Len(AnyChars)   
      FirstAnyChar = 1
      Do While (FirstAnyChar .le. LenString)
	  Do II = 1, LenAnyChars
              If (String(FirstAnyChar:FirstAnyChar) .eq. 
     >           AnyChars(II:II)) Return
	  End Do
          FirstAnyChar = FirstAnyChar +1
      End Do                                                                    
      FirstAnyChar = 0    ! character not found!!
      Return                                                                    
      End ! Integer*4 Function FirstAnyChar ends
