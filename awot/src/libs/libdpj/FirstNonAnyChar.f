C  "<25-Jul-1995 20:33:12><UTC>"
*************************************************************************       
      Integer*4 Function FirstNonAnyChar(String,AnyChars)
                                                                                
c  Find first non-occurence of any character from AnyChars in character string.
c   Ex.     FirstNonAnyChar('001234','29') will equal 1
c   Ex.     FirstNonAnyChar('2212','24') will equal 3   
c   Ex.     FirstNonAnyChar('333456','34') will equal 5  
c   Ex.     FirstNonAnyChar('001234','2') will equal 1
c   Ex.     FirstNonAnyChar('2212','2') will equal 3   
c   Ex.     FirstNonAnyChar('333456','3') will equal 4  
c   Ex.     FirstNonAnyChar('333','3') will equal 0,
c      i.e., zero will be returned if character is the only 
c      character present in string (although one might have logically
c      guessed the string_length + 1 would be returned).
c   Ex.     FirstNonAnyChar('30303','30') will equal 0
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *         
      Implicit none                                                             
      Character String*(*),AnyChars*(*)
      Integer*4 LenString, LenAnyChars,II
      Logical  CharMatch
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *         
      LenString = Len(String)   
      LenAnyChars = Len(AnyChars)   
      FirstNonAnyChar = 1
      Do While (FirstNonAnyChar .le. LenString)
	  CharMatch = .False.
	  Do II = 1,LenAnyChars
              If (String(FirstNonAnyChar:FirstNonAnyChar)
     >           .eq. AnyChars(II:II)) CharMatch = .True.
	  End Do
          If (.Not. CharMatch) Return
          FirstNonAnyChar = FirstNonAnyChar +1
      End Do                                                                    
      FirstNonAnyChar = 0    ! non-character not found!!
      Return                                                                    
      End ! Integer*4 Function FirstNonAnyChar ends
