C  "<08-Oct-1993 21:54:29><UTC>"
*************************************************************************       
      Integer*4 Function FirstNonChar(String,Char)
                                                                                
c  Find first non-occurence of character in character string.               
c   Ex.     FirstNonChar('001234','2') will equal 1
c   Ex.     FirstNonChar('2212','2') will equal 3   
c   Ex.     FirstNonChar('333456','3') will equal 4  
c   Ex.     FirstNonChar('333','3') will equal 0,
c      i.e., zero will be returned if character is the only 
c      character present in string (although one might have logically
c      guessed the string_length + 1 would be returned).
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *         
      Implicit none                                                             
      Character String*(*),Char*1
      Integer*4 LenString
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *         
      LenString = Len(String)   
      FirstNonChar = 1
      Do While (FirstNonChar .le. LenString)
          If (String(FirstNonChar:FirstNonChar) .ne. Char) Return
          FirstNonChar = FirstNonChar +1
      End Do                                                                    
      FirstNonChar = 0    ! non-character not found!!
      Return                                                                    
      End ! Integer*4 Function FirstNonChar ends
