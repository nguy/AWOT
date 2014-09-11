C  "<08-Oct-1993 21:54:29><UTC>"
*************************************************************************       
      Integer*4 Function FirstChar(String,Char)
                                                                                
c  Find first occurence of character in character string.               
c   Ex.     FirstChar('001234','2') will equal 4
c   Ex.     FirstChar('0012','2') will equal 4   
c   Ex.     FirstChar('223456','2') will equal 1  
c   Ex.     FirstChar('123','4') will equal 0,
c      i.e., zero will be returned if character is not present in string
c      (although one might have logically guessed the string_length + 1
c      would be returned).
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *         
      Implicit none                                                             
      Character String*(*),Char*1
      Integer*4 LenString
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *         
      LenString = Len(String)   
      FirstChar = 1
      Do While (FirstChar .le. LenString)
          If (String(FirstChar:FirstChar) .eq. Char) Return
          FirstChar = FirstChar +1
      End Do                                                                    
      FirstChar = 0    ! character not found!!
      Return                                                                    
      End ! Integer*4 Function FirstChar ends
