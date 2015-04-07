C  "<08-Oct-1993 21:54:32><UTC>"
************************************************************************
      Integer*4 Function LastNonChar(String,Char)
                                                                                
c  Find last non-occurence of character in character string.               
c   Ex.     LastNonChar('001234','2') will equal 6
c   Ex.     LastNonChar('0012','2') will equal 3   
c   Ex.     LastNonChar('222222','2') will equal 0,
c      i.e, zero will be returned if character matches all string characters.
c   Ex.     LastNonChar('135','2') will equal 3   
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *         
      Implicit none                                                             
      Character String*(*),Char*1
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *         
      LastNonChar = Len(String)   
      Do While (LastNonChar .gt. 0)
          If (String(LastNonChar:LastNonChar) .ne. Char) Return
          LastNonChar = LastNonChar -1
      End Do                                                                    
      Return                                                                    
      End ! Integer*4 Function LastNonChar ends
