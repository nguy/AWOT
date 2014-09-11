C  "<08-Oct-1993 21:54:31><UTC>"
*************************************************************************       
      Integer*4 Function LastChar(String,Char)
                                                                                
c  Find last occurence of character in character string.               
c   Ex.     LastChar('001234','2') will equal 4
c   Ex.     LastChar('0012','2') will equal 4   
c   Ex.     LastChar('223456','2') will equal 2  
c   Ex.     LastChar('34567','2') will equal 0,
c      i.e., zero will be returned if character is not in string.
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *         
      Implicit none                                                             
      Character String*(*),Char*1
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *         
      LastChar = Len(String)   
      Do While (LastChar .gt. 0)
          If (String(LastChar:LastChar) .eq. Char) Return
          LastChar = LastChar -1
      End Do                                                                    
      Return                                                                    
      End ! Integer*4 Function LastChar ends
