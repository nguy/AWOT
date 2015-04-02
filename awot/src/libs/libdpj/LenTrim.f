C  "<08-Oct-1993 21:54:32><UTC>"
*************************************************************************       
      Integer*4 Function Lentrim (String)                                       
                                                                                
c  Find last nonblank character of character string.                            
c   Ex.     LenTrim('  12  ') will equal 4                                      
c   Ex.     LenTrim('  12') will equal 4                                        
c   Ex.     LenTrim('12    ') will equal 2                                      
c   Ex.     LenTrim('      ') will equal 0                                      
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *         
      Implicit none                                                             
      Character String*(*)                                                      
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *         
      LenTrim = Len(String)                                                     
      Do While (LenTrim .gt. 0)                                                 
          If (String(LenTrim:LenTrim) .ne. ' ') Return                          
          LenTrim = LenTrim -1                                                  
      End Do                                                                    
      Return                                                                    
      End ! Integer*4 Function LenTrim ends
