C  "<08-Oct-1993 21:54:24><UTC>"
*************************************************************************       
      Subroutine ASCII32_126(String,SubChar,Icount)                             
                                                                                
c  Search through string, and all characters with ASCII <32 or                  
c	ASCII > 126 are replaced by the character SubChar.                            
c  Icount is Integer*4. On exit it will equal number of substitutions made.     
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *         
      Implicit none                                                             
      Character String*(*),SubChar*1                                            
      Integer*4 Ipoint,Iasc, Icount                                             
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *         
      Icount = 0                                                                
      Ipoint  = Len(String)                                                     
      Do While (Ipoint .gt. 0)                                                  
          Iasc = Ichar(String(Ipoint:Ipoint))                                   
          If((Iasc.lt.32).or.(Iasc.gt.126))Then                                 
              Icount = Icount + 1                                               
              String(Ipoint:Ipoint)=SubChar                                     
          End If                                                                
          Ipoint = Ipoint -1                                                    
      End Do                                                                    
      Return                                                                    
      End !  Subroutine ASCII32_126  ends                                     
