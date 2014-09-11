*********************************************************************
      Integer*4 Function FindPattern(pattern,string)
!    "<20-May-1994 19:06:20><UTC>"
c
c This is function to check if a pattern can be found in the string 
c (NO case folding).
c Input conditions:
c    pattern: character string pattern to search for in string.
c    string: character string
c Output of function:
c    0 is returned if pattern was not found.
c    >0 = n will be returned to mean the first occurrence of the
c       pattern begins in the n-th character position in the string.
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
      Implicit none
      Character pattern*(*),String*(*)
      Integer*4 ipatlen,istringlen,ii
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
      ipatlen = len(pattern)
      istringlen = len (string)
      FindPattern = 0
      ii = 0
      Do While (ii + ipatlen .le. istringlen)
	  ii = ii + 1
	  if (pattern .eq. string(ii:ii+ipatlen-1))then
	      FindPattern = ii
	      return
	  end if
      end do
      return
      end ! FindPattern ends
*********************************************************************
