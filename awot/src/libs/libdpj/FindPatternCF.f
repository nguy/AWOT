*********************************************************************
      Integer*4 Function FindPatternCF(pattern,string)
!    "<20-May-1994 19:06:26><UTC>"
c
c This is function to check if a pattern can be found in the string 
c (Case Folding for a-z to A-Z, ASCII).
c Input conditions:
c    pattern: character string pattern to search for in string.
c    string: character string
c Output of function:
c    0 is returned if pattern was not found.
c    >0 = n will be returned to mean the first occurrence of the
c       pattern begins in the n-th character position in the string.
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
      Implicit none
      Character pattern*(*),String*(*), charp*1,chars*1
      Integer*4 ipatlen,istringlen,ii,jj,kk
      logical Match
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
      ipatlen = len(pattern)
      istringlen = len (string)
      FindPatternCF = 0
      ii = 0
      Do While (ii + ipatlen .le. istringlen)
	  ii = ii + 1
	  match = .true.
	  jj = 0
	  kk = ii - 1
	  Do While (match .and. (jj .lt.ipatlen))
	      jj = jj + 1
	      kk = kk + 1
              charp = pattern(jj:jj)
	      call casefold(charp)
              chars = string(kk:kk)
	      call casefold(chars)
	      match = (charp .eq. chars)
	  end do
	  if (match)Then
	      FindPatternCF = ii
	      return
	  end if
      end do
      return
      end ! FindPatternCF ends
*********************************************************************
