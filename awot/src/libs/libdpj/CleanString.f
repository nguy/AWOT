C  "<08-Oct-1993 21:54:29><UTC>"
*************************************************************************       
      Subroutine CleanString(String,Marker,Clean)
c  Clean up a string (usually blanking out) all characters from the
c     "Marker" character onward.
c  Input conditions:
c     String: character string to clean
c     Marker: character string of length one; first occurence in the
c         string will cause from that position to end of string to
c         be replace by the Clean character.
c     Clean: character string of length one to use as a replacement.
c  Example of call
c     Call CleanString(LineRead,'!',' ') ! ignore after '!' in line
c  Exit Conditions:
c     Marker and Clean will be unchanged.
c     String will be changed if-and-only-if the Marker character was
c        found in the string.
      Implicit None
      Character String*(*), Marker*1, Clean*1
      Integer*4 Ilength, ii
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      Ilength = Len(String)
      ii = 1
      Do While (ii .le. Ilength)
          if (String(ii:ii) .eq. Marker)then
              Do While (ii .le. Ilength)
		  String(ii:ii) = Clean
	          ii = ii + 1
	      end do
	  end if
	  ii = ii + 1
      end do
      return
      end !  CleanString ends
