C  "<08-Oct-1993 21:54:34><UTC>"
*********************************************************************
      Integer*4 Function String_len_c(string)
c Find number of characters in string that are before a null byte which
c  is used in c for string terminator.
c Return value is >=0 and <= length of string in Fortran sense.
      implicit none
      integer*4 length
      character string*(*)
      character NULL_BYTE*(1)
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      NULL_BYTE = char(0)
      length = len (string) ! get Fortran string length 
      string_len_c = 1
      do while (string_len_c .le. length)
	  if (string(string_len_c:string_len_c) .eq. NULL_BYTE)then
	      string_len_c = string_len_c - 1
	      return
	  end if
	  string_len_c = string_len_c + 1
      end do
      string_len_c = string_len_c - 1
      return
      end !  Integer*4 Function String_len_c ends
