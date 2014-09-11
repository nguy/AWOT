C  "<08-Oct-1993 21:54:27><UTC>"
************************************************************************
      Subroutine Casefold(String)

c  This subroutine will convert all lower case letters (ASCII 97-122)
c    in string to upper case (ASCII 65 - 90)

      Implicit none
      Character String*(*)
      Integer*4 Length, II, IASC

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      Length = Len(String)
      Do II = 1,Length
          Iasc = Ichar(String(II:II))
          If ((Iasc .ge. 97) .and. (Iasc .le. 122))
     >        String(II:II) = Char(Iasc-32) ! convert to capital letter
      End Do
      Return
      End !     subroutine Casefold ends
