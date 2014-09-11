C  "<26-Jul-1995 17:36:22><UTC>"
**************************************************************************
      Subroutine GetToken(String,NumToken,DelimChar,IPosition) 
c Get starting position of token in string.
c
c Input conditions
c   String: string to search
c   NumToken: which token to get, >=1.
c   Delimchar: characters which can be delimiters.
c   IPosition: undefined.
c Exit conditions:
c   String, NumToken, and DelimChar are not changed.
c   IPosition:
c       if = 0, then NumToken is <1 or Token could not be found.
c       if >0, then the Token begins in String(Iposition:Iposition).

c  Examples:
c    Call GetToken(',abc d', 2, ', ',Iposition) sets IPosition to 6
c    Call GetToken(',abc d',-1, ', ',Iposition) sets IPosition to 0
c    Call GetToken(',abc d', 2, ' ' ,Iposition) sets IPosition to 6
c    Call GetToken(',abc d', 3, ' ' ,Iposition) sets IPosition to 0
c    Call GetToken('hi 30 ',2,' ',Iposition) sets IPosition to 4
c    Call GetToken(' hi 30 ',2,' ',Iposition) sets IPosition to 5
c    Call GetToken(' hi 30 ',1,' ',Iposition) sets IPosition to 2
c    Call GetToken(',',      1,',',Iposition) sets IPosition to 0
c    Call GetToken('A',      1,' ',Iposition) sets IPosition to 1
c    Call GetToken('avd30def442gh9',3,'0123456789',IPosition) sets IPosition to 12
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      Implicit None
      Character String*(*),DelimChar*(*)
      integer *4 FirstNonAnyChar,FirstAnyChar
      Integer*4   LenString,NumToken,IPosition,II,ia,ib,iaprime,ibprime
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      IPosition = 0  ! in case the desired token not found
      if (NumToken .lt. 1) Return
      LenString = Len(String)

c ia initialized to start of string; each iteration of loop changes it
c    to point just past the current token being counted.
c ib will be set to start of current token being counted.
c iaprime and ibprime are pointers relative to ib and ia, respectively.

      ia = 1
      Do II = 1,NumToken
	 ibprime = FirstNonAnyChar(String(ia:),DelimChar)
	 if (ibprime .le. 0) return  ! nothing left but delim chars
	 ib = ibprime + ia -1
	 iaprime = FirstAnyChar(String(ib:),DelimChar)
	 if (II.ne.NumToken) Then
	     If (iaprime .le. 0) return ! no more delim, hence no more Tokens
	 end if
	 ia = iaprime + ib -1   ! start of delim chars after current token
      end do
      IPosition = ib
      return
      end ! GetToken ends
