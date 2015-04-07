C  "<08-Oct-1993 21:54:26><UTC>"
****************************************************************************
      Real*4 Function BinAng_Deg(I_Ang)

c This function takes a 16-bit integer and converts to a real number.
c  Please note: "P3 Prog. Des. Manual 5 July 1988" describes seveal items
c     in the ray header as binary angles (16-bit).  If a binary format had
c     been used, then the integer would have been interpreted as:
c         MSB: sign bit
c         2ns MSB: 180 degrees
c         3rd MSB: 90 degrees
c         and so on.  The real number would then be -359.989 to 359.989.
c     However, the Sigmet routines brought back by Jorgensen  July 13, 1988,
c     use a different scheme for what they label binary angles.  Therefore
c     this routine will provide the inverse for the Sigmet routines, the
c     interpretation of the 16-bit integer being:
c         -32768 to 32767 are mapped to  -180.0000 to 179.9945 degrees.
c
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      Integer*2 I_Ang
      Integer*4 J_Ang
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  
      J_Ang = I_Ang
      BinAng_Deg = 90.0 * (Float (J_ang)/16384.0)

      Return

      End ! Function BinAng_Deg ends
