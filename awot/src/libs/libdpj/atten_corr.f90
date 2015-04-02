Subroutine atten_corr(Dbz, Zdr, PhiDP, Cdbz, Czdr, dbz_corr, zdr_corr, Nsteps, Nmax, Bad)

  ! Subroutine to calculate the attenuation correction for dbz and zdr based on PhiDP

  Real, Dimension(Nsteps) :: Dbz, Zdr, PhiDP, Cdbz, Czdr, dbz_corr, zdr_corr, Phi_Smooth
  
  ! Find the first good gate of PhiDP as the minimum of the first 20 gates
  
  Phi_Start = 999.0

  Do i = 1, 20
     If (PhiDP(i) .ne. Bad) Then
        Phi_Start = Amin1(Phi_Start,PhiDP(i))
     End If
  End Do

  ! no good first start gate - exit without atten correction

  If (Phi_Start .eq. 999.0) Then
     Do i = 1, Nmax
        Cdbz(i) = Dbz(i)
        Czdr(i) = Zdr(i)
        dbz_corr(i) = 0.0
        zdr_corr(i) = 0.0
     End Do

     Return
  End If

  ! Replace bad points

  Phi_Good = Phi_Start

  Do i = 1, Nmax
     If (PhiDP(i) .eq. Bad) Cycle
     
     If (PhiDP(i) .lt. 45.0) Then
        PhiDP(i) = Phi_Good
        Cycle
     End If

     Phi_Good = PhiDP(i)
  End Do

  Prev_PhiDP = Phi_Start

  ! if PhiDP is missing use the previous good gate

  Do i = 1, nmax
     If (PhiDP(i) .eq. Bad) Then
        PhiDP_good = Prev_PhiDP
     Else
        PhiDP_good = PhiDP(i)
     End If

     !  The coeffs 0.08 and 0.015 come from Bingi and Chandra's book

     Prev_PhiDP = PhiDP_good
     dbz_corr(i) = 0.08  * (PhiDP_good - Phi_Start)
     zdr_corr(i) = 0.015 * (PhiDP_good - Phi_Start)

     If (dbz(i) .ne. Bad) Then
        Cdbz(i) = dbz(i) + dbz_corr(i)
     Else
        Cdbz(i) = Bad
     End If

     If (Zdr(i) .ne. Bad) Then
        Czdr(i) = Zdr(i) + zdr_corr(i)
     Else
        Czdr(i) = Bad
     End If
  End Do

  Return
End Subroutine atten_corr
