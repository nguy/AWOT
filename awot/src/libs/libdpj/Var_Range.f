C  "<08-Oct-1993 21:54:36><UTC>"
*************************************************************************       
      Subroutine Var_range (Iarray,Nbins)
  
c  This subroutine will take the 512 points in Integer*2 array Iarray and use
c     either filling or discarding to achieve equal spacing.
c         The original data was
c             256 samples at 75 m: 19.2 km
c             128 samples at 150 m: 19.2 km
c             128 samples at 300 m: 38.4 km
c         There are three possible values for Nbins
c             Nbins= 256  (discards 192 samples in first section and
c                 64 samples in middle secton) with 300 m spacing.
c             Nbins= 512  (discards 128 samples in first section and
c                 fills for 128 new bins in last section) at 150 m.
c             Nbins=1024 (fills for 128 new bins in second
c                 section and for 384 new bins in last section.
c         If Nbins is not one of these three values, then the default
c         of 512 will be used.
c     If filling is done between two values, one of which is
c         not defined (=0), then all the filled values will
c         be set =0.
c     Filling is done instead of interpolation because of problems
c         where adjacent bins do not have values close to each other,
c         e.g., with velocity folding.
c
c  On exit the array Iarray will have had the equal spacing
c     put in; instead of 512 points being defined for each array Nbins
c     will be defined (If Nbins is not one of the three values, then 512
c     bins at equal spacing will be defined).  Iarray will be
c     undefined for bins beyond the equal spaced bins.
* *  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  
      Integer*2  Iarray(*)
  
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  
      Imode = 2 ! default mode, output = 512 bins
      If (Nbins .eq. 256) Then
          Imode = 1 ! output = 256 bins
      Else If (Nbins .eq. 1024) Then
          Imode = 3 ! output =1024 bins
      End If
  
      If (Imode .eq.1) Then
  
c  IMode 1: output is 256 points at 300 m. ************************
  
c  Use only 1 of 4 points for the 256 data at 75m.
  
          Ind = 0
          Call Var_Range_Shift(Iarray,Ind,4,256,4)
  
c  Use 1 of 2 points for the 128 data at 150m.
  
          Call Var_range_shift(Iarray,Ind,258,384,2)
  
c  Use all points for the 128 data at 300m, just shift them over.
  
          Call Var_range_shift(Iarray,Ind,385,512,1)
 
      Else If (Imode .eq. 2) then
  
c  IMode 2: output is 512 points at 150 m. ************************
  
c  Use only 1 of 2 points for the 256 data at 75m.
  
          Ind = 0
          Call Var_Range_Shift(Iarray,Ind,2,256,2)
 
c  Use all points for the 128 data at 150m, just shift them over.
  
          Call Var_range_shift(Iarray,Ind,257,384,1)
  
c  Use all points for the 128 data at 300m plus fill for another 128.
  
          Call Var_Range_Fill(Iarray,Ind,385,512,1)
  
      Else
 
c  IMode 3: output is 1024 points at 75 m. ************************
  
c  Use all of the 256 data at 75m and no shifting is needed.
  
  
c  Use all points for the 128 data at 300m plus fill for another
c     3*128 = 384 points.  This must be done before the middle section,
c     because the middle secton will expand over the old location of
c     the third section.
  
          Ind = 512
          Call Var_Range_Fill(Iarray,Ind,385,512,3)
  
c  Use all points for the 128 data at 150m plus fill for another
c     128 points.  This is for the middle section. The order of filling has to
c     be switched so that the data is not erased before it is used.
  
          Ind = 513
          Call Var_Range_Fill(Iarray,Ind,257,384,-1)
  
      End If
      If((Nbins .ne. 512) .and. (Imode .eq.2)) Then
  
c  Since this occurred by forcing the default of 512 points, zero out the
c     remainder of the array, as calling routine may be a little wild.
  
          Do II = 513,1024
              Iarray(II) = 0
          End Do
      End If
      Return
      End !  subroutine Var_range ends
