      SUBROUTINE DOP_SOLN_DRIVER_3_SIMPLE(N_OBS,USE,DIRCOS_X1,DIRCOS_X2,
     $DIRCOS_X3,VR,SD_VR,BADDATA,VP1,SD_VP1,VP2,SD_VP2,VP3,SD_VP3)

C  Thomas Matejka NOAA/NSSL 24 August 2001

C  This subroutine calculates Cartesian components of target velocity
C  and their standard deviations from colocated Doppler radar data at
C  one point using a three-component solution.

C  Input:

C  N_OBS is an integer variable that specifies the number of colocated
C  velocity observations.

C  USE is a one-dimensional logical array.  USE(I) should be .TRUE. if
C  and only if the Ith velocity observation is to be used in the
C  calculation.

C  DIRCOS_X1 is a one-dimensional real array.  DIRCOS_X1(I) specifies
C  the direction cosine between the radar beam and the first Cartesian
C  coordinate axis for the Ith velocity observation.

C  DIRCOS_X2 is a one-dimensional real array.  DIRCOS_X2(I) specifies
C  the direction cosine between the radar beam and the second Cartesian
C  coordinate axis for the Ith velocity observation.

C  DIRCOS_X3 is a one-dimensional real array.  DIRCOS_X3(I) specifies
C  the direction cosine between the radar beam and the third Cartesian
C  coordinate axis for the Ith velocity observation.

C  VR is a one-dimensional real array.  VR(I) specifies the Ith velocity
C  observation.

C  SD_VR is a one-dimensional real array.  SD_VR(I) specifies the
C  standard deviation of VR(I).

C  BADDATA is a real variable that indicates a missing value as
C  described.

C  Output:

C  VP1 is a real variable that returns the first Cartesian component of
C  target velocity.  If it is not calculated, it is returned as BADDATA.

C  SD_VP1 is a real variable that returns the standard deviation of VP1.
C  If it is not calculated, it is returned as BADDATA.

C  VP2 is a real variable that returns the second Cartesian component of
C  target velocity.  If it is not calculated, it is returned as BADDATA.

C  SD_VP2 is a real variable that returns the standard deviation of VP2.
C  If it is not calculated, it is returned as BADDATA.

C  VP3 is a real variable that returns the third Cartesian component of
C  target velocity.  If it is not calculated, it is returned as BADDATA.

C  SD_VP3 is a real variable that returns the standard deviation of VP3.
C  If it is not calculated, it is returned as BADDATA.

      IMPLICIT NONE
      INTEGER::N_OBS
      LOGICAL,DIMENSION(N_OBS)::USE
      REAL::VP1_3,VP2_3,VP3_3,SD_VP1_3,SD_VP2_3,SD_VP3_3,BADDATA,VP1,
     $VP2,VP3,SD_VP1,SD_VP2,SD_VP3
      REAL,DIMENSION(N_OBS)::DIRCOS_X1,DIRCOS_X2,DIRCOS_X3,VR,SD_VR

C  Perform a three-component solution.
      CALL DOP_3_COMP_SOLN_SIMPLE(N_OBS,USE,DIRCOS_X1,DIRCOS_X2,
     $DIRCOS_X3,VR,SD_VR,BADDATA,VP1_3,SD_VP1_3,VP2_3,SD_VP2_3,VP3_3,
     $SD_VP3_3)
      VP1=VP1_3
      SD_VP1=SD_VP1_3
      VP2=VP2_3
      SD_VP2=SD_VP2_3
      VP3=VP3_3
      SD_VP3=SD_VP3_3

      END SUBROUTINE DOP_SOLN_DRIVER_3_SIMPLE
