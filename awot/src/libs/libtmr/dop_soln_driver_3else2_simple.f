      SUBROUTINE DOP_SOLN_DRIVER_3ELSE2_SIMPLE(SUBSTITUTE_VP3,N_OBS,USE,
     $DIRCOS_X1,DIRCOS_X2,DIRCOS_X3,VR,SD_VR,VP3_ASSUMED,SD_VP3_ASSUMED,
     $BADDATA,VP1,SD_VP1,VP2,SD_VP2,VP3,SD_VP3)

C  Thomas Matejka NOAA/NSSL 24 August 2001

C  This subroutine calculates Cartesian components of target velocity
C  and their standard deviations from colocated Doppler radar data at
C  one point using a three-component solution if possible and otherwise
C  using a two-component solution and an assumed third Cartesian
C  component.

C  Input:

C  SUBSTITUTE_VP3 is a logical variable.  If it is .FALSE., when the
C  third Cartesian component is not solved for, the third Cartesian
C  component and its standard deviation are returned as BADDATA.  If it
C  is .TRUE., when the third Cartesian component is not solved for, and
C  when the first or second Cartesian component has been solved for, the
C  third Cartesian component and its standard deviation are returned as
C  the assumed values.

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

C  VP3_ASSUMED is a real variable that specifies an assumed third
C  Cartesian component of the target velocity.  If it is missing, it
C  should equal BADDATA.

C  SD_VP3_ASSUMED is a real variable that specifies the standard
C  deviation of VP3_ASSUMED.

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
      LOGICAL::SUBSTITUTE_VP3
      LOGICAL,DIMENSION(N_OBS)::USE
      INTEGER::I_OBS
      REAL::VP1_3,VP2_3,VP3_3,VP1_2,VP2_2,SD_VP1_3,SD_VP2_3,SD_VP3_3,
     $SD_VP1_2,SD_VP2_2,BADDATA,VP3_ASSUMED,SD_VP3_ASSUMED,VP1,VP2,VP3,
     $SD_VP1,SD_VP2,SD_VP3
      REAL,DIMENSION(N_OBS)::DIRCOS_X1,DIRCOS_X2,DIRCOS_X3,VR,SD_VR,VR2,
     $SD_VR2

C  Perform a three-component solution.
      CALL DOP_3_COMP_SOLN_SIMPLE(N_OBS,USE,DIRCOS_X1,DIRCOS_X2,
     $DIRCOS_X3,VR,SD_VR,BADDATA,VP1_3,SD_VP1_3,VP2_3,SD_VP2_3,VP3_3,
     $SD_VP3_3)

C  If there is an assumed third Cartesian component, reduce the
C  observations by its contribution and calculate a new standard
C  deviation.
      IF(VP3_ASSUMED.NE.BADDATA)THEN
         DO I_OBS=1,N_OBS
            VR2(I_OBS)=VR(I_OBS)-VP3_ASSUMED*DIRCOS_X3(I_OBS)
            SD_VR2(I_OBS)=SQRT(SD_VR(I_OBS)**2+DIRCOS_X3(I_OBS)**2*
     $      SD_VP3_ASSUMED**2)
         ENDDO

C  Perform a two-component solution.
         CALL DOP_2_COMP_SOLN_SIMPLE(N_OBS,USE,DIRCOS_X1,DIRCOS_X2,VR2,
     $   SD_VR2,BADDATA,VP1_2,SD_VP1_2,VP2_2,SD_VP2_2)

C  Use the three-component solution for each component if available.
C  Otherwise, use the two-component solution.
         IF(VP1_3.NE.BADDATA)THEN
            VP1=VP1_3
            SD_VP1=SD_VP1_3
         ELSE
            VP1=VP1_2
            SD_VP1=SD_VP1_2
         ENDIF
         IF(VP2_3.NE.BADDATA)THEN
            VP2=VP2_3
            SD_VP2=SD_VP2_3
         ELSE
            VP2=VP2_2
            SD_VP2=SD_VP2_2
         ENDIF
         IF(VP3_3.NE.BADDATA)THEN
            VP3=VP3_3
            SD_VP3=SD_VP3_3
         ELSE
            VP3=BADDATA
            SD_VP3=BADDATA
         ENDIF

C  Substitute the assumed values for the third Cartesian component and
C  its standard deviation into the solution if they are missing in the
C  solution and if requested.
         IF(SUBSTITUTE_VP3)THEN
            IF(VP3.EQ.BADDATA.AND.
     $      (VP1.NE.BADDATA.OR.VP2.NE.BADDATA))THEN
               VP3=VP3_ASSUMED
               SD_VP3=SD_VP3_ASSUMED
            ENDIF
         ENDIF

C  If there is not an assumed third Cartesian component, use the
C  three-component solution.
      ELSE
         VP1=VP1_3
         SD_VP1=SD_VP1_3
         VP2=VP2_3
         SD_VP2=SD_VP2_3
         VP3=VP3_3
         SD_VP3=SD_VP3_3
      ENDIF

      END SUBROUTINE DOP_SOLN_DRIVER_3ELSE2_SIMPLE
