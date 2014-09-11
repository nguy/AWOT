      SUBROUTINE DOP_2_COMP_SOLN_SIMPLE(N_OBS,USE,DIRCOS_X1,DIRCOS_X2,
     $VR2,SD_VR2,BADDATA,VP1,SD_VP1,VP2,SD_VP2)

C  Thomas Matejka NOAA/NSSL 24 August 2001

C  This subroutine calculates Cartesian components of target velocity
C  and their standard deviations from colocated Doppler radar data at
C  one point using a two-component solution.

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

C  VR2 is a one-dimensional real array.  VR2(I) specifies the Ith
C  velocity observation reduced by an assumed contribution of the third
C  Cartesian component of the target velocity.

C  SD_VR2 is a one-dimensional real array.  SD_VR2(I) specifies the
C  standard deviation of VR2(I).

C  BADDATA is a real variable that indicates a missing value as
C  described.

C  Output:

C  VP1 is a real variable that returns the first Cartesian component of
C  target velocity.  If it cannot be calculated, it is returned as
C  BADDATA.

C  SD_VP1 is a real variable that returns the standard deviation of VP1.
C  If it cannot be calculated, it is returned as BADDATA.

C  VP2 is a real variable that returns the second Cartesian component of
C  target velocity.  If it cannot be calculated, it is returned as
C  BADDATA.

C  SD_VP2 is a real variable that returns the standard deviation of VP2.
C  If it cannot be calculated, it is returned as BADDATA.

      IMPLICIT NONE
      INCLUDE 'tmrlib.inc'
      INTEGER::N_OBS
      LOGICAL::SUCCESS
      LOGICAL,DIMENSION(N_OBS)::USE
      INTEGER::N_GOOD,I_OBS,EFFECTIVE_NDATA
      REAL::VP1,VP2,SD_VP1,SD_VP2,BADDATA,A11,A22,SSQ
      REAL,DIMENSION(2)::VAR_COEFF,COEFF
      REAL,DIMENSION(N_OBS)::DIRCOS_X1,DIRCOS_X2,VR2,SD_VR2,
     $RESPONSE_VECTOR,VAR_RESPONSE_VECTOR,GROUP_SIZE
      REAL,DIMENSION(2,2)::A_INV
      REAL,DIMENSION(N_OBS,2)::PREDICTOR_ARRAY

      N_GOOD=0
      A11=0.
      A22=0.
      IF(N_OBS.GE.1)THEN
         DO I_OBS=1,N_OBS
            IF(USE(I_OBS))THEN
               N_GOOD=N_GOOD+1
               PREDICTOR_ARRAY(N_GOOD,1)=DIRCOS_X1(I_OBS)
               PREDICTOR_ARRAY(N_GOOD,2)=DIRCOS_X2(I_OBS)
               RESPONSE_VECTOR(N_GOOD)=VR2(I_OBS)
               VAR_RESPONSE_VECTOR(N_GOOD)=SD_VR2(I_OBS)**2
               A11=A11+DIRCOS_X1(I_OBS)**2
               A22=A22+DIRCOS_X2(I_OBS)**2
            ENDIF
         ENDDO
      ENDIF

      IF(N_GOOD.GE.2)THEN
         CALL LLS_VAR(2,N_GOOD,PREDICTOR_ARRAY,N_OBS,RESPONSE_VECTOR,
     $   VAR_RESPONSE_VECTOR,.FALSE.,GROUP_SIZE,.FALSE.,1.,
     $   .TRUE.,DOP_SINGULAR_THRESHOLD,MAT_WRITE_MODE,2,
     $   EFFECTIVE_NDATA,COEFF,VAR_COEFF,SSQ,A_INV,SUCCESS)
         IF(SUCCESS)THEN
            VP1=COEFF(1)
            SD_VP1=SQRT(VAR_COEFF(1))
            VP2=COEFF(2)
            SD_VP2=SQRT(VAR_COEFF(2))
         ELSE
            VP1=BADDATA
            SD_VP1=BADDATA
            VP2=BADDATA
            SD_VP2=BADDATA
         ENDIF
      ELSE
         VP1=BADDATA
         SD_VP1=BADDATA
         VP2=BADDATA
         SD_VP2=BADDATA
      ENDIF

      END SUBROUTINE DOP_2_COMP_SOLN_SIMPLE
