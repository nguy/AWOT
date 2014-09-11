      SUBROUTINE DOP_3_COMP_SOLN(N_OBS,USE,DIRCOS_X1,DIRCOS_X2,
     $DIRCOS_X3,VR,SD_VR,BADDATA,VP1,SD_VP1,VP2,SD_VP2,VP3,SD_VP3)

C  Thomas Matejka NOAA/NSSL 13 September 1995

C  This subroutine calculates Cartesian components of target velocity
C  and their standard deviations from colocated Doppler radar data at
C  one point using a three-component solution.

C  A two-component solution is used if the data are degenerate in one
C  coordinate.

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
C  target velocity.  If it cannot be calculated, it is returned as
C  BADDATA.

C  SD_VP1 is a real variable that returns the standard deviation of VP1.
C  If it cannot be calculated, it is returned as BADDATA.

C  VP2 is a real variable that returns the second Cartesian component of
C  target velocity.  If it cannot be calculated, it is returned as
C  BADDATA.

C  SD_VP2 is a real variable that returns the standard deviation of VP2.
C  If it cannot be calculated, it is returned as BADDATA.

C  VP3 is a real variable that returns the third Cartesian component of
C  target velocity.  If it cannot be calculated, it is returned as
C  BADDATA.

C  SD_VP3 is a real variable that returns the standard deviation of VP3.
C  If it cannot be calculated, it is returned as BADDATA.

      IMPLICIT NONE
      INCLUDE 'tmrlib.inc'
      INTEGER N_OBS
      LOGICAL USE(N_OBS),SUCCESS
      INTEGER N_GOOD,I_OBS,EFFECTIVE_NDATA
      REAL VP1,VP2,VP3,SD_VP1,SD_VP2,SD_VP3,BADDATA,A11,A22,A33,SSQ
      REAL VAR_COEFF(3),COEFF(3)
      REAL DIRCOS_X1(N_OBS),DIRCOS_X2(N_OBS),DIRCOS_X3(N_OBS),VR(N_OBS),
     $SD_VR(N_OBS),RESPONSE_VECTOR(N_OBS),VAR_RESPONSE_VECTOR(N_OBS),
     $GROUP_SIZE(N_OBS)
      REAL A_INV(3,3)
      REAL PREDICTOR_ARRAY(N_OBS,3)

C  Load the data into arrays for the regression analysis.
      N_GOOD=0
      A11=0.
      A22=0.
      A33=0.
      IF(N_OBS.GE.1)THEN
         DO I_OBS=1,N_OBS
            IF(USE(I_OBS))THEN
               N_GOOD=N_GOOD+1
               PREDICTOR_ARRAY(N_GOOD,1)=DIRCOS_X1(I_OBS)
               PREDICTOR_ARRAY(N_GOOD,2)=DIRCOS_X2(I_OBS)
               PREDICTOR_ARRAY(N_GOOD,3)=DIRCOS_X3(I_OBS)
               RESPONSE_VECTOR(N_GOOD)=VR(I_OBS)
               VAR_RESPONSE_VECTOR(N_GOOD)=SD_VR(I_OBS)**2
               A11=A11+DIRCOS_X1(I_OBS)**2
               A22=A22+DIRCOS_X2(I_OBS)**2
               A33=A33+DIRCOS_X3(I_OBS)**2
            ENDIF
         ENDDO
      ENDIF

C  Check whether there are any data to use.
      IF(N_GOOD.GE.1)THEN

C  Check whether the data are degenerate in the third Cartesian
C  coordinate.  If they are, perform a two-component solution.
         IF(A33.EQ.0.)THEN
            CALL DOP_2_COMP_SOLN(N_OBS,USE,DIRCOS_X1,DIRCOS_X2,VR,SD_VR,
     $      BADDATA,VP1,SD_VP1,VP2,SD_VP2)
            VP3=BADDATA
            SD_VP3=BADDATA

C  Check whether the data are degenerate in the second Cartesian
C  coordinate.  If they are, perform a two-component solution.
         ELSEIF(A22.EQ.0.)THEN
            CALL DOP_2_COMP_SOLN(N_OBS,USE,DIRCOS_X1,DIRCOS_X3,VR,SD_VR,
     $      BADDATA,VP1,SD_VP1,VP3,SD_VP3)
            VP2=BADDATA
            SD_VP2=BADDATA

C  Check whether the data are degenerate in the first Cartesian
C  coordinate.  If they are, perform a two-component solution.
         ELSEIF(A11.EQ.0.)THEN
            CALL DOP_2_COMP_SOLN(N_OBS,USE,DIRCOS_X2,DIRCOS_X3,VR,SD_VR,
     $      BADDATA,VP2,SD_VP2,VP3,SD_VP3)
            VP1=BADDATA
            SD_VP1=BADDATA

C  The data are not degenerate.  Perform a three-parameter regression.
         ELSE
            CALL LLS_VAR(3,N_GOOD,PREDICTOR_ARRAY,N_OBS,RESPONSE_VECTOR,
     $      VAR_RESPONSE_VECTOR,.FALSE.,GROUP_SIZE,.FALSE.,1.,
     $      .TRUE.,DOP_SINGULAR_THRESHOLD,MAT_WRITE_MODE,3,
     $      EFFECTIVE_NDATA,COEFF,VAR_COEFF,SSQ,A_INV,SUCCESS)
            IF(SUCCESS)THEN
               VP1=COEFF(1)
               SD_VP1=SQRT(VAR_COEFF(1))
               VP2=COEFF(2)
               SD_VP2=SQRT(VAR_COEFF(2))
               VP3=COEFF(3)
               SD_VP3=SQRT(VAR_COEFF(3))

C  The regression was not successful.
            ELSE
               VP1=BADDATA
               SD_VP1=BADDATA
               VP2=BADDATA
               SD_VP2=BADDATA
               VP3=BADDATA
               SD_VP3=BADDATA
            ENDIF
         ENDIF

C  There are no data to use.
      ELSE
         VP1=BADDATA
         SD_VP1=BADDATA
         VP2=BADDATA
         SD_VP2=BADDATA
         VP3=BADDATA
         SD_VP3=BADDATA
      ENDIF

C  Done.
      RETURN
      END
