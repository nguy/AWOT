      SUBROUTINE DOP_1_COMP_SOLN(N_OBS,USE,DIRCOS_X1,VR1,SD_VR1,BADDATA,
     $VP1,SD_VP1)

C  Thomas Matejka NOAA/NSSL 13 September 1995

C  This subroutine calculates Cartesian components of target velocity
C  and their standard deviations from colocated Doppler radar data at
C  one point using a one-component solution.
      
C  Input:

C  N_OBS is an integer variable that specifies the number of colocated
C  velocity observations.

C  USE is a one-dimensional logical array.  USE(I) should be .TRUE. if
C  and only if the Ith velocity observation is to be used in the
C  calculation.

C  DIRCOS_X1 is a one-dimensional real array.  DIRCOS_X1(I) specifies
C  the direction cosine between the radar beam and the first Cartesian
C  coordinate axis for the Ith velocity observation.

C  VR1 is a one-dimensional real array.  VR1(I) specifies the Ith
C  velocity observation reduced by assumed contributions of the second
C  and third Cartesian components of the target velocity.

C  SD_VR1 is a one-dimensional real array.  SD_VR1(I) specifies the
C  standard deviation of VR1(I).

C  BADDATA is a real variable that indicates a missing value as
C  described.

C  Output:

C  VP1 is a real variable that returns the first Cartesian component of
C  target velocity.  If it cannot be calculated, it is returned as
C  BADDATA.

C  SD_VP1 is a real variable that returns the standard deviation of VP1.
C  If it cannot be calculated, it is returned as BADDATA.

      IMPLICIT NONE
      INCLUDE 'tmrlib.inc'
      INTEGER N_OBS
      LOGICAL SUCCESS
      LOGICAL USE(N_OBS)
      INTEGER N_GOOD,I_OBS,EFFECTIVE_NDATA
      REAL VP1,SD_VP1,BADDATA,A11,SSQ
      REAL VAR_COEFF(1),COEFF(1)
      REAL DIRCOS_X1(N_OBS),VR1(N_OBS),SD_VR1(N_OBS),
     $RESPONSE_VECTOR(N_OBS),VAR_RESPONSE_VECTOR(N_OBS),
     $GROUP_SIZE(N_OBS)
      REAL A_INV(1,1)
      REAL PREDICTOR_ARRAY(N_OBS,1)

C  Load the data into arrays for the regression analysis.
      N_GOOD=0
      A11=0.
      IF(N_OBS.GE.1)THEN
         DO I_OBS=1,N_OBS
            IF(USE(I_OBS))THEN
               N_GOOD=N_GOOD+1
               PREDICTOR_ARRAY(N_GOOD,1)=DIRCOS_X1(I_OBS)
               RESPONSE_VECTOR(N_GOOD)=VR1(I_OBS)
               VAR_RESPONSE_VECTOR(N_GOOD)=SD_VR1(I_OBS)**2
               A11=A11+DIRCOS_X1(I_OBS)**2
            ENDIF
         ENDDO
      ENDIF

C  Check whether there are any data to use.
      IF(N_GOOD.GE.1)THEN

C  Check whether the data are degenerate in the first Cartesian 
C  coordinate.
         IF(A11.EQ.0.)THEN
            VP1=BADDATA
            SD_VP1=BADDATA

C  The data are not degenerate.  Perform a one-parameter regression.
         ELSE
            CALL LLS_VAR(1,N_GOOD,PREDICTOR_ARRAY,N_OBS,RESPONSE_VECTOR,
     $      VAR_RESPONSE_VECTOR,.FALSE.,GROUP_SIZE,.FALSE.,1.,
     $      .TRUE.,DOP_SINGULAR_THRESHOLD,MAT_WRITE_MODE,1,
     $      EFFECTIVE_NDATA,COEFF,VAR_COEFF,SSQ,A_INV,SUCCESS)
            IF(SUCCESS)THEN
               VP1=COEFF(1)
               SD_VP1=SQRT(VAR_COEFF(1))

C  The regression was not successful.
            ELSE
               VP1=BADDATA
               SD_VP1=BADDATA
            ENDIF
         ENDIF

C  There are no data to use.
      ELSE
         VP1=BADDATA
         SD_VP1=BADDATA
      ENDIF

C  Done.
      RETURN
      END
