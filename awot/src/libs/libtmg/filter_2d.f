      SUBROUTINE FILTER_2D(A_IN,SD_A_IN,MAXX_IN,MAXY_IN,NX,NY,BADDATA,
     $UNCORRELATED,WT_X,MAX_REACH,REACH_X,WT_Y,REACH_Y,MAXX_OUT,
     $MAXY_OUT,A_OUT,SD_A_OUT)

C  Thomas Matejka NOAA/NSSL 28 May 1999

C  This subroutine applies a two-dimensional filter to a two-dimensional
C  data field.  Filtered values are produced only at points where the
C  original data field is not missing.  A field of the standard
C  deviation of the filtered data is also produced, either under the
C  assumption that errors in the original data are uncorrelated or under
C  the assumption that errors in the original data are fully correlated.

C  Input:

C  A_IN is a two-dimensional real array.  A_IN(I,J) specifies the value,
C  at the Ith grid point in the first dimension and the Jth grid point
C  in the second dimension, of the data field to be filtered.  If it is
C  missing, it should equal BADDATA.  Errors in the values in A_IN are
C  assumed to be uncorrelated.

C  SD_A_IN is a two-dimensional real array.  SD_A_IN(I,J) specifies the
C  standard deviation of A_IN(I,J).  If it is missing, it should equal
C  BADDATA.  SD_A_IN(I,J) must not be missing for A_IN(I,J) to be
C  usable.

C  MAXX_IN is an integer variable that specifies the first dimension of
C  A_IN and SD_A_IN in the calling program.

C  MAXY_IN is an integer variable that specifies the second dimension of
C  A_IN and SD_A_IN in the calling program.

C  NX is an integer variable that specifies the number of grid points in
C  the first dimension for A_IN, SD_A_IN, A_OUT, and SD_A_OUT.

C  NY is an integer variable that specifies the number of grid points in
C  the second dimension for A_IN, SD_A_IN, A_OUT, and SD_A_OUT.

C  BADDATA is a real variable that indicates a missing value as
C  described.

C  UNCORRELATED is a logical variable.  If UNCORRELATED is .TRUE., then
C  the original data are assumed to be uncorrelated.  If UNCORRELATED is
C  .FALSE., then the original data are assumed to be fully correlated.

C  WT_X is a one-dimensional real array.  WT_X(I) specifies the weight
C  at the Ith grid point relative to the point of application for the
C  one-dimensional filter to be applied in the first dimension.  The
C  weights in WT_X and WT_Y can be specified to within a constant
C  factor, since the filter weights are normalized in the subroutine.

C  MAX_REACH is an integer variable that specifies that WT_X and WT_Y
C  are dimensioned from -MAX_REACH to MAX_REACH in the calling program.

C  REACH_X is an integer variable that specifies that weights in WT_X
C  run from -REACH_X to REACH_X grid points relative to the point of
C  application.

C  WT_Y is a one-dimensional real array.  WT_Y(I) specifies the weight
C  at the Ith grid point relative to the point of application for the
C  one-dimensional filter to be applied in the second dimension.  The
C  weights in WT_X and WT_Y can be specified to within a constant
C  factor, since the filter weights are normalized in the subroutine.

C  REACH_Y is an integer variable that specifies that weights in WT_Y
C  run from -REACH_Y to REACH_Y grid points relative to the point of
C  application.

C  MAXX_OUT is an integer variable that specifies the first dimension of
C  A_OUT and SD_A_OUT in the calling program.

C  MAXY_OUT is an integer variable that specifies the second dimension
C  of A_OUT and SD_A_OUT in the calling program.

C  Output:

C  A_OUT is a two-dimensional real array.  A_OUT(I,J) returns the value,
C  at the Ith grid point in the first dimension and the Jth grid point
C  in the second dimension, of the data field that has been filtered.
C  If it is missing, it is returned as BADDATA.  A_OUT may be identical
C  to A_IN, in which case A_IN will be overwritten.

C  SD_A_OUT is a two-dimensional real array.  SD_A_OUT(I,J) returns the
C  standard deviation of A_OUT(I,J).  If it is missing, it is returned
C  as BADDATA.  SD_A_OUT may be identical to SD_A_IN, in which case
C  SD_A_IN will be overwritten.

      IMPLICIT NONE
      LOGICAL UNCORRELATED
      INTEGER MAXX_IN,MAXY_IN,MAXX_OUT,MAXY_OUT,NX,NY,REACH_X,REACH_Y,
     $IX,IY,JX,JY,MAX_REACH
      REAL BADDATA,SUM,SUM_WT,SUM_VAR,SUM_SD
      REAL WT_X(-MAX_REACH:MAX_REACH),WT_Y(-MAX_REACH:MAX_REACH)
      REAL WT_XY(-REACH_X:REACH_X,-REACH_Y:REACH_Y)
      REAL A_IN(MAXX_IN,MAXY_IN),SD_A_IN(MAXX_IN,MAXY_IN)
      REAL A_OUT(MAXX_OUT,MAXY_OUT),SD_A_OUT(MAXX_OUT,MAXY_OUT)
      REAL A_TEMP(NX,NY),SD_A_TEMP(NX,NY)

C  Calculate the weights for the two-dimensional filter.
      DO JY=-REACH_Y,REACH_Y
         DO JX=-REACH_X,REACH_X
            WT_XY(JX,JY)=WT_X(JX)*WT_Y(JY)
         ENDDO
      ENDDO

C  Copy the input fields in case they are to be overwritten.
      DO IY=1,NY
         DO IX=1,NX
            A_TEMP(IX,IY)=A_IN(IX,IY)
            SD_A_TEMP(IX,IY)=SD_A_IN(IX,IY)
         ENDDO
      ENDDO

C  Loop through the points on the plane.
      DO IY=1,NY
         DO IX=1,NX

C  Check whether a value exists.
            IF(A_TEMP(IX,IY).NE.BADDATA.AND.
     $      SD_A_TEMP(IX,IY).NE.BADDATA)THEN
            
C  Loop through the points within the reach of the filter and calculate
C  sums.
               SUM=0.
               SUM_VAR=0.
               SUM_SD=0.
               SUM_WT=0.
               DO JY=-REACH_Y,REACH_Y
                  IF(IY+JY.GE.1.AND.
     $            IY+JY.LE.NY)THEN
                     DO JX=-REACH_X,REACH_X
                        IF(IX+JX.GE.1.AND.
     $                  IX+JX.LE.NX)THEN
                           IF(A_TEMP(IX+JX,IY+JY).NE.BADDATA.AND.
     $                     SD_A_TEMP(IX+JX,IY+JY).NE.BADDATA)THEN
                              SUM=SUM+A_TEMP(IX+JX,IY+JY)*WT_XY(JX,JY)
                              SUM_VAR=SUM_VAR+(SD_A_TEMP(IX+JX,IY+JY)*
     $                        WT_XY(JX,JY))**2
                              SUM_SD=SUM_SD+SD_A_TEMP(IX+JX,IY+JY)*
     $                        WT_XY(JX,JY)
                              SUM_WT=SUM_WT+WT_XY(JX,JY)
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO

C  Calculate the filtered value and its standard deviation.
               IF(SUM_WT.GE.WT_XY(0,0))THEN
                  A_OUT(IX,IY)=SUM/SUM_WT
                  IF(UNCORRELATED)THEN
                     SD_A_OUT(IX,IY)=SQRT(SUM_VAR)/SUM_WT
                  ELSE
                     SD_A_OUT(IX,IY)=SUM_SD/SUM_WT
                  ENDIF
               ELSE
                  A_OUT(IX,IY)=BADDATA
                  SD_A_OUT(IX,IY)=BADDATA
               ENDIF
            ELSE
               A_OUT(IX,IY)=BADDATA
               SD_A_OUT(IX,IY)=BADDATA
            ENDIF
         ENDDO
      ENDDO

      END SUBROUTINE FILTER_2D
