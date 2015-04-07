      SUBROUTINE SPIRAL_OUT_2D(IX_START,IY_START,NX,NY,RMAX_IN,IX,IY,
     $DONE)

C  Thomas Matejka NOAA/NSSL 16 February 1995

C  This subroutine advances to the next point on an ever-expanding
C  square spiral on a plane.  All points on the plane are eventually hit
C  out to a maximum specified radius.

C  Input:

C  IX_START is an integer variable that specifies the grid point number
C  in the first dimension of the central point of the spiral.

C  IY_START is an integer variable that specifies the grid point number
C  in the second dimension of the central point of the spiral.

C  NX is an integer variable that specifies the number of grid points on
C  the plane in the first dimension.

C  NY is an integer variable that specifies the number of grid points on
C  the plane in the second dimension.

C  RMAX_IN is an integer variable that specifies the maximum spiral
C  radius, in terms of grid points, to proceed to.  Thus, if RMAX_IN is
C  0, only one point, the central point, is contained in the spiral, and
C  no points can be advanced to.  If RMAX_IN is 1, the 8 points
C  surrounding the central point are also contained in the spiral and
C  these 8 points can be advanced to.  Etc.

C  Input and output:

C  IX is an integer variable that specifies the current grid point
C  number in the first dimension on the spiral and, if and only if DONE
C  is returned .FALSE., returns the next grid point number in the first
C  dimension on the spiral.

C  IY is an integer variable that specifies the current grid point
C  number in the second dimension on the spiral and, if and only if DONE
C  is returned .FALSE., returns the next grid point number in the second
C  dimension on the spiral.  Normally, IX will be IX_START and IY will
C  be IY_START the first time this subroutine is called, so that the
C  march along the spiral starts at its center.

C  Output:

C  DONE is a logical variable that returns .TRUE. if and only if there
C  are no more points along the spiral on the plane and within the
C  specified maximum spiral radius.

      IMPLICIT NONE
      LOGICAL DONE
      INTEGER IX_START,IY_START,NX,NY,IX,IY,JX,JY,R,RMAX,RMAX_IN

C  Transform to coordinates relative to the center of the spiral.
      JX=IX-IX_START
      JY=IY-IY_START

C  Calculate the maximum desired and necessary radius of the spiral.
      RMAX=MIN0(MAX0(IX_START-1,NX-IX_START,IY_START-1,NY-IY_START),
     $RMAX_IN)

C  Calculate the current radius of the spiral.
      R=MAX0(IABS(JX),IABS(JY))

C  Proceed with the next point on the spiral.  Check whether this point
C  is contained within the plane.  If so, then transform the coordinates
C  back and return them.  If not, then proceed with the next point until
C  there are no more left.
      DO
         IF(JX.EQ.0.AND.
     $   JY.EQ.0)THEN
            JY=-1
            R=R+1
         ELSEIF(JX.EQ.1.AND.
     $   JY.EQ.-R)THEN
            JX=0
            JY=JY-1
            R=R+1
         ELSEIF(JX.EQ.-R.AND.
     $   JY.LT.R)THEN
            JY=JY+1
         ELSEIF(JX.EQ.R.AND.
     $   JY.GT.-R)THEN
            JY=JY-1
         ELSEIF(JY.EQ.R)THEN
            JX=JX+1
         ELSE
            JX=JX-1
         ENDIF
         IF(R.GT.RMAX)THEN
            DONE=.TRUE.
            RETURN
         ENDIF
         IX=IX_START+JX
         IY=IY_START+JY
         IF(IX.GE.1.AND.
     $   IX.LE.NX.AND.
     $   IY.GE.1.AND.
     $   IY.LE.NY)THEN
            DONE=.FALSE.
            RETURN
         ENDIF    
      ENDDO
      END
