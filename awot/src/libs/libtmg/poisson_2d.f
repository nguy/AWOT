      SUBROUTINE POISSON_2D(F,G,MAXX,MAXY,NX,NY,DELX,DELY,BADDATA,P)

C  Thomas Matejka NOAA/NSSL 17 May 1996

      IMPLICIT NONE
      INTEGER,PARAMETER::MAX_ITERATIONS=1000
      REAL,PARAMETER::CONVERGENCE_FACTOR=0.05
      REAL,PARAMETER::RELAX=1.78
      LOGICAL::CHANGED
      LOGICAL,DIMENSION(6,NX,NY)::CONDITION
      INTEGER::MAXX,MAXY,
     $NX,NY,IX,IY,
     $I,J,ITERATION,IR,IS,N,
     $COUNT
      INTEGER,DIMENSION(NX,NY)::CODE
      REAL::DELX,DELY,
     $BADDATA,
     $A,B,C,D,
     $SUM1,TEMP,
     $CONVERGENCE_THRESHOLD,MAXDIFF
      REAL,DIMENSION(MAXX,MAXY)::F,G,P
      REAL,DIMENSION(NX,NY)::RHS,
     $P_OLD,
     $SUM_ARRAY,
     $DPDX,DPDY
      REAL,DIMENSION(NX,NY,-1:1)::XWEIGHT,YWEIGHT

C  Initially assume that no points are suitable for a solution.
      DO IY=1,NY
         DO IX=1,NX
            CODE(IX,IY)=0
         ENDDO
      ENDDO

C  Identify the points suitable for an interior solution.
      DO IY=1,NY
         DO IX=1,NX
            IF(IX.NE.1.AND.
     $      IX.NE.NX.AND.
     $      IY.NE.1.AND.
     $      IY.NE.NY.AND.
     $      F(IX+1,IY).NE.BADDATA.AND.
     $      F(IX-1,IY).NE.BADDATA.AND.
     $      G(IX,IY+1).NE.BADDATA.AND.
     $      G(IX,IY-1).NE.BADDATA)THEN
               CODE(IX,IY)=1
            ENDIF
         ENDDO
      ENDDO

C  Identify the points suitable for a boundary solution.
      DO IY=1,NY
         DO IX=1,NX
            IF(CODE(IX,IY).NE.1.AND.
     $      (IX.NE.1.AND.CODE(IX-1,IY).EQ.1).OR.
     $      (IX.NE.NX.AND.CODE(IX+1,IY).EQ.1).OR.
     $      (IY.NE.1.AND.CODE(IX,IY-1).EQ.1).OR.
     $      (IY.NE.NY.AND.CODE(IX,IY+1).EQ.1))THEN
               CODE(IX,IY)=2
            ENDIF
         ENDDO
      ENDDO

C  Identify the points suitable for a corner or peninsula solution.
      CHANGED=.TRUE.
      DO
         IF(.NOT.CHANGED)THEN
            EXIT
         ENDIF
         CHANGED=.FALSE.
         DO IY=1,NY
            DO IX=1,NX
               IF(CODE(IX,IY).EQ.0)THEN
                  IF(F(IX,IY).NE.BADDATA)THEN
                     IF((IX.NE.1.AND.CODE(IX-1,IY).NE.0).OR.
     $               (IX.NE.NX.AND.CODE(IX+1,IY).NE.0))THEN
                        CODE(IX,IY)=3
                        CHANGED=.TRUE.
                        CYCLE
                     ENDIF
                  ENDIF
                  IF(G(IX,IY).NE.BADDATA)THEN
                     IF((IY.NE.1.AND.CODE(IX,IY-1).NE.0).OR.
     $               (IY.NE.NY.AND.CODE(IX,IY+1).NE.0))THEN
                        CODE(IX,IY)=3
                        CHANGED=.TRUE.
                        CYCLE
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO

C  Initialize the solution to 0. where a solution is possible and to
C  missing elsewhere.
      N=0
      DO IY=1,NY
         DO IX=1,NX
            IF(CODE(IX,IY).NE.0)THEN
               P(IX,IY)=0.
               N=N+1
            ELSE
               P(IX,IY)=BADDATA
            ENDIF
            P_OLD(IX,IY)=P(IX,IY)
         ENDDO
      ENDDO

C  Check whether any points are suitable for a solution.
      IF(N.EQ.0)THEN
         RETURN
      ENDIF

C  Calculate the maximum change threshold.
      SUM1=0.
      COUNT=0
      DO IY=1,NY
         DO IX=1,NX
            IF(F(IX,IY).NE.BADDATA)THEN
               IF(G(IX,IY).NE.BADDATA)THEN
                  COUNT=COUNT+1
                  SUM1=SUM1+SQRT((F(IX,IY)*DELX)**2+(G(IX,IY)*DELY)**2)
               ELSE
                  COUNT=COUNT+1
                  SUM1=SUM1+ABS(F(IX,IY)*DELX)
               ENDIF
            ELSE
               IF(G(IX,IY).NE.BADDATA)THEN
                  COUNT=COUNT+1
                  SUM1=SUM1+ABS(G(IX,IY)*DELY)
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      IF(COUNT.NE.0)THEN
         CONVERGENCE_THRESHOLD=CONVERGENCE_FACTOR*SUM1/FLOAT(COUNT)
      ENDIF

C  Set the weights.
      DO IR=-1,1
         DO IY=1,NY
            DO IX=1,NX
               XWEIGHT(IX,IY,IR)=0.
               YWEIGHT(IX,IY,IR)=0.
            ENDDO
         ENDDO
      ENDDO
      DO IY=1,NY
         DO IX=1,NX
            IF(CODE(IX,IY).NE.0)THEN
               IF(IX.NE.1.AND.
     $         IX.NE.NX.AND.
     $         CODE(IX-1,IY).NE.0.AND.
     $         CODE(IX+1,IY).NE.0)THEN
                  XWEIGHT(IX,IY,-1)=-1./(2.*DELX)
                  XWEIGHT(IX,IY,1)=1./(2.*DELX)
               ELSEIF(IX.NE.1.AND.
     $         CODE(IX-1,IY).NE.0)THEN
                  XWEIGHT(IX,IY,-1)=-1./DELX
                  XWEIGHT(IX,IY,0)=1./DELX
               ELSEIF(IX.NE.NX.AND.
     $         CODE(IX+1,IY).NE.0)THEN
                  XWEIGHT(IX,IY,0)=-1./DELX
                  XWEIGHT(IX,IY,1)=1./DELX
               ENDIF
               IF(IY.NE.1.AND.
     $         IY.NE.NY.AND.
     $         CODE(IX,IY-1).NE.0.AND.
     $         CODE(IX,IY+1).NE.0)THEN
                  YWEIGHT(IX,IY,-1)=-1./(2.*DELY)
                  YWEIGHT(IX,IY,1)=1./(2.*DELY)
               ELSEIF(IY.NE.1.AND.
     $         CODE(IX,IY-1).NE.0)THEN
                  YWEIGHT(IX,IY,-1)=-1./DELY
                  YWEIGHT(IX,IY,0)=1./DELY
               ELSEIF(IY.NE.NY.AND.
     $         CODE(IX,IY+1).NE.0)THEN
                  YWEIGHT(IX,IY,0)=-1./DELY
                  YWEIGHT(IX,IY,1)=1./DELY
               ENDIF
            ENDIF
         ENDDO
      ENDDO

C  Calculate the right hand side.
      DO IY=1,NY
         DO IX=1,NX
            IF(CODE(IX,IY).EQ.1)THEN
               RHS(IX,IY)=XWEIGHT(IX-1,IY,1)*F(IX-1,IY)+
     $         XWEIGHT(IX+1,IY,-1)*F(IX+1,IY)+
     $         YWEIGHT(IX,IY-1,1)*G(IX,IY-1)+
     $         YWEIGHT(IX,IY+1,-1)*G(IX,IY+1)
            ELSE
               RHS(IX,IY)=0.
            ENDIF
         ENDDO
      ENDDO

C  Perform calculations that will not vary between iterations.
      DO I=1,6
         DO IY=1,NY
            DO IX=1,NX
               CONDITION(I,IX,IY)=.FALSE.
            ENDDO
         ENDDO
      ENDDO
      DO IY=1,NY
         DO IX=1,NX
            IF(CODE(IX,IY).EQ.2.OR.
     $      CODE(IX,IY).EQ.3)THEN
               SUM_ARRAY(IX,IY)=0.
               IF(F(IX,IY).NE.BADDATA.AND.
     $         XWEIGHT(IX,IY,0).NE.0.)THEN
                  CONDITION(1,IX,IY)=.TRUE.
                  SUM_ARRAY(IX,IY)=SUM_ARRAY(IX,IY)+
     $            ABS(XWEIGHT(IX,IY,0))
               ENDIF
               IF(G(IX,IY).NE.BADDATA.AND.
     $         YWEIGHT(IX,IY,0).NE.0.)THEN
                  CONDITION(2,IX,IY)=.TRUE.
                  SUM_ARRAY(IX,IY)=SUM_ARRAY(IX,IY)+
     $            ABS(YWEIGHT(IX,IY,0))
               ENDIF
               IF(IY.NE.NY.AND.
     $         G(IX,IY+1).NE.BADDATA.AND.
     $         YWEIGHT(IX,IY+1,-1).NE.0.)THEN
                  CONDITION(3,IX,IY)=.TRUE.
                  SUM_ARRAY(IX,IY)=SUM_ARRAY(IX,IY)+
     $            ABS(YWEIGHT(IX,IY+1,-1))
               ENDIF
               IF(IX.NE.1.AND.
     $         F(IX-1,IY).NE.BADDATA.AND.
     $         XWEIGHT(IX-1,IY,1).NE.0.)THEN
                  CONDITION(4,IX,IY)=.TRUE.
                  SUM_ARRAY(IX,IY)=SUM_ARRAY(IX,IY)+
     $            ABS(XWEIGHT(IX-1,IY,1))
               ENDIF
               IF(IX.NE.NX.AND.
     $         F(IX+1,IY).NE.BADDATA.AND.
     $         XWEIGHT(IX+1,IY,-1).NE.0.)THEN
                  CONDITION(5,IX,IY)=.TRUE.
                  SUM_ARRAY(IX,IY)=SUM_ARRAY(IX,IY)+
     $            ABS(XWEIGHT(IX+1,IY,-1))
               ENDIF
               IF(IY.NE.1.AND.
     $         G(IX,IY-1).NE.BADDATA.AND.
     $         YWEIGHT(IX,IY-1,1).NE.0.)THEN
                  CONDITION(6,IX,IY)=.TRUE.
                  SUM_ARRAY(IX,IY)=SUM_ARRAY(IX,IY)+
     $            ABS(YWEIGHT(IX,IY-1,1))
               ENDIF
            ENDIF
         ENDDO
      ENDDO

C  Iterate to a solution.
      ITERATION=0
      DO
         ITERATION=ITERATION+1
         DO IY=1,NY
            DO IX=1,NX
               IF(CODE(IX,IY).EQ.1)THEN
                  A=0.
                  B=0.
                  C=0.
                  D=0.
                  DO IR=-1,1,1
                     DO I=IX-1,IX+1,1
                        IF(I+IR.EQ.IX)THEN
                           A=A+XWEIGHT(I,IY,IR)*XWEIGHT(I,IY,IX-I)
                        ELSEIF(I+IR.GE.1.AND.
     $                  I+IR.LE.NX)THEN
                           B=B+P(I+IR,IY)*XWEIGHT(I,IY,IR)*
     $                     XWEIGHT(I,IY,IX-I)
                        ENDIF
                     ENDDO
                  ENDDO
                  DO IS=-1,1,1
                     DO J=IY-1,IY+1,1
                        IF(J+IS.EQ.IY)THEN
                           C=C+YWEIGHT(IX,J,IS)*YWEIGHT(IX,J,IY-J)
                        ELSEIF(J+IS.GE.1.AND.
     $                  J+IS.LE.NY)THEN
                           D=D+P(IX,J+IS)*YWEIGHT(IX,J,IS)*
     $                     YWEIGHT(IX,J,IY-J)
                        ENDIF
                     ENDDO
                  ENDDO
                  P(IX,IY)=P(IX,IY)+
     $            (((RHS(IX,IY)-(B+D))/(A+C))-P(IX,IY))*RELAX
               ELSEIF(CODE(IX,IY).EQ.2.OR.
     $         CODE(IX,IY).EQ.3)THEN
                  SUM1=0.
                  IF(CONDITION(1,IX,IY))THEN
                     TEMP=F(IX,IY)
                     IF(IX.NE.1.AND.
     $               P(IX-1,IY).NE.BADDATA)THEN
                        TEMP=TEMP-P(IX-1,IY)*XWEIGHT(IX,IY,-1)
                     ENDIF
                     IF(IX.NE.NX.AND.
     $               P(IX+1,IY).NE.BADDATA)THEN
                        TEMP=TEMP-P(IX+1,IY)*XWEIGHT(IX,IY,1)
                     ENDIF
                     SUM1=SUM1+TEMP*SIGN(1.,XWEIGHT(IX,IY,0))
                  ENDIF
                  IF(CONDITION(2,IX,IY))THEN
                     TEMP=G(IX,IY)
                     IF(IY.NE.1.AND.
     $               P(IX,IY-1).NE.BADDATA)THEN
                        TEMP=TEMP-P(IX,IY-1)*YWEIGHT(IX,IY,-1)
                     ENDIF
                     IF(IY.NE.NY.AND.
     $               P(IX,IY+1).NE.BADDATA)THEN
                        TEMP=TEMP-P(IX,IY+1)*YWEIGHT(IX,IY,1)
                     ENDIF
                     SUM1=SUM1+TEMP*SIGN(1.,YWEIGHT(IX,IY,0))
                  ENDIF
                  IF(CONDITION(3,IX,IY))THEN
                     TEMP=G(IX,IY+1)
                     IF(P(IX,IY+1).NE.BADDATA)THEN
                        TEMP=TEMP-P(IX,IY+1)*YWEIGHT(IX,IY+1,0)
                     ENDIF
                     IF(IY+2.LE.NY.AND.
     $               P(IX,IY+2).NE.BADDATA)THEN
                        TEMP=TEMP-P(IX,IY+2)*YWEIGHT(IX,IY+1,1)
                     ENDIF
                     SUM1=SUM1+TEMP*SIGN(1.,YWEIGHT(IX,IY+1,-1))
                  ENDIF
                  IF(CONDITION(4,IX,IY))THEN
                     TEMP=F(IX-1,IY)
                     IF(IX-2.GE.1.AND.
     $               P(IX-2,IY).NE.BADDATA)THEN
                        TEMP=TEMP-P(IX-2,IY)*XWEIGHT(IX-1,IY,-1)
                     ENDIF
                     IF(P(IX-1,IY).NE.BADDATA)THEN
                        TEMP=TEMP-P(IX-1,IY)*XWEIGHT(IX-1,IY,0)
                     ENDIF
                     SUM1=SUM1+TEMP*SIGN(1.,XWEIGHT(IX-1,IY,1))
                  ENDIF
                  IF(CONDITION(5,IX,IY))THEN
                     TEMP=F(IX+1,IY)
                     IF(P(IX+1,IY).NE.BADDATA)THEN
                        TEMP=TEMP-P(IX+1,IY)*XWEIGHT(IX+1,IY,0)
                     ENDIF
                     IF(IX+2.LE.NX.AND.
     $               P(IX+2,IY).NE.BADDATA)THEN
                        TEMP=TEMP-P(IX+2,IY)*XWEIGHT(IX+1,IY,1)
                     ENDIF
                     SUM1=SUM1+TEMP*SIGN(1.,XWEIGHT(IX+1,IY,-1))
                  ENDIF
                  IF(CONDITION(6,IX,IY))THEN
                     TEMP=G(IX,IY-1)
                     IF(IY-2.GE.1.AND.
     $               P(IX,IY-2).NE.BADDATA)THEN
                        TEMP=TEMP-P(IX,IY-2)*YWEIGHT(IX,IY-1,-1)
                     ENDIF
                     IF(IY-1.GE.1.AND.
     $               P(IX,IY-1).NE.BADDATA)THEN
                        TEMP=TEMP-P(IX,IY-1)*YWEIGHT(IX,IY-1,0)
                     ENDIF
                     SUM1=SUM1+TEMP*SIGN(1.,YWEIGHT(IX,IY-1,1))
                  ENDIF
                  IF(SUM_ARRAY(IX,IY).NE.0.)THEN
                     P(IX,IY)=P(IX,IY)+
     $               (SUM1/SUM_ARRAY(IX,IY)-P(IX,IY))*RELAX
                  ELSE
                     SUM1=0.
                     COUNT=0
                     IF(IX.GT.1.AND.
     $               CODE(IX-1,IY).NE.0)THEN
                        SUM1=SUM1+P(IX-1,IY)
                        COUNT=COUNT+1
                     ENDIF
                     IF(IX.LT.NX.AND.
     $               CODE(IX+1,IY).NE.0)THEN
                        SUM1=SUM1+P(IX+1,IY)
                        COUNT=COUNT+1
                     ENDIF
                     IF(IY.GT.1.AND.
     $               CODE(IX,IY-1).NE.0)THEN
                        SUM1=SUM1+P(IX,IY-1)
                        COUNT=COUNT+1
                     ENDIF
                     IF(IY.LT.NY.AND.
     $               CODE(IX,IY+1).NE.0)THEN
                        SUM1=SUM1+P(IX,IY+1)
                        COUNT=COUNT+1
                     ENDIF
                     P(IX,IY)=SUM1/FLOAT(COUNT)
                  ENDIF
               ENDIF
            ENDDO
         ENDDO

C  Subtract the mean from the solution.
         CALL SUBTRACT_MEAN_3D(P,MAXX,MAXY,1,NX,NY,1,BADDATA,P)

C  Check whether the procedure has converged.
         IF(MOD(ITERATION,10).EQ.0)THEN
            MAXDIFF=-1.0
            DO IY=1,NY
               DO IX=1,NX
                  IF(P(IX,IY).NE.BADDATA)THEN
                     TEMP=ABS(P_OLD(IX,IY)-P(IX,IY))
                     IF(TEMP.GT.MAXDIFF)THEN
                        MAXDIFF=TEMP
                     ENDIF
                  ENDIF
                  P_OLD(IX,IY)=P(IX,IY)
               ENDDO
            ENDDO
            IF(ITERATION.GE.MAX_ITERATIONS)THEN
               EXIT
            ENDIF
            IF(MAXDIFF.LE.CONVERGENCE_THRESHOLD)THEN
               EXIT
            ENDIF
         ENDIF

C  Return for the next iteration.
      ENDDO

C  Done.
      END SUBROUTINE POISSON_2D
