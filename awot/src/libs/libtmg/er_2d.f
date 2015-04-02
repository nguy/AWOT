      SUBROUTINE ER_2D(P,F,G,MAXX_IN,MAXY_IN,NX,NY,DELX,DELY,BADDATA,ER)

C  Thomas Matejka NOAA/NSSL 17 May 1996

      IMPLICIT NONE
      INTEGER::MAXX_IN,MAXY_IN,
     $NX,NY,IX,IY,
     $COUNT
      REAL::DELX,DELY,
     $BADDATA,
     $ER,
     $SUM1,SUM2
      REAL,DIMENSION(MAXX_IN,MAXY_IN)::P,F,G
      REAL,DIMENSION(NX,NY)::DDX_P,DDY_P

      CALL DDX_3D(P,MAXX_IN,MAXY_IN,1,NX,NY,1,DELX,BADDATA,NX,NY,1,
     $DDX_P)
      CALL DDY_3D(P,MAXX_IN,MAXY_IN,1,NX,NY,1,DELY,BADDATA,NX,NY,1,
     $DDY_P)
      SUM1=0.
      SUM2=0.
      COUNT=0
      DO IY=1,NY
         DO IX=1,NX
            IF(DDX_P(IX,IY).NE.BADDATA.AND.
     $      F(IX,IY).NE.BADDATA)THEN
               IF(DDY_P(IX,IY).NE.BADDATA.AND.
     $         G(IX,IY).NE.BADDATA)THEN
                  SUM1=SUM1+(DDX_P(IX,IY)-F(IX,IY))**2+
     $            (DDY_P(IX,IY)-G(IX,IY))**2 
                  SUM2=SUM2+(F(IX,IY)**2+G(IX,IY)**2)
                  COUNT=COUNT+2
               ELSE
                  SUM1=SUM1+(DDX_P(IX,IY)-F(IX,IY))**2 
                  SUM2=SUM2+F(IX,IY)**2
                  COUNT=COUNT+1
               ENDIF
            ELSE
               IF(DDY_P(IX,IY).NE.BADDATA.AND.
     $         G(IX,IY).NE.BADDATA)THEN
                  SUM1=SUM1+(DDY_P(IX,IY)-G(IX,IY))**2
                  SUM2=SUM2+G(IX,IY)**2
                  COUNT=COUNT+1
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      IF(SUM2.NE.0.)THEN
         ER=SUM1/SUM2
      ELSE
         ER=BADDATA
      ENDIF
      END
