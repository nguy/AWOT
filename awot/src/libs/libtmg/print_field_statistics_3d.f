      SUBROUTINE PRINT_FIELD_STATISTICS_3D(A,MAXX,MAXY,MAXZ,IX_START,
     $IX_END,IY_START,IY_END,IZ_START,IZ_END,BADDATA,LU,LABEL)

C  Thomas Matejka NOAA/NSSL 18 November 1996

C  This subroutine writes formatted statistics for each plane in the
C  first and second dimensions and for the entire volume of a
C  rectangular subset of a three-dimensional data field.  No carriage
C  controls are specified.  122 columns are required.

C  Input:

C  A is a three-dimensional real array.  A(I,J,K) specifies the value,
C  at the Ith grid point in the first dimension, the Jth grid point in
C  the second dimension, and the Kth grid point in the third dimension,
C  of the data field whose statistics are sought.  If it is missing, it
C  should equal BADDATA.

C  MAXX is an integer variable that specifies the first dimension of A
C  in the calling program.

C  MAXY is an integer variable that specifies the second dimension of A
C  in the calling program.

C  MAXZ is an integer variable that specifies the third dimension of A
C  in the calling program.

C  IX_START is an integer variable that specifies the first grid point
C  number in the first dimension of the rectangular subset of A for
C  which statistics will be compiled.

C  IX_END is an integer variable that specifies the last grid point
C  number in the first dimension of the rectangular subset of A for
C  which statistics will be compiled.

C  IY_START is an integer variable that specifies the first grid point
C  number in the second dimension of the rectangular subset of A for
C  which statistics will be compiled.

C  IY_END is an integer variable that specifies the last grid point
C  number in the second dimension of the rectangular subset of A for
C  which statistics will be compiled.

C  IZ_START is an integer variable that specifies the first grid point
C  number in the third dimension of the rectangular subset of A for
C  which statistics will be compiled.

C  IZ_END is an integer variable that specifies the last grid point
C  number in the third dimension of the rectangular subset of A for
C  which statistics will be compiled.

C  BADDATA is a real variable that indicates a missing value as
C  described.

C  LU is an integer variable that specifies the logical unit number to
C  which to write the statistics.

C  LABEL is an identifying character string that will be printed at the
C  top of the table of statistics.

      IMPLICIT NONE

      INTEGER S_L

      CHARACTER*(*) LABEL
      INTEGER MAXX,MAXY,MAXZ,IX_START,IX_END,IY_START,IY_END,IZ_START,
     $IZ_END,IZ,N_BAD,N_GOOD,IX_A_MIN,IY_A_MIN,IZ_A_MIN,IX_A_MAX,
     $IY_A_MAX,IZ_A_MAX,LU,IX_MIN,IX_MAX,IY_MIN,IY_MAX,IZ_MIN,IZ_MAX
      REAL BADDATA,MEAN,SD,A_MIN,A_MAX
      REAL A(MAXX,MAXY,MAXZ)

C  Write table headers.
      WRITE(LU,"(/124('=')/
     $A//
     $'starting (ix,iy,iz) = ','(',I3,',',I3,',',I3,')'/
     $'ending (ix,iy,iz)   = ','(',I3,',',I3,',',I3,')'//
     $3X,'iz',4X,'good',5X,'bad',6X,'(ix,iy)-(ix,iy)',9X,'mean',11X,
     $'sd',10X,'min',5X,'(ix,iy,iz)',10X,'max',5X,'(ix,iy,iz)'/)")
     $LABEL(1:S_L(LABEL)),
     $IX_START,IY_START,IZ_START,
     $IX_END,IY_END,IZ_END

C  Calculate and write statistics for each plane in the first and second
C  dimensions of the rectangular subset of the three-dimensional data
C  field.
      DO IZ=IZ_END,IZ_START,-1
         CALL FIELD_STATISTICS_3D(A,MAXX,MAXY,MAXZ,IX_START,IX_END,
     $ IY_START,IY_END,IZ,IZ,BADDATA,N_GOOD,N_BAD,MEAN,SD,A_MIN,
     $ IX_A_MIN,IY_A_MIN,IZ_A_MIN,A_MAX,IX_A_MAX,IY_A_MAX,IZ_A_MAX,
     $ IX_MIN,IX_MAX,IY_MIN,IY_MAX,IZ_MIN,IZ_MAX)
         WRITE(LU,"(I5,I8,I8,2X,'(',I3,',',I3,')-(',I3,',',I3,')',E13.5,
     $ E13.5,E13.5,2X,'(',I3,',',I3,',',I3,')',E13.5,2X,'(',I3,',',I3,
     $ ',',I3,')')")
     $ IZ,N_GOOD,N_BAD,IX_MIN,IY_MIN,IX_MAX,IY_MAX,MEAN,SD,A_MIN,
     $ IX_A_MIN,IY_A_MIN,IZ_A_MIN,A_MAX,IX_A_MAX,IY_A_MAX,IZ_A_MAX
      ENDDO

C  Calculate and write statistics for the entire volume of the
C  rectangular subset of the three-dimensional data field.
      CALL FIELD_STATISTICS_3D(A,MAXX,MAXY,MAXZ,IX_START,IX_END,
     $IY_START,IY_END,IZ_START,IZ_END,BADDATA,N_GOOD,N_BAD,MEAN,SD,
     $A_MIN,IX_A_MIN,IY_A_MIN,IZ_A_MIN,A_MAX,IX_A_MAX,IY_A_MAX,IZ_A_MAX,
     $IX_MIN,IX_MAX,IY_MIN,IY_MAX,IZ_MIN,IZ_MAX)
      WRITE(LU,"(/'total',I8,I8,2X,'(',I3,',',I3,')-(',I3,',',I3,')',
     $E13.5,E13.5,E13.5,2X,'(',I3,',',I3,',',I3,')',E13.5,2X,'(',I3,',',
     $I3,',',I3,')'/
     $124('=')/)")
     $N_GOOD,N_BAD,IX_MIN,IY_MIN,IX_MAX,IY_MAX,MEAN,SD,A_MIN,IX_A_MIN,
     $IY_A_MIN,IZ_A_MIN,A_MAX,IX_A_MAX,IY_A_MAX,IZ_A_MAX

C  Done.
      END
