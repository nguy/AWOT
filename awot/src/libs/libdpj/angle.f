      Function ANGLE(ACAR)

C   Returns THE TRIGONOMETRIC ANGLE (DEG) GIVEN THE CARDINAL ANGLE (DEG)

      A = 450.0 - ACAR
      If (A.gt.360.0) A = A - 360.0
      ANGLE = A
      Return
      End
