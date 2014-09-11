C  tmtlib.inc

C  Thomas Matejka NOAA/NSSL 31 March 1993

C  This INCLUDE file defines constants for TMTLIB.

C  Logical device number for TMTLIB messages.
      INTEGER TMTLIB_MESSAGE_UNIT
      PARAMETER(TMTLIB_MESSAGE_UNIT=7)

C  Value to return when a parameter is undefined.
      REAL TMTLIB_BADFLAG
      PARAMETER(TMTLIB_BADFLAG=-9999.)

C  Temperature tolerance in iterative procedures.
      REAL TTOL
      PARAMETER(TTOL=0.01)

C  Pressure tolerance in iterative procedures.
      REAL PTOL
      PARAMETER(PTOL=5.)

C  Mixing ratio tolerance in iterative procedures.
      REAL WTOL
      PARAMETER(WTOL=1.E-7)

C  Maximum temperature in iterative procedures.
      REAL T_MAX
      PARAMETER(T_MAX=400.)

C  Maximum dew point temperature in iterative procedures.
      REAL TD_MAX
      PARAMETER(TD_MAX=400.)

C  Maximum pressure in iterative procedures.
      REAL P_MAX
      PARAMETER(P_MAX=120000.)

C  Maximum number of iterations in iterative procedures.
      INTEGER MAX_ITERS
      PARAMETER(MAX_ITERS=1000)

C  Temperature increment in incremental procedures.
      REAL DELTA_T
      PARAMETER(DELTA_T=4.)

C  Pressure increment in incremental procedures.
      REAL DELTA_P
      PARAMETER(DELTA_P=1000.)

C  Maximum number of increments in incremental procedures.
      INTEGER MAX_INCS
      PARAMETER(MAX_INCS=1000)
