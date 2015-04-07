C  tmmlib.inc

C  Thomas Matejka NOAA/NSSL 25 July 2002

C  This INCLUDE file defines constants and memory sizes for TMMLIB.

C  Logical device number for TMMLIB messages.
      INTEGER TMMLIB_MESSAGE_UNIT
      PARAMETER(TMMLIB_MESSAGE_UNIT=7)

C  Maximum number of iterations to allow when adjusting wind component
C  fields to a divergence field.
      INTEGER ADJUST_WINDS_TO_DIV_MAX_ITERATIONS
      PARAMETER(ADJUST_WINDS_TO_DIV_MAX_ITERATIONS=200)

C  Maximum amount of adjustment of a wind component to allow when
C  terminating the iterative adjustment of the wind component fields to
C  a divergence field.
      REAL ADJUST_WINDS_TO_DIV_THRESH
      PARAMETER(ADJUST_WINDS_TO_DIV_THRESH=0.05)

C  Threshold for determining nearly singular matrices for the HOLE_FILL*
C  family of routines.
      REAL HOLE_FILL_SINGULAR_THRESHOLD
      PARAMETER(HOLE_FILL_SINGULAR_THRESHOLD=0.001)

C  Threshold for determining whether two vectors are distinct enough to
C  determine a plane.
      REAL PLANE_THRESHOLD
      PARAMETER(PLANE_THRESHOLD=0.0001)

C  Maximum number of significant digits for a real number.
      INTEGER MAX_SIGNIF
      PARAMETER(MAX_SIGNIF=8)

C  Maximum character string length.
      INTEGER MAX_STRING
      PARAMETER(MAX_STRING=200)

C  Maximum character string length for numbers.
      INTEGER MAX_NUMBER_STRING
      PARAMETER(MAX_NUMBER_STRING=25)
