HOMDIR	= $(HOME)/software/bin/

#  NCARG compiler

NCARGC	= ncargf90
FC	= gfortran

OPTS	= -ffixed-line-length-132 -fbounds-check

OPTIMIZE = -O2
DEBUG = 

##############################################################################
#
# Home grown libraries
#
##############################################################################

DPJLIBRARIES	= $(HOME)/software/lib
DPJFLIB		= -L$(DPJLIBRARIES)
LIBNAME		= -ldpj #-libf95getdata #-ldpj

##############################################################################
#
# Programs to build
#
##############################################################################

#PROGRAMS        = gage_interp grnd_interp hybrid_interp dualp_interp fast_interp2
PROGRAMS        = fast_interp_so_out
RM              = rm

##############################################################################
#
# Targets to build
#
##############################################################################

all : $(PROGRAMS)

#gage_interp:  gage_interp.f
#	$(NCARGC) $(OPTS) $(OPTIMIZE) $(DEBUG) -o $(HOMDIR)gage_interp gage_interp.f $#(DPJFLIB) $(LIBNAME) -lngmath

#grnd_interp:  grnd_interp.f90
#	$(FC) $(OPTS) $(OPTIMIZE) $(DEBUG) -o $(HOMDIR)grnd_interp grnd_interp.f90 $(D#PJFLIB) $(LIBNAME)

#hybrid_interp:  hybrid_interp.f90
#	$(FC) $(OPTS) $(OPTIMIZE) $(DEBUG) -o $(HOMDIR)hybrid_interp hybrid_interp.f90 $(DPJFLIB) $(LIBNAME)

#dualp_interp:  dualp_interp.f90
#	$(FC) $(OPTS) $(OPTIMIZE) $(DEBUG) -o $(HOMDIR)dualp_interp dualp_interp.f90 $(DPJFLIB) $(LIBNAME)

#fast_interp2:  fast_interp2.f90
#	$(FC) $(OPTS) $(OPTIMIZE) $(DEBUG) -o $(HOMDIR)fast_interp2 fast_interp2.f90 $#(DPJFLIB) $(LIBNAME)

fast_interp_so_out:  fast_interp_so_out.f90
	$(FC) $(OPTS) $(OPTIMIZE) $(DEBUG) -o $(HOMDIR)fast_interp_so_out fast_interp_so_out.f90 $(DPJFLIB) $(LIBNAME)

##############################################################################
clean :
	$(RM) *.o

