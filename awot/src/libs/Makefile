##############################################################################
#
# Makefile for archive libraries
#
##############################################################################

OPT       = -O
INSTALLDIR = /Users/nickguy/software/lib
##############################################################################
#
# Defines, Include directories, and flags
#
##############################################################################

DEFINES 	=  
INCLUDES	=  

##############################################################################
#
# Compilers
#
##############################################################################

FC          = gfortran
FFLAGS      = -ffixed-line-length-132

CC          = gcc -w
CFLAGS      = -I/Users/nickguy/software/p3_analysis/davej_radar/libs/librah_C -I/usr/local/include -I/usr/X11/include -I/usr/OpenMotif/include

##############################################################################
#
# LIBRARIES will be the libraries to be made
#
#############################################################################

LIBRARIES       =  libdpj.a librah_C.a libtmr.a libtmg.a libtmt.a

##############################################################################
#
# Formulas for making the libraries.  The object files are removed after
#    the archives are made.
#
#     For any library libXYZ.a the Fortran files are in directory libXYZ,
#     and there should be no .o files in the subdirectory.
#
##############################################################################

MAKEFLIBS = \
	cd  $(@:.a= ) ; \
	$(FC) $(OPT) $(FFLAGS) -c *.f* ; \
	rm  -f ../$@  ; \
	ar rv $(INSTALLDIR)/$@   *.o ; \
	rm -f *.o ; \
        echo Library $@ is now ready 

MAKECLIBS = \
        cd  $(@:.a= ) ; \
        $(CC) $(OPT) $(CFLAGS) -c *.c ; \
        rm  -f ../$@  ; \
        ar rv $(INSTALLDIR)/$@   *.o ; \
        rm -f *.o ; \
        echo Librariy $@ is now ready

##############################################################################
#
# Targets to build
#
##############################################################################

all : $(LIBRARIES)

libdpj.a : libdpj/*.f
	$(MAKEFLIBS)
libtmr.a : libtmr/*.f
	$(MAKEFLIBS)
libtmg.a : libtmg/*.f
	$(MAKEFLIBS)
libtmt.a : libtmt/*.f
	$(MAKEFLIBS)
librah_C.a : librah_C/*.c
	$(MAKECLIBS)
