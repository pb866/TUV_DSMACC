# Makefile for TUV 5.2.1
# Use with disort (discrete ordinate), or ps2str (2 stream approximation,
# pseudo-spherical correction)
#----------
# EXC      : name of executable
# INCLUDES : required include files
# USE_INCL : object files referencing include file (params)
# FOBJS    : all required object files that do not use the include file
#

intel := $(shell command -v ifort 2> /dev/null)
EXC = tuv
INCLUDES = params

USE_INCL = SRC/TUV.o \
		   SRC/grids.o \
		   SRC/rdinp.o SRC/rdetfl.o SRC/rdxs.o \
		   SRC/swphys.o SRC/swbiol.o SRC/RXN/swchem.o SRC/RXN/mcmext.o\
		   SRC/RXN/rxn.o SRC/RXN/rxn_ald.o SRC/RXN/rxn_ket.o SRC/RXN/rxn_dicar.o\
		   SRC/RXN/rxn_nit.o SRC/RXN/rxn_dinit.o SRC/RXN/rxn_rooh.o\
		   SRC/RXN/rxn_mult.o SRC/RXN/rxn_sci.o SRC/RXN/rxn_halo.o SRC/qys.o \
		   SRC/wshift.o \
		   SRC/vpair.o SRC/vptmp.o SRC/vpo3.o \
		   SRC/odrl.o SRC/odo3.o \
		   SRC/setaer.o SRC/setalb.o SRC/setcld.o SRC/setsnw.o \
		   SRC/setno2.o SRC/seto2.o SRC/setso2.o \
		   SRC/sphers.o  \
		   SRC/la_srb.o \
		   SRC/rtrans.o \
		   SRC/savout.o \
			 SRC/options.o

FOBJS = SRC/numer.o SRC/functs.o SRC/orbit.o

#----------
# FC   : FORTRAN compiler
#        Linux users:  try FC = g77
#        Cray users :  try FC = f90
# FC = pgf90
# FC = gfortran

# FFLAGS : command line options to compiler call (if not set, default is
#          probably some basic optimization level)
# FFLAGS = -fcheck=all #-Wall
FFLAGS = -cpp -fpp -fp-model strict -O3 -no-prec-div -static -xHost

# LIBS  : libraries required
# LIBS =

#----------

$(EXC):	$(FOBJS) $(USE_INCL)
	$(FC) $(FFLAGS) $(FOBJS) $(USE_INCL) $(LIBS) -o $@

$(USE_INCL):	$(INCLUDES)

clean:
		rm -f $(EXC) $(USE_INCL) $(FOBJS)

tidy: clean
	rm -f core *~ fort.*
