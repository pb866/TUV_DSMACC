# Makefile for TUV 5.2.1 as subroutine of DSMACC
# Use with disort (discrete ordinate), or ps2str (2 stream approximation,
# pseudo-spherical correction)
#----------
# EXC      : name of executable
# INCLUDES : required include files
# USE_INCL : object files referencing include file (params)
# FOBJS    : all required object files that do not use the include file
#

 #-fp-stack-check -check bou    nds -check arg_temp_created -check all #-warn all # -openmp



EXC = tuv

INCLUDES = params

USE_INCL = TUV521.o \
		   grids.o \
		   rdinp.o rdetfl.o rdxs.o \
		   swphys.o swbiol.o swchem.o mcmext.o gcext.o\
		   rxn_mcm.o rxn_ald.o rxn_ket.o rxn_dicar.o\
		   rxn_nit.o rxn_dinit.o rxn_rooh.o\
		   rxn_mult.o rxn_rad.o rxn_gc11.o qys.o \
		   wshift.o \
		   vpair.o vptmp.o vpo3.o \
		   odrl.o odo3.o \
		   setaer.o setalb.o setcld.o setsnw.o \
		   setno2.o seto2.o setso2.o \
		   sphers.o  \
		   la_srb.o \
		   rtrans.o \
		   savout.o
# Include in list above, if needed: rxn_test.o

FOBJS = numer.o functs.o orbit.o

#----------
# FC   : FORTRAN compiler
#        Linux users:  try FC = g77
#        Cray users :  try FC = f90
# FC = pgf90
# FC = gfortran
#FC = ifort

# FFLAGS : command line options to compiler call (if not set, default is
#          probably some basic optimization level)
# FFLAGS = -fcheck=all #-Wall -cpp -mcmodel medium -fpp
# FFLAGS = -cpp -fpp -fp-model strict -O3 -no-prec-div -static -xHost

# LIBS  : libraries required
# LIBS =

#----------

$(EXC):	compiler	$(FOBJS) $(USE_INCL)
		$(FC) $(FFLAGS) $(FOBJS) $(USE_INCL) $(LIBS) -o $@

$(USE_INCL):	$(INCLUDES)

.f.o: compiler
		$(FC) $(FFLAGS) -c $*.f

clean:
		rm -f core $(EXC) $(USE_INCL) $(FOBJS)

tidy: clean
		rm -f *~ fort.*


compiler:
ifndef intel
        @echo 'using gfortran'
        $(eval export FC=gfortran)
        $(eval export F90=gfortran)
        $(eval export FFLAGS=-fcheck=all)
else
        @echo 'using ifort'
        $(eval export FC=ifort)
        $(eval export F90=ifort)
        $(eval export FFLAGS = -cpp -fpp -fp-model strict -O3 -no-prec-div -static -xHost)
endif

