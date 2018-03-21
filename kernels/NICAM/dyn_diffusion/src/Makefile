################################################################################
#
# Makefile grand central
#
################################################################################

TOPDIR   = ../../../..
BUILDDIR = .

VPATH = $(BUILDDIR)

include $(TOPDIR)/sysdep/Makedef.$(IAB_SYS)
include $(TOPDIR)/kernels/NICAM/Mkinclude

all:
	$(MAKE) makedir
	$(MAKE) makebin

allclean: clean cleandir

PRJ = dyn_diffusion.exe

PRG = main

MODS = \
	dumpio.o \
	mod_precision.o \
	mod_misc.o      \
	mod_oprt.o

makedir:
	mkdir -p $(BUILDDIR)
	mkdir -p $(BINDIR)

makebin: makedir $(PRJ)



dumpio.o        : dumpio.c
mod_precision.o : mod_precision.f90
mod_misc.o      : mod_misc.f90 problem_size.inc mod_precision.o
mod_oprt.o      : mod_oprt.f90 mod_misc.o mod_precision.o

main.o : main.f90 mod_precision.o dumpio.o mod_misc.o mod_oprt.o

$(BUILDDIR)/$(PRJ) : $(PRG).o $(BUILDDIR)/$(MODS)
	$(LD) $(LDFLAGS) -o $@ $^
	$(INSTALL) $@ $(BINDIR)
	@echo; echo "##### PROGRAM $(notdir $@) is installed Successfully #####"; echo



cleandir:
	-rm $(BINDIR)/*

clean:
	rm -f $(PRJ)
	rm -f *~ *.mod *.o *.lst *.L
	cd $(BUILDDIR); rm -f $(PRJ); rm -f *~ *.mod *.o *.lst *.L

.SUFFIXES: .o .f90 .c .mod
%.o : %.f90
	$(FC) $(FFLAGS) -o $(BUILDDIR)/$@ -c $<

%.o : %.c
	$(CC) $(CFLAGS) -o $(BUILDDIR)/$@ -c $<

%.o: %.mod

%.mod: %.f90
	$(make) $(patsubst %.f90,%.o,$<)

.PHONY : clean cleandir allclean
