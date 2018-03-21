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

PRJ = communication.exe

PRG = main

MODS = \
	mod_precision.o \
	mod_stdio.o     \
	mod_prof.o      \
	mod_process.o   \
	mod_adm.o       \
	mod_comm.o

makedir:
	mkdir -p $(BUILDDIR)
	mkdir -p $(BINDIR)

makebin: makedir $(PRJ)



mod_precision.o : mod_precision.f90
mod_stdio.o     : mod_stdio.f90
mod_process.o   : mod_process.f90 mod_precision.o mod_stdio.o
mod_prof.o      : mod_prof.f90    mod_precision.o mod_stdio.o mod_process.o
mod_adm.o       : mod_adm.f90     mod_precision.o mod_stdio.o mod_process.o mod_prof.o
mod_comm.o      : mod_comm.f90    mod_precision.o mod_stdio.o mod_process.o mod_prof.o mod_adm.o

main.o : main.f90 mod_precision.o mod_stdio.o mod_process.o mod_prof.o mod_adm.o mod_comm.o

$(BUILDDIR)/$(PRJ) : $(PRG).o $(BUILDDIR)/$(MODS)
	$(MPILD) $(LDFLAGS) -o $@ $^
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
	$(MPIFC) $(FFLAGS) -o $(BUILDDIR)/$@ -c $<

%.o : %.c
	$(MPICC) $(CFLAGS) -o $(BUILDDIR)/$@ -c $<

%.o: %.mod

%.mod: %.f90
	$(make) $(patsubst %.f90,%.o,$<)

.PHONY : clean cleandir allclean
