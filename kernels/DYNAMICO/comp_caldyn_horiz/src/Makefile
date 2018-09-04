################################################################################
#
# Makefile grand central
#
################################################################################

TOPDIR   = ../../../..
BUILDDIR = .

VPATH = $(BUILDDIR)

include $(TOPDIR)/sysdep/Makedef.$(IAB_SYS)
include $(TOPDIR)/kernels/DYNAMICO/Mkinclude

all:
	$(MAKE) makedir
	$(MAKE) makebin

allclean: clean cleandir

PRJ = comp_caldyn_horiz.exe

PRG = main

MODS = \
	dumpio.o \
	prec.o \
	mod_precision.o \
	mod_misc.o \
	caldyn_gcm.o 

makedir:
	mkdir -p $(BUILDDIR)
	mkdir -p $(BINDIR)

makebin: makedir $(PRJ)



mod_misc.o : problem_size.inc mod_precision.o
caldyn_gcm.o : prec.o mod_misc.o 
main.o     : prec.o mod_misc.o caldyn_gcm.o

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
