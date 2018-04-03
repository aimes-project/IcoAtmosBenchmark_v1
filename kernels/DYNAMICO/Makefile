################################################################################
#
# Build all DYNAMICO kernels
#
################################################################################
TOPDIR   = ../..

include $(TOPDIR)/sysdep/Makedef.$(IAB_SYS)

KERNELS = \
	comp_caldyn_horiz \
	comp_caldyn_vert \
	comp_geopot \
	comp_pvort



runshell = run-$(IAB_SYS).sh

ifeq ($(IAB_SYS),Linux64-gnu-ompi)
   runshell = run-linux.sh
endif

ifeq ($(IAB_SYS),Linux64-intel-impi)
   runshell = run-linux.sh
endif

ifeq ($(IAB_SYS),Linux64-intelmic-impi)
   runshell = run-linux.sh
endif

ifeq ($(IAB_SYS),Linux64-pgi-ompi)
   runshell = run-linux.sh
endif

ifeq ($(IAB_SYS),Ubuntu-gnu-ompi)
   runshell = run-linux.sh
endif

ifeq ($(IAB_SYS),Ubuntu-intel-ompi)
   runshell = run-linux.sh
endif

ifeq ($(IAB_SYS),Ubuntu-pgi-ompi)
   runshell = run-linux.sh
endif

ifeq ($(IAB_SYS),MacOSX-gnu-ompi)
   runshell = run-linux.sh
endif



all:
	for d in $(KERNELS); do  $(MAKE) -C $$d/src $@; done

clean:
	for d in $(KERNELS); do  $(MAKE) -C $$d/src $@; done

data:
	for d in $(KERNELS); do  cd $$d/data; bash download.sh ;cd -; done

run:
	for d in $(KERNELS); do  cd $$d/run; $(JOBSUB) $(runshell) ;cd -; done
