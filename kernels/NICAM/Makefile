################################################################################
#
# Build all NICAM kernels
#
################################################################################
TOPDIR   = ../..

include $(TOPDIR)/sysdep/Makedef.$(IAB_SYS)

KERNELS = \
	dyn_diffusion \
	dyn_divdamp \
	dyn_vi_rhow_solver \
	dyn_vert_adv_limiter \
	dyn_horiz_adv_flux \
	dyn_horiz_adv_limiter \
	dyn_metrics \
	communication



runshell = run-$(IAB_SYS).sh

ifeq ($(IAB_SYS),Linux64-gnu-ompi)
   runshell = run-linux-torque.sh
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
