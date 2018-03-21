################################################################################
#
# ------ For K computer -----
#
################################################################################

##### Fortran setting
#FFLAGS is set in Mkinclude. FFLAGS_DEBUG is used if IAB_DEBUG=T

FC  = mpifrtpx
MPIFC = $(FC)

FFLAGS_FAST  = -Kfast,parallel,ocl,preex,array_private,noalias=s,mfunc=2 \
               -Qi -Qt -X03 -Ncompdisp -Koptmsg=1 -Cpp                   \
               -Kdynamic_iteration -Ksimd                                \
               -Kauto,threadsafe

#               -Qm -Nquickdbg -NRtrap

FFLAGS_DEBUG = -O0                                                   \
               -Qi -Qt -X03 -Ncompdisp -Koptmsg=2 -Cpp               \
               -Ec -Eg -Ha -He -Hf -Ho -Hs -Hu -Hx -v03s -v03d -v03o \
               -Ncheck_cache_arraysize -Ncheck_global



##### C setting
#CFLAGS is set in Mkinclude. CFLAGS_DEBUG is used if IAB_DEBUG=T

CC = mpifccpx
MPICC = $(CC)

CFLAGS_FAST  = -Kfast,parallel,ocl,preex,array_private,region_extension,restp=all -Ksimd
CFLAGS_DEBUG = -O0



##### Special setting

FFLAGS_OPENMP = -Kopenmp
# enable + disable parallel
# FFLAGS_OPENMP = -Kopenmp,noparallel

### Performance monitor
# disable
#PERF_MONIT = -Ntl_notrt -UFIPP -UFAPP -UFINEPA
# fipp
PERF_MONIT = -Ntl_trt -DFIPP -UFAPP -UFINEPA
# fapp
#PERF_MONIT = -Ntl_trt -UFIPP -DFAPP -UFINEPA
# fine PA
#PERF_MONIT = -Ntl_trt -UFIPP -UFAPP -DFINEPA

FFLAGS_SYSDEPEND = $(PERF_MONIT)
CFLAGS_SYSDEPEND = $(PERF_MONIT)

##### Linker setting

LD      = $(FC)
MPILD   = $(MPIFC)
LDFLAGS = $(FFLAGS)



################################################################################

##### for frontend
MAKE    = make
INSTALL = cp -u
AR      = ar
ARFLAGS = r
RANLIB  = ranlib
JOBSUB  = pjsub
