!-------------------------------------------------------------------------------
!
! Problem size and global parameters
!
!-------------------------------------------------------------------------------

!!$   integer, public :: iim = 23
!!$   integer, public :: jjm = 25
!!$   integer, public :: llm = 1
!!$   integer, public :: ij_begin = 48
!!$   integer, public :: ij_end = 528
!!$   integer, public :: ij_begin_ext = 24
!!$   integer, public :: ij_end_ext = 552
!!$   integer, public :: ll_begin = 1
!!$   integer, public :: ll_end = 1

  !! below are read from input data file.
  integer, public :: iim = -huge(1)
  integer, public :: jjm = -huge(1)
  integer, public :: llm = -huge(1)

  integer, public :: ij_begin     = -huge(1)
  integer, public :: ij_end       =  huge(1)
  integer, public :: ij_begin_ext = -huge(1)
  integer, public :: ij_end_ext   =  huge(1)
  integer, public :: ll_begin     = -huge(1)
  integer, public :: ll_end       =  huge(1)

  integer, public :: ll_beginp1   = -huge(1)
  integer, public :: ll_endm1     =  huge(1)

  integer, public :: ij_omp_begin = -huge(1)
  integer, public :: ij_omp_end   =  huge(1)

  INTEGER, public :: t_right = -huge(1)
  INTEGER, public :: t_rup   = -huge(1)
  INTEGER, public :: t_lup   = -huge(1)
  INTEGER, public :: t_left  = +huge(1)
  INTEGER, public :: t_ldown = +huge(1)
  INTEGER, public :: t_rdown = +huge(1)

  INTEGER, public :: u_right = -huge(1)
  INTEGER, public :: u_rup   = -huge(1)
  INTEGER, public :: u_lup   = -huge(1)
  INTEGER, public :: u_left  = +huge(1)
  INTEGER, public :: u_ldown = +huge(1)
  INTEGER, public :: u_rdown = +huge(1)

  integer, public :: z_rup   = -huge(1)
  integer, public :: z_up    = -huge(1)
  integer, public :: z_lup   = -huge(1)
  integer, public :: z_ldown = +huge(1)
  integer, public :: z_down  = +huge(1)
  integer, public :: z_rdown = +huge(1)

!!$  REAL(rstd),public :: g=9.80616
  REAL(rstd),public :: g

  INTEGER, public, PARAMETER :: energy=1, enstrophy=2
!!$  INTEGER, public :: caldyn_conserv = energy
  INTEGER, public :: caldyn_conserv

!!$  logical, public :: boussinesq = .true.
  logical, public :: boussinesq

  logical, public :: is_omp_first_level = .true.

  !! parameters from geometry.f90
  integer, public, parameter :: ne_right=1
  integer, public, parameter :: ne_rup=-1
  integer, public, parameter :: ne_lup=1
  integer, public, parameter :: ne_left=-1
  integer, public, parameter :: ne_ldown=1
  integer, public, parameter :: ne_rdown=-1

  !! parameters from metric.f90
  integer, public, parameter :: vrup=1
  integer, public, parameter :: vup=2
  integer, public, parameter :: vlup=3
  integer, public, parameter :: vldown=4
  integer, public, parameter :: vdown=5
  integer, public, parameter :: vrdown=6

  integer, public, parameter :: IO_FREAD = 0



  integer, public :: SET_prc_me    = 1
  integer, public :: SET_iteration = 1000



  integer,  public, parameter :: ADM_NSYS     = 32
  integer,  public            :: ADM_LOG_FID  = 6

  !--- constant parameters
  real(RP), public            :: CONST_PI    = 3.14159265358979_RP !< pi
  real(RP), public            :: CONST_EPS   = 1.E-16_RP           !< small number
  real(RP), public            :: CONST_HUGE  = 1.E+30_RP           !< huge  number



!!$ Local Variables:
!!$ mode: f90
!!$ f90-do-indent: 2
!!$ f90-if-indent: 2
!!$ f90-program-indent: 2
!!$ f90-type-indent: 2
!!$ f90-directive-comment-re: "!\$omp "
!!$ end:


