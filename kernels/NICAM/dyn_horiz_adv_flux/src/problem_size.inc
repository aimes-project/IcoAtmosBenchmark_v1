!-------------------------------------------------------------------------------
!
! Problem size and global parameters
!
!-------------------------------------------------------------------------------

  integer,  public, parameter :: ADM_NSYS     = 32
  integer,  public, parameter :: ADM_MAXFNAME = 1024
  integer,  public            :: ADM_LOG_FID  = 6

  integer,  public, parameter :: IO_FREAD = 0

  !--- Identifier of triangle element (i-axis-side or j-axis side)
  integer,  public, parameter :: TI  = 1
  integer,  public, parameter :: TJ  = 2

  !--- Identifier of line element (i-axis-side, ij-axis side, or j-axis side)
  integer,  public, parameter :: AI  = 1
  integer,  public, parameter :: AIJ = 2
  integer,  public, parameter :: AJ  = 3

  !--- Identifier of 1 variable
  integer,  public, parameter :: K0  = 1

  integer,  public, parameter :: ADM_nxyz = 3 ! dimension of the spacial vector

  !--- region
  integer,  public            :: ADM_lall      = 1     ! number of regular region per process
  integer,  public, parameter :: ADM_lall_pl   = 2     ! number of pole    region per process

  !--- horizontal grid
  integer,  public            :: ADM_gall      = 16900 ! number of horizontal grid per regular region
  integer,  public            :: ADM_gall_1d   = 130   ! number of horizontal grid (1D)
  integer,  public            :: ADM_gmin      = 2     ! start index of 1D horizontal grid
  integer,  public            :: ADM_gmax      = 129   ! end   index of 1D horizontal grid

  integer,  public            :: ADM_gall_pl   = 6     ! number of horizontal grid for pole region
  integer,  public, parameter :: ADM_gslf_pl   = 1     ! index for pole point
  integer,  public, parameter :: ADM_gmin_pl   = 2     ! start index of grid around the pole point
  integer,  public            :: ADM_gmax_pl   = 6     ! end   index of grid around the pole point
  integer,  public, parameter :: ADM_vlink     = 5     ! number of grid around the pole point

  !--- vertical grid
  integer,  public            :: ADM_vlayer    = 40    ! number of vertical layer
  integer,  public            :: ADM_kall      = 42    ! number of vertical grid
  integer,  public            :: ADM_kmin      = 2     ! start index of vertical grid
  integer,  public            :: ADM_kmax      = 41    ! end   index of vertical grid

  logical,  public            :: ADM_have_pl     = .true. ! this ID manages pole region?
  logical,  public            :: ADM_have_sgp(1) = .true. ! region have singlar point?

  !--- constant parameters
  real(RP), public            :: CONST_PI    = 3.14159265358979_RP !< pi
  real(RP), public            :: CONST_EPS   = 1.E-16_RP           !< small number
  real(RP), public            :: CONST_HUGE  = 1.E+30_RP           !< huge  number
  real(RP), public            :: CONST_GRAV  = 9.80616_RP          !< Gravitational accerlaration of the Earth [m/s2]
  real(RP), public            :: CONST_Rdry  =   287.0_RP          !< Gas constant of air
  real(RP), public            :: CONST_CVdry =   717.5_RP          !< Specific heat of air (consant volume)
  integer,  public            :: NON_HYDRO_ALPHA = 1               !< Nonhydrostatic/hydrostatic flag

  !--- mod_grd
  integer,  public, parameter :: XDIR = 1
  integer,  public, parameter :: YDIR = 2
  integer,  public, parameter :: ZDIR = 3

  integer,  public, parameter :: GRD_LAT = 1
  integer,  public, parameter :: GRD_LON = 2

  real(RP), public            :: GRD_rscale = 6.37122E+6_RP !< radius of the planet [m]

  integer,  public, parameter :: GRD_grid_type_on_sphere = 1
  integer,  public, parameter :: GRD_grid_type_on_plane  = 2
  integer,  public            :: GRD_grid_type           = GRD_grid_type_on_sphere

  character(len=ADM_MAXFNAME), public :: vgrid_fname   ='./vgrid40_600m_24km.dat'

  !--- mod_gmtr
  character(len=ADM_NSYS),     public :: GMTR_polygon_type = 'ON_SPHERE'

  integer,  public, parameter :: GMTR_p_nmax = 8

  integer,  public, parameter :: GMTR_p_AREA  = 1
  integer,  public, parameter :: GMTR_p_RAREA = 2
  integer,  public, parameter :: GMTR_p_IX    = 3
  integer,  public, parameter :: GMTR_p_IY    = 4
  integer,  public, parameter :: GMTR_p_IZ    = 5
  integer,  public, parameter :: GMTR_p_JX    = 6
  integer,  public, parameter :: GMTR_p_JY    = 7
  integer,  public, parameter :: GMTR_p_JZ    = 8

  integer,  public, parameter :: GMTR_t_nmax = 5

  integer,  public, parameter :: GMTR_t_AREA  = 1
  integer,  public, parameter :: GMTR_t_RAREA = 2
  integer,  public, parameter :: GMTR_t_W1    = 3
  integer,  public, parameter :: GMTR_t_W2    = 4
  integer,  public, parameter :: GMTR_t_W3    = 5

  integer,  public, parameter :: GMTR_a_nmax    = 12
  integer,  public, parameter :: GMTR_a_nmax_pl = 18

  integer,  public, parameter :: GMTR_a_HNX  = 1
  integer,  public, parameter :: GMTR_a_HNY  = 2
  integer,  public, parameter :: GMTR_a_HNZ  = 3
  integer,  public, parameter :: GMTR_a_HTX  = 4
  integer,  public, parameter :: GMTR_a_HTY  = 5
  integer,  public, parameter :: GMTR_a_HTZ  = 6
  integer,  public, parameter :: GMTR_a_TNX  = 7
  integer,  public, parameter :: GMTR_a_TNY  = 8
  integer,  public, parameter :: GMTR_a_TNZ  = 9
  integer,  public, parameter :: GMTR_a_TTX  = 10
  integer,  public, parameter :: GMTR_a_TTY  = 11
  integer,  public, parameter :: GMTR_a_TTZ  = 12

  integer,  public, parameter :: GMTR_a_TN2X = 13
  integer,  public, parameter :: GMTR_a_TN2Y = 14
  integer,  public, parameter :: GMTR_a_TN2Z = 15
  integer,  public, parameter :: GMTR_a_TT2X = 16
  integer,  public, parameter :: GMTR_a_TT2Y = 17
  integer,  public, parameter :: GMTR_a_TT2Z = 18

  integer,  public, parameter :: P_RAREA = GMTR_p_RAREA
  integer,  public, parameter :: T_RAREA = GMTR_t_RAREA
  integer,  public, parameter :: W1      = GMTR_t_W1
  integer,  public, parameter :: W2      = GMTR_t_W2
  integer,  public, parameter :: W3      = GMTR_t_W3
  integer,  public, parameter :: HNX     = GMTR_a_HNX
  integer,  public, parameter :: HNY     = GMTR_a_HNY
  integer,  public, parameter :: HNZ     = GMTR_a_HNZ
  integer,  public, parameter :: HTX     = GMTR_a_HTX
  integer,  public, parameter :: TNX     = GMTR_a_TNX
  integer,  public, parameter :: TN2X    = GMTR_a_TN2X

  integer,  public :: SET_iteration = 1
  integer,  public :: SET_prc_me    = 1
  real(DP), public :: SET_dt_large  = 360.0_DP
  real(DP), public :: SET_dt_small  =  60.0_DP
