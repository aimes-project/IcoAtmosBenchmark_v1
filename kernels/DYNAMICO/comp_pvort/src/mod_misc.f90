module mod_misc
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use prec
  use mod_precision
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public, including parameters
  !
  include 'problem_size.inc'

  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: DEBUG_rapstart
  public :: DEBUG_rapend
  public :: DEBUG_rapreport
  public :: DEBUG_valuecheck

  interface DEBUG_valuecheck
     module procedure DEBUG_valuecheck_1D
     module procedure DEBUG_valuecheck_2D
     module procedure DEBUG_valuecheck_3D
     module procedure DEBUG_valuecheck_4D
     module procedure DEBUG_valuecheck_5D
     module procedure DEBUG_valuecheck_6D
  end interface DEBUG_valuecheck

  public :: MISC_make_idstr        !--- make file name with a number
  public :: MISC_get_available_fid !--- get an available file ID

  public :: ADM_proc_stop

  public :: CONST_setup

  public :: trace_start, trace_end
  public :: wait_message, test_message


  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer,             public :: EX_CSTEP_diffusion = 0
  integer,             public :: EX_TSTEP_diffusion = 60
  integer,             public :: EX_CSTEP_divdamp3d = 0
  integer,             public :: EX_TSTEP_divdamp3d = 140
  integer,             public :: EX_fid
  integer,             public :: EX_err
  character(len=1024), public :: EX_fname
  character(len=16),   public :: EX_item
  real(RP),            public :: EX_max
  real(RP),            public :: EX_min
  real(RP),            public :: EX_sum

  type t_message
    !! dummy here
  end type t_message
  public t_message

  REAL(rstd), public, allocatable :: Av(:) ! area of dual mesk cell
  real(rstd), public, allocatable :: de(:) ! distance from a neighbour == lenght of an edge of the dual mesh
  real(rstd), public, allocatable :: Riv2(:,:) ! weigt
  real(rstd), public, allocatable :: fv(:) ! coriolis (evaluted on a vertex)

  REAL(rstd), public, allocatable :: mass_dak(:) ! llm
  REAL(rstd), public, allocatable :: mass_dbk(:) ! llm
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: DEBUG_rapid

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,                 private, parameter :: DEBUG_rapnlimit = 100
  integer,                 private            :: DEBUG_rapnmax   = 0
  character(len=ADM_NSYS), private            :: DEBUG_rapname(DEBUG_rapnlimit)
  real(8),                 private            :: DEBUG_raptstr(DEBUG_rapnlimit)
  real(8),                 private            :: DEBUG_rapttot(DEBUG_rapnlimit)
  integer,                 private            :: DEBUG_rapnstr(DEBUG_rapnlimit)
  integer,                 private            :: DEBUG_rapnend(DEBUG_rapnlimit)

!!$#ifdef _FIXEDINDEX_
!!$  real(DP), public              :: GRD_gz   (ADM_kall)
!!$  real(DP), public              :: GRD_gzh  (ADM_kall)
!!$  real(DP), public              :: GRD_dgz  (ADM_kall)
!!$  real(DP), public              :: GRD_dgzh (ADM_kall)
!!$  real(DP), public              :: GRD_rdgz (ADM_kall)
!!$  real(DP), public              :: GRD_rdgzh(ADM_kall)
!!$
!!$  real(DP), public              :: GRD_afact(ADM_kall)
!!$  real(DP), public              :: GRD_bfact(ADM_kall)
!!$  real(DP), public              :: GRD_cfact(ADM_kall)
!!$  real(DP), public              :: GRD_dfact(ADM_kall)
!!$
!!$  real(RP), public              :: GRD_xr   (ADM_gall   ,ADM_kall,ADM_lall   ,ADM_AI:ADM_AJ,XDIR:ZDIR)
!!$  real(RP), public              :: GRD_xr_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,              XDIR:ZDIR)
!!$
!!$  real(RP), public              :: GMTR_p   (ADM_gall   ,ADM_KNONE,ADM_lall   ,              GMTR_p_nmax   )
!!$  real(RP), public              :: GMTR_p_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,              GMTR_p_nmax   )
!!$  real(RP), public              :: GMTR_t   (ADM_gall   ,ADM_KNONE,ADM_lall   ,ADM_TI:ADM_TJ,GMTR_t_nmax   )
!!$  real(RP), public              :: GMTR_t_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,              GMTR_t_nmax   )
!!$  real(RP), public              :: GMTR_a   (ADM_gall   ,ADM_KNONE,ADM_lall   ,ADM_AI:ADM_AJ,GMTR_a_nmax   )
!!$  real(RP), public              :: GMTR_a_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,              GMTR_a_nmax_pl)
!!$
!!$  real(RP), public              :: VMTR_RGSQRTH     (ADM_gall   ,ADM_kall,ADM_lall   )
!!$  real(RP), public              :: VMTR_RGSQRTH_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
!!$  real(RP), public              :: VMTR_RGAM        (ADM_gall   ,ADM_kall,ADM_lall   )
!!$  real(RP), public              :: VMTR_RGAM_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
!!$  real(RP), public              :: VMTR_RGAMH       (ADM_gall   ,ADM_kall,ADM_lall   )
!!$  real(RP), public              :: VMTR_RGAMH_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
!!$  real(RP), public              :: VMTR_RGSGAM2     (ADM_gall   ,ADM_kall,ADM_lall   )
!!$  real(RP), public              :: VMTR_RGSGAM2_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
!!$  real(RP), public              :: VMTR_RGSGAM2H    (ADM_gall   ,ADM_kall,ADM_lall   )
!!$  real(RP), public              :: VMTR_RGSGAM2H_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
!!$  real(RP), public              :: VMTR_GSGAM2H     (ADM_gall   ,ADM_kall,ADM_lall   )
!!$  real(RP), public              :: VMTR_GSGAM2H_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
!!$  real(RP), public              :: VMTR_C2WfactGz   (ADM_gall   ,ADM_kall,6,ADM_lall   )
!!$  real(RP), public              :: VMTR_C2WfactGz_pl(ADM_gall_pl,ADM_kall,6,ADM_lall_pl)
!!$#else
!!$  real(DP), public, allocatable :: GRD_gz   (:) ! gsi (z-star) coordinate
!!$  real(DP), public, allocatable :: GRD_gzh  (:) ! gsi (z-star) coordinate at the half point
!!$  real(DP), public, allocatable :: GRD_dgz  (:) ! d(gsi)
!!$  real(DP), public, allocatable :: GRD_dgzh (:) ! d(gsi) at the half point
!!$  real(DP), public, allocatable :: GRD_rdgz (:)
!!$  real(DP), public, allocatable :: GRD_rdgzh(:)
!!$
!!$  real(DP), public, allocatable :: GRD_afact(:) ! From the cell center value to the cell wall value
!!$  real(DP), public, allocatable :: GRD_bfact(:) !    A(k-1/2) = ( afac(k) A(k) + bfac(k) * A(k-1) ) / 2
!!$  real(DP), public, allocatable :: GRD_cfact(:) ! From the cell wall value to the cell center value
!!$  real(DP), public, allocatable :: GRD_dfact(:) !    A(k) = ( cfac(k) A(k+1/2) + dfac(k) * A(k-1/2) ) / 2
!!$
!!$  real(RP), public, allocatable :: GRD_xr   (:,:,:,:,:) ! mass centroid position
!!$  real(RP), public, allocatable :: GRD_xr_pl(:,:,:,:)
!!$
!!$  real(RP), public, allocatable :: GMTR_p   (:,:,:,:)   ! geometrics for the cell point
!!$  real(RP), public, allocatable :: GMTR_p_pl(:,:,:,:)
!!$  real(RP), public, allocatable :: GMTR_t   (:,:,:,:,:) ! geometrics for the cell vertex
!!$  real(RP), public, allocatable :: GMTR_t_pl(:,:,:,:)
!!$  real(RP), public, allocatable :: GMTR_a   (:,:,:,:,:) ! geometrics for the cell arc
!!$  real(RP), public, allocatable :: GMTR_a_pl(:,:,:,:)
!!$
!!$  real(RP), public, allocatable :: VMTR_RGSQRTH     (:,:,:)
!!$  real(RP), public, allocatable :: VMTR_RGSQRTH_pl  (:,:,:)
!!$  real(RP), public, allocatable :: VMTR_RGAM        (:,:,:)
!!$  real(RP), public, allocatable :: VMTR_RGAM_pl     (:,:,:)
!!$  real(RP), public, allocatable :: VMTR_RGAMH       (:,:,:)
!!$  real(RP), public, allocatable :: VMTR_RGAMH_pl    (:,:,:)
!!$  real(RP), public, allocatable :: VMTR_RGSGAM2     (:,:,:)
!!$  real(RP), public, allocatable :: VMTR_RGSGAM2_pl  (:,:,:)
!!$  real(RP), public, allocatable :: VMTR_RGSGAM2H    (:,:,:)
!!$  real(RP), public, allocatable :: VMTR_RGSGAM2H_pl (:,:,:)
!!$  real(RP), public, allocatable :: VMTR_GSGAM2H     (:,:,:)
!!$  real(RP), public, allocatable :: VMTR_GSGAM2H_pl  (:,:,:)
!!$  real(RP), public, allocatable :: VMTR_C2WfactGz   (:,:,:,:)
!!$  real(RP), public, allocatable :: VMTR_C2WfactGz_pl(:,:,:,:)
!!$#endif

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  function DEBUG_rapid( rapname ) result(id)
    implicit none

    character(len=*), intent(in) :: rapname

    integer :: id
    !---------------------------------------------------------------------------

    if ( DEBUG_rapnmax >= 1 ) then
       do id = 1, DEBUG_rapnmax
          if( trim(rapname) == trim(DEBUG_rapname(id)) ) return
       enddo
    endif

    DEBUG_rapnmax     = DEBUG_rapnmax + 1
    id                = DEBUG_rapnmax
    DEBUG_rapname(id) = trim(rapname)
    DEBUG_raptstr(id) = 0.D0
    DEBUG_rapttot(id) = 0.D0
    DEBUG_rapnstr(id) = 0
    DEBUG_rapnend(id) = 0

  end function DEBUG_rapid

  !-----------------------------------------------------------------------------
  subroutine DEBUG_rapstart( rapname )
    implicit none

    character(len=*), intent(in) :: rapname

    real(8) :: time
    integer :: id
    !---------------------------------------------------------------------------

    id = DEBUG_rapid( rapname )

    call CPU_TIME(time)

    DEBUG_raptstr(id) = time
    DEBUG_rapnstr(id) = DEBUG_rapnstr(id) + 1

    !write(ADM_LOG_FID,*) rapname, DEBUG_rapnstr(id)

#ifdef _FAPP_
    call fapp_start( rapname, id, 1 )
#endif

    return
  end subroutine DEBUG_rapstart

  !-----------------------------------------------------------------------------
  subroutine DEBUG_rapend( rapname )
    implicit none

    character(len=*), intent(in) :: rapname

    real(8) :: time
    integer :: id
    !---------------------------------------------------------------------------

    id = DEBUG_rapid( rapname )

    call CPU_TIME(time)

    DEBUG_rapttot(id) = DEBUG_rapttot(id) + ( time-DEBUG_raptstr(id) )
    DEBUG_rapnend(id) = DEBUG_rapnend(id) + 1

#ifdef _FAPP_
    call fapp_stop( rapname, id, 1 )
#endif

    return
  end subroutine DEBUG_rapend

  !-----------------------------------------------------------------------------
  subroutine DEBUG_rapreport
    implicit none

    integer :: id
    !---------------------------------------------------------------------------

    if ( DEBUG_rapnmax >= 1 ) then

       do id = 1, DEBUG_rapnmax
          if ( DEBUG_rapnstr(id) /= DEBUG_rapnend(id) ) then
              write(*,*) '*** Mismatch Report',id,DEBUG_rapname(id),DEBUG_rapnstr(id),DEBUG_rapnend(id)
          endif
       enddo

       write(ADM_LOG_FID,*)
       write(ADM_LOG_FID,*) '*** Computational Time Report'

       do id = 1, DEBUG_rapnmax
          write(ADM_LOG_FID,'(1x,A,I3.3,A,A,A,F10.3,A,I7)') &
          '*** ID=',id,' : ',DEBUG_rapname(id),' T=',DEBUG_rapttot(id),' N=',DEBUG_rapnstr(id)
       enddo

    else
       write(ADM_LOG_FID,*)
       write(ADM_LOG_FID,*) '*** Computational Time Report: NO item.'
    endif

    return
  end subroutine DEBUG_rapreport

  !-----------------------------------------------------------------------------
  subroutine DEBUG_valuecheck_1D( &
       varname, &
       var      )
    implicit none

    character(len=*),  intent(in)  :: varname
    real(RP),          intent(in)  :: var(:)
    !---------------------------------------------------------------------------

     EX_item = trim  (varname)
     EX_max  = maxval(var)
     EX_min  = minval(var)
     EX_sum  = sum   (var)
     write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum

    return
  end subroutine DEBUG_valuecheck_1D

  !-----------------------------------------------------------------------------
  subroutine DEBUG_valuecheck_2D( &
       varname, &
       var      )
    implicit none

    character(len=*),  intent(in)  :: varname
    real(RP),          intent(in)  :: var(:,:)
    !---------------------------------------------------------------------------

     EX_item = trim  (varname)
     EX_max  = maxval(var)
     EX_min  = minval(var)
     EX_sum  = sum   (var)
     write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum

    return
  end subroutine DEBUG_valuecheck_2D

  !-----------------------------------------------------------------------------
  subroutine DEBUG_valuecheck_3D( &
       varname, &
       var      )
    implicit none

    character(len=*),  intent(in)  :: varname
    real(RP),          intent(in)  :: var(:,:,:)
    !---------------------------------------------------------------------------

     EX_item = trim  (varname)
     EX_max  = maxval(var)
     EX_min  = minval(var)
     EX_sum  = sum   (var)
     write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum

    return
  end subroutine DEBUG_valuecheck_3D

  !-----------------------------------------------------------------------------
  subroutine DEBUG_valuecheck_4D( &
       varname, &
       var      )
    implicit none

    character(len=*),  intent(in)  :: varname
    real(RP),          intent(in)  :: var(:,:,:,:)
    !---------------------------------------------------------------------------

     EX_item = trim  (varname)
     EX_max  = maxval(var)
     EX_min  = minval(var)
     EX_sum  = sum   (var)
     write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum

    return
  end subroutine DEBUG_valuecheck_4D

  !-----------------------------------------------------------------------------
  subroutine DEBUG_valuecheck_5D( &
       varname, &
       var      )
    implicit none

    character(len=*),  intent(in)  :: varname
    real(RP),          intent(in)  :: var(:,:,:,:,:)
    !---------------------------------------------------------------------------

     EX_item = trim  (varname)
     EX_max  = maxval(var)
     EX_min  = minval(var)
     EX_sum  = sum   (var)
     write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum

    return
  end subroutine DEBUG_valuecheck_5D

  !-----------------------------------------------------------------------------
  subroutine DEBUG_valuecheck_6D( &
       varname, &
       var      )
    implicit none

    character(len=*),  intent(in)  :: varname
    real(RP),          intent(in)  :: var(:,:,:,:,:,:)
    !---------------------------------------------------------------------------

     EX_item = trim  (varname)
     EX_max  = maxval(var)
     EX_min  = minval(var)
     EX_sum  = sum   (var)
     write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum

    return
  end subroutine DEBUG_valuecheck_6D

  !-----------------------------------------------------------------------------
  subroutine MISC_make_idstr( &
       str,    &
       prefix, &
       ext,    &
       numID,  &
       digit   )
    implicit none

    character(len=*),  intent(out) :: str    !< combined extention string
    character(len=*),  intent(in)  :: prefix !< prefix
    character(len=*),  intent(in)  :: ext    !< extention ( e.g. .rgn )
    integer,           intent(in)  :: numID  !< number
    integer, optional, intent(in)  :: digit  !< digit

    logical, parameter            :: NSTR_ZERO_START = .true. ! number of separated file starts from 0 ?
    integer, parameter            :: NSTR_MAX_DIGIT  = 6      ! digit of separated file

    character(len=128) :: rankstr
    integer            :: setdigit
    !---------------------------------------------------------------------------

    if ( NSTR_ZERO_START ) then
       write(rankstr,'(I128.128)') numID-1
    else
       write(rankstr,'(I128.128)') numID
    endif

    if ( present(digit) ) then
       setdigit = digit
    else
       setdigit = NSTR_MAX_DIGIT
    endif

    rankstr(1:setdigit) = rankstr(128-(setdigit-1):128)
    rankstr(setdigit+1:128) = ' '

    str = trim(prefix)//'.'//trim(ext)//trim(rankstr) ! -> prefix.ext00000

    return
  end subroutine MISC_make_idstr

  !-----------------------------------------------------------------------------
  !> Search and get available machine id
  !> @return fid
  function MISC_get_available_fid() result(fid)
    implicit none

    integer :: fid

    integer, parameter :: min_fid =  7 !< minimum available fid
    integer, parameter :: max_fid = 99 !< maximum available fid

    logical :: i_opened
    !---------------------------------------------------------------------------

    do fid = min_fid, max_fid
       inquire(fid,opened=i_opened)
       if( .NOT. i_opened ) return
    enddo

  end function MISC_get_available_fid

  !-----------------------------------------------------------------------------
  subroutine ADM_proc_stop

    stop

  end subroutine ADM_proc_stop

  !-----------------------------------------------------------------------------
  subroutine CONST_setup
    implicit none

    write(ADM_LOG_FID,*) '*** setup constants'

    CONST_PI   = 4.E0_RP * atan( 1.0_RP )
    CONST_EPS  = epsilon(0.0_RP)
    CONST_HUGE =    huge(0.0_RP)

    write(ADM_LOG_FID,*) '*** PI   = ', CONST_PI
    write(ADM_LOG_FID,*) '*** EPS  = ', CONST_EPS
    write(ADM_LOG_FID,*) '*** HUGE = ', CONST_HUGE

    return
  end subroutine CONST_setup


  !! from trace.F90
  SUBROUTINE trace_start(name)
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: name 
#ifdef VTRACE
#include <vt_user.inc>
#endif 

!$OMP MASTER
#ifdef VTRACE
     VT_USER_START(name)
#endif
!$OMP END MASTER

  END SUBROUTINE trace_start    


  !! from trace.F90
  SUBROUTINE trace_end(name)
  IMPLICIT NONE
#ifdef VTRACE
#include <vt_user.inc>
#endif 

    CHARACTER(LEN=*),INTENT(IN) :: name

!$OMP MASTER
#ifdef VTRACE
     VT_USER_END(name)
#endif
!$OMP END MASTER

  END SUBROUTINE trace_end    


  !! from transfert_mpi.f90 with _seq version.
  SUBROUTINE wait_message(message)
  IMPLICIT NONE
    TYPE(t_message) :: message
    
  END SUBROUTINE wait_message

  !! from transfert_mpi.f90 with _seq version.
  SUBROUTINE test_message(message)
  IMPLICIT NONE
    TYPE(t_message) :: message
  END SUBROUTINE  test_message

end module mod_misc

!!$ Local Variables:
!!$ f90-do-indent: 2
!!$ f90-if-indent: 2
!!$ f90-program-indent: 2
!!$ f90-type-indent: 2
!!$ f90-directive-comment-re: "!\$omp "
!!$ end:

