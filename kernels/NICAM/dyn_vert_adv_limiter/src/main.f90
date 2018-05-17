!-------------------------------------------------------------------------------
!
!+  Program dynamics kernel driver (flux limiter for vertical tracer advection)
!
!-------------------------------------------------------------------------------
program dyn_vert_adv_limiter
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_misc
  use mod_src_tracer, only: &
     vertical_limiter_thuburn
  !-----------------------------------------------------------------------------
  implicit none

  real(DP), allocatable :: ORG_q_h_prev   (:,:,:)
  real(DP), allocatable :: ORG_q_h_prev_pl(:,:,:)
  real(DP), allocatable :: ORG_q_h        (:,:,:)
  real(DP), allocatable :: ORG_q_h_pl     (:,:,:)
  real(DP), allocatable :: ORG_q          (:,:,:)
  real(DP), allocatable :: ORG_q_pl       (:,:,:)
  real(DP), allocatable :: ORG_d          (:,:,:)
  real(DP), allocatable :: ORG_d_pl       (:,:,:)
  real(DP), allocatable :: ORG_ck         (:,:,:,:)
  real(DP), allocatable :: ORG_ck_pl      (:,:,:,:)

  real(RP), allocatable :: q_h_prev   (:,:,:)
  real(RP), allocatable :: q_h_prev_pl(:,:,:)
  real(RP), allocatable :: q_h        (:,:,:)
  real(RP), allocatable :: q_h_pl     (:,:,:)
  real(RP), allocatable :: q          (:,:,:)
  real(RP), allocatable :: q_pl       (:,:,:)
  real(RP), allocatable :: d          (:,:,:)
  real(RP), allocatable :: d_pl       (:,:,:)
  real(RP), allocatable :: ck         (:,:,:,:)
  real(RP), allocatable :: ck_pl      (:,:,:,:)

  real(RP), allocatable :: check_q_h    (:,:,:)
  real(RP), allocatable :: check_q_h_pl (:,:,:)

  integer :: iteration
  !=============================================================================

  write(*,*) "[KERNEL] dyn_vert_adv_limiter"

  !###############################################################################

  write(*,*) "*** Start  initialize"

  allocate( q_h_prev   (ADM_gall   ,ADM_kall,ADM_lall   )   )
  allocate( q_h_prev_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)   )
  allocate( q_h        (ADM_gall   ,ADM_kall,ADM_lall   )   )
  allocate( q_h_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)   )
  allocate( q          (ADM_gall   ,ADM_kall,ADM_lall   )   )
  allocate( q_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl)   )
  allocate( d          (ADM_gall   ,ADM_kall,ADM_lall   )   )
  allocate( d_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl)   )
  allocate( ck         (ADM_gall   ,ADM_kall,ADM_lall   ,2) )
  allocate( ck_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl,2) )
  q_h_prev   (:,:,:)   = 0.0_RP
  q_h_prev_pl(:,:,:)   = 0.0_RP
  q_h        (:,:,:)   = 0.0_RP
  q_h_pl     (:,:,:)   = 0.0_RP
  q          (:,:,:)   = 0.0_RP
  q_pl       (:,:,:)   = 0.0_RP
  d          (:,:,:)   = 0.0_RP
  d_pl       (:,:,:)   = 0.0_RP
  ck         (:,:,:,:) = 0.0_RP
  ck_pl      (:,:,:,:) = 0.0_RP

  allocate( check_q_h   (ADM_gall   ,ADM_kall,ADM_lall   ) )
  allocate( check_q_h_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
  check_q_h   (:,:,:) = 0.0_RP
  check_q_h_pl(:,:,:) = 0.0_RP

  !###############################################################################

  !---< read input data >---
  allocate( ORG_q_h_prev   (ADM_gall   ,ADM_kall,ADM_lall   )   )
  allocate( ORG_q_h_prev_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)   )
  allocate( ORG_q_h        (ADM_gall   ,ADM_kall,ADM_lall   )   )
  allocate( ORG_q_h_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)   )
  allocate( ORG_q          (ADM_gall   ,ADM_kall,ADM_lall   )   )
  allocate( ORG_q_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl)   )
  allocate( ORG_d          (ADM_gall   ,ADM_kall,ADM_lall   )   )
  allocate( ORG_d_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl)   )
  allocate( ORG_ck         (ADM_gall   ,ADM_kall,ADM_lall   ,2) )
  allocate( ORG_ck_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl,2) )

  call dumpio_syscheck
  call dumpio_mk_fname(EX_fname,'snapshot.dyn_vert_adv_limiter','pe',SET_prc_me-1,6)
  call dumpio_fopen(EX_fid,EX_fname,IO_FREAD)

  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*ADM_lall     , ORG_q_h_prev   (:,:,:)   )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*ADM_lall_pl  , ORG_q_h_prev_pl(:,:,:)   )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*ADM_lall     , ORG_q_h        (:,:,:)   )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*ADM_lall_pl  , ORG_q_h_pl     (:,:,:)   )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*ADM_lall     , ORG_q          (:,:,:)   )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*ADM_lall_pl  , ORG_q_pl       (:,:,:)   )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*ADM_lall     , ORG_d          (:,:,:)   )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*ADM_lall_pl  , ORG_d_pl       (:,:,:)   )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*ADM_lall   *2, ORG_ck         (:,:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*ADM_lall_pl*2, ORG_ck_pl      (:,:,:,:) )

  call dumpio_fclose(EX_fid)

  q_h_prev   (:,:,:)   = real( ORG_q_h_prev   (:,:,:)  , kind=RP )
  q_h_prev_pl(:,:,:)   = real( ORG_q_h_prev_pl(:,:,:)  , kind=RP )
  q_h        (:,:,:)   = real( ORG_q_h        (:,:,:)  , kind=RP )
  q_h_pl     (:,:,:)   = real( ORG_q_h_pl     (:,:,:)  , kind=RP )
  q          (:,:,:)   = real( ORG_q          (:,:,:)  , kind=RP )
  q_pl       (:,:,:)   = real( ORG_q_pl       (:,:,:)  , kind=RP )
  d          (:,:,:)   = real( ORG_d          (:,:,:)  , kind=RP )
  d_pl       (:,:,:)   = real( ORG_d_pl       (:,:,:)  , kind=RP )
  ck         (:,:,:,:) = real( ORG_ck         (:,:,:,:), kind=RP )
  ck_pl      (:,:,:,:) = real( ORG_ck_pl      (:,:,:,:), kind=RP )

  deallocate( ORG_q_h_prev    )
  deallocate( ORG_q_h_prev_pl )
  deallocate( ORG_q_h         )
  deallocate( ORG_q_h_pl      )
  deallocate( ORG_q           )
  deallocate( ORG_q_pl        )
  deallocate( ORG_d           )
  deallocate( ORG_d_pl        )
  deallocate( ORG_ck          )
  deallocate( ORG_ck_pl       )

  !###############################################################################

  call CONST_setup

  check_q_h   (:,:,:) = q_h   (:,:,:)
  check_q_h_pl(:,:,:) = q_h_pl(:,:,:)

  write(*,*) "*** Finish initialize"

  !###############################################################################

  write(*,*) "*** Start  kernel"

  do iteration = 1, SET_iteration

     call DEBUG_rapstart('MAIN_dyn_vert_adv_limiter')

     ! restore previous value for multiple iteration
     q_h   (:,:,:) = q_h_prev   (:,:,:)
     q_h_pl(:,:,:) = q_h_prev_pl(:,:,:)

     call vertical_limiter_thuburn( q_h(:,:,:),   q_h_pl(:,:,:),  & ! [INOUT]
                                    q  (:,:,:),   q_pl  (:,:,:),  & ! [IN]
                                    d  (:,:,:),   d_pl  (:,:,:),  & ! [IN]
                                    ck (:,:,:,:), ck_pl (:,:,:,:) ) ! [IN]

     call DEBUG_rapend  ('MAIN_dyn_vert_adv_limiter')

     write(ADM_LOG_FID,*) '### Input ###'
     call DEBUG_valuecheck( 'q_h_prev    ', q_h_prev    (:,:,:)   )
     call DEBUG_valuecheck( 'q_h_prev_pl ', q_h_prev_pl (:,:,:)   )
     call DEBUG_valuecheck( 'check_q_h   ', check_q_h   (:,:,:)   )
     call DEBUG_valuecheck( 'check_q_h_pl', check_q_h_pl(:,:,:)   )
     call DEBUG_valuecheck( 'q           ', q           (:,:,:)   )
     call DEBUG_valuecheck( 'q_pl        ', q_pl        (:,:,:)   )
     call DEBUG_valuecheck( 'd           ', d           (:,:,:)   )
     call DEBUG_valuecheck( 'd_pl        ', d_pl        (:,:,:)   )
     call DEBUG_valuecheck( 'ck          ', ck          (:,:,:,:) )
     call DEBUG_valuecheck( 'ck_pl       ', ck_pl       (:,:,:,:) )
     write(ADM_LOG_FID,*) '### Output ###'
     call DEBUG_valuecheck( 'q_h         ', q_h         (:,:,:)   )
     call DEBUG_valuecheck( 'q_h_pl      ', q_h_pl      (:,:,:)   )
  enddo

  write(ADM_LOG_FID,*) '### Validation : point-by-point diff ###'
  check_q_h   (:,:,:) = check_q_h   (:,:,:) - q_h   (:,:,:)
  check_q_h_pl(:,:,:) = check_q_h_pl(:,:,:) - q_h_pl(:,:,:)
  call DEBUG_valuecheck( 'check_q_h   ', check_q_h   (:,:,:)   )
  call DEBUG_valuecheck( 'check_q_h_pl', check_q_h_pl(:,:,:)   )

  write(*,*) "*** Finish kernel"

  !###############################################################################

  call DEBUG_rapreport

  stop
end program dyn_vert_adv_limiter
!-------------------------------------------------------------------------------
