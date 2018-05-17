!-------------------------------------------------------------------------------
!
!+  Program dynamics kernel driver (flux limiter for horizontal tracer advection)
!
!-------------------------------------------------------------------------------
program dyn_horiz_adv_limiter
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_misc
  use mod_src_tracer, only: &
     horizontal_limiter_thuburn
  !-----------------------------------------------------------------------------
  implicit none

  real(DP), allocatable :: ORG_q_a_prev    (:,:,:,:)
  real(DP), allocatable :: ORG_q_a_prev_pl (:,:,:)
  real(DP), allocatable :: ORG_q_a         (:,:,:,:)
  real(DP), allocatable :: ORG_q_a_pl      (:,:,:)
  real(DP), allocatable :: ORG_q           (:,:,:)
  real(DP), allocatable :: ORG_q_pl        (:,:,:)
  real(DP), allocatable :: ORG_d           (:,:,:)
  real(DP), allocatable :: ORG_d_pl        (:,:,:)
  real(DP), allocatable :: ORG_ch          (:,:,:,:)
  real(DP), allocatable :: ORG_ch_pl       (:,:,:)
  real(DP), allocatable :: ORG_cmask       (:,:,:,:)
  real(DP), allocatable :: ORG_cmask_pl    (:,:,:)
  real(DP), allocatable :: ORG_Qout_prev   (:,:,:,:)
  real(DP), allocatable :: ORG_Qout_prev_pl(:,:,:,:)
  real(DP), allocatable :: ORG_Qout_post   (:,:,:,:)
  real(DP), allocatable :: ORG_Qout_post_pl(:,:,:,:)

  real(RP), allocatable :: q_a_prev    (:,:,:,:)
  real(RP), allocatable :: q_a_prev_pl (:,:,:)
  real(RP), allocatable :: q_a         (:,:,:,:)
  real(RP), allocatable :: q_a_pl      (:,:,:)
  real(RP), allocatable :: q           (:,:,:)
  real(RP), allocatable :: q_pl        (:,:,:)
  real(RP), allocatable :: d           (:,:,:)
  real(RP), allocatable :: d_pl        (:,:,:)
  real(RP), allocatable :: ch          (:,:,:,:)
  real(RP), allocatable :: ch_pl       (:,:,:)
  real(RP), allocatable :: cmask       (:,:,:,:)
  real(RP), allocatable :: cmask_pl    (:,:,:)
  real(RP), allocatable :: Qout_prev   (:,:,:,:)
  real(RP), allocatable :: Qout_prev_pl(:,:,:,:)
  real(RP), allocatable :: Qout_post   (:,:,:,:)
  real(RP), allocatable :: Qout_post_pl(:,:,:,:)

  real(RP), allocatable :: check_q_a         (:,:,:,:)
  real(RP), allocatable :: check_q_a_pl      (:,:,:)
  real(RP), allocatable :: check_Qout_prev   (:,:,:,:)
  real(RP), allocatable :: check_Qout_prev_pl(:,:,:,:)

  integer :: iteration
  !=============================================================================

  write(*,*) "[KERNEL] dyn_horiz_adv_limiter"

  !###############################################################################

  write(*,*) "*** Start  initialize"

  allocate( q_a_prev    (ADM_gall   ,ADM_kall,ADM_lall   ,6) )
  allocate( q_a_prev_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl  ) )
  allocate( q_a         (ADM_gall   ,ADM_kall,ADM_lall   ,6) )
  allocate( q_a_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl  ) )
  allocate( q           (ADM_gall   ,ADM_kall,ADM_lall     ) )
  allocate( q_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl  ) )
  allocate( d           (ADM_gall   ,ADM_kall,ADM_lall     ) )
  allocate( d_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl  ) )
  allocate( ch          (ADM_gall   ,ADM_kall,ADM_lall   ,6) )
  allocate( ch_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl  ) )
  allocate( cmask       (ADM_gall   ,ADM_kall,ADM_lall   ,6) )
  allocate( cmask_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl  ) )
  allocate( Qout_prev   (ADM_gall   ,ADM_kall,ADM_lall   ,2) )
  allocate( Qout_prev_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,2) )
  allocate( Qout_post   (ADM_gall   ,ADM_kall,ADM_lall   ,2) )
  allocate( Qout_post_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,2) )
  q_a_prev    (:,:,:,:) = 0.0_RP
  q_a_prev_pl (:,:,:)   = 0.0_RP
  q_a         (:,:,:,:) = 0.0_RP
  q_a_pl      (:,:,:)   = 0.0_RP
  q           (:,:,:)   = 0.0_RP
  q_pl        (:,:,:)   = 0.0_RP
  d           (:,:,:)   = 0.0_RP
  d_pl        (:,:,:)   = 0.0_RP
  ch          (:,:,:,:) = 0.0_RP
  ch_pl       (:,:,:)   = 0.0_RP
  cmask       (:,:,:,:) = 0.0_RP
  cmask_pl    (:,:,:)   = 0.0_RP
  Qout_prev   (:,:,:,:) = 0.0_RP
  Qout_prev_pl(:,:,:,:) = 0.0_RP
  Qout_post   (:,:,:,:) = 0.0_RP
  Qout_post_pl(:,:,:,:) = 0.0_RP

  allocate( check_q_a         (ADM_gall   ,ADM_kall,ADM_lall   ,6) )
  allocate( check_q_a_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl  ) )
  allocate( check_Qout_prev   (ADM_gall   ,ADM_kall,ADM_lall   ,2) )
  allocate( check_Qout_prev_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,2) )
  check_q_a         (:,:,:,:) = 0.0_RP
  check_q_a_pl      (:,:,:)   = 0.0_RP
  check_Qout_prev   (:,:,:,:) = 0.0_RP
  check_Qout_prev_pl(:,:,:,:) = 0.0_RP

  !###############################################################################

  !---< read input data >---
  allocate( ORG_q_a_prev    (ADM_gall   ,ADM_kall,ADM_lall   ,6) )
  allocate( ORG_q_a_prev_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl  ) )
  allocate( ORG_q_a         (ADM_gall   ,ADM_kall,ADM_lall   ,6) )
  allocate( ORG_q_a_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl  ) )
  allocate( ORG_q           (ADM_gall   ,ADM_kall,ADM_lall     ) )
  allocate( ORG_q_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl  ) )
  allocate( ORG_d           (ADM_gall   ,ADM_kall,ADM_lall     ) )
  allocate( ORG_d_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl  ) )
  allocate( ORG_ch          (ADM_gall   ,ADM_kall,ADM_lall   ,6) )
  allocate( ORG_ch_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl  ) )
  allocate( ORG_cmask       (ADM_gall   ,ADM_kall,ADM_lall   ,6) )
  allocate( ORG_cmask_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl  ) )
  allocate( ORG_Qout_prev   (ADM_gall   ,ADM_kall,ADM_lall   ,2) )
  allocate( ORG_Qout_prev_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,2) )
  allocate( ORG_Qout_post   (ADM_gall   ,ADM_kall,ADM_lall   ,2) )
  allocate( ORG_Qout_post_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,2) )

  call dumpio_syscheck
  call dumpio_mk_fname(EX_fname,'snapshot.dyn_horiz_adv_limiter','pe',SET_prc_me-1,6)
  call dumpio_fopen(EX_fid,EX_fname,IO_FREAD)

  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*ADM_lall   *6, ORG_q_a_prev    (:,:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*ADM_lall_pl  , ORG_q_a_prev_pl (:,:,:)   )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*ADM_lall   *6, ORG_q_a         (:,:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*ADM_lall_pl  , ORG_q_a_pl      (:,:,:)   )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*ADM_lall     , ORG_q           (:,:,:)   )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*ADM_lall_pl  , ORG_q_pl        (:,:,:)   )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*ADM_lall     , ORG_d           (:,:,:)   )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*ADM_lall_pl  , ORG_d_pl        (:,:,:)   )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*ADM_lall   *6, ORG_ch          (:,:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*ADM_lall_pl  , ORG_ch_pl       (:,:,:)   )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*ADM_lall   *6, ORG_cmask       (:,:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*ADM_lall_pl  , ORG_cmask_pl    (:,:,:)   )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*ADM_lall   *2, ORG_Qout_prev   (:,:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*ADM_lall_pl*2, ORG_Qout_prev_pl(:,:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*ADM_lall   *2, ORG_Qout_post   (:,:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*ADM_lall_pl*2, ORG_Qout_post_pl(:,:,:,:) )

  call dumpio_fclose(EX_fid)

  q_a_prev    (:,:,:,:) = real( ORG_q_a_prev    (:,:,:,:), kind=RP )
  q_a_prev_pl (:,:,:)   = real( ORG_q_a_prev_pl (:,:,:)  , kind=RP )
  q_a         (:,:,:,:) = real( ORG_q_a         (:,:,:,:), kind=RP )
  q_a_pl      (:,:,:)   = real( ORG_q_a_pl      (:,:,:)  , kind=RP )
  q           (:,:,:)   = real( ORG_q           (:,:,:)  , kind=RP )
  q_pl        (:,:,:)   = real( ORG_q_pl        (:,:,:)  , kind=RP )
  d           (:,:,:)   = real( ORG_d           (:,:,:)  , kind=RP )
  d_pl        (:,:,:)   = real( ORG_d_pl        (:,:,:)  , kind=RP )
  ch          (:,:,:,:) = real( ORG_ch          (:,:,:,:), kind=RP )
  ch_pl       (:,:,:)   = real( ORG_ch_pl       (:,:,:)  , kind=RP )
  cmask       (:,:,:,:) = real( ORG_cmask       (:,:,:,:), kind=RP )
  cmask_pl    (:,:,:)   = real( ORG_cmask_pl    (:,:,:)  , kind=RP )
  Qout_prev   (:,:,:,:) = real( ORG_Qout_prev   (:,:,:,:), kind=RP )
  Qout_prev_pl(:,:,:,:) = real( ORG_Qout_prev_pl(:,:,:,:), kind=RP )
  Qout_post   (:,:,:,:) = real( ORG_Qout_post   (:,:,:,:), kind=RP )
  Qout_post_pl(:,:,:,:) = real( ORG_Qout_post_pl(:,:,:,:), kind=RP )

  deallocate( ORG_q_a_prev     )
  deallocate( ORG_q_a_prev_pl  )
  deallocate( ORG_q_a          )
  deallocate( ORG_q_a_pl       )
  deallocate( ORG_q            )
  deallocate( ORG_q_pl         )
  deallocate( ORG_d            )
  deallocate( ORG_d_pl         )
  deallocate( ORG_ch           )
  deallocate( ORG_ch_pl        )
  deallocate( ORG_cmask        )
  deallocate( ORG_cmask_pl     )
  deallocate( ORG_Qout_prev    )
  deallocate( ORG_Qout_prev_pl )
  deallocate( ORG_Qout_post    )
  deallocate( ORG_Qout_post_pl )

  !###############################################################################

  call CONST_setup

  check_q_a         (:,:,:,:) = q_a         (:,:,:,:)
  check_q_a_pl      (:,:,:)   = q_a_pl      (:,:,:)
  check_Qout_prev   (:,:,:,:) = Qout_prev   (:,:,:,:)
  check_Qout_prev_pl(:,:,:,:) = Qout_prev_pl(:,:,:,:)

  write(*,*) "*** Finish initialize"

  !###############################################################################

  write(*,*) "*** Start  kernel"

  do iteration = 1, SET_iteration

     call DEBUG_rapstart('MAIN_dyn_horiz_adv_limiter')

     ! restore previous value for multiple iteration
     q_a   (:,:,:,:) = q_a_prev   (:,:,:,:)
     q_a_pl(:,:,:)   = q_a_prev_pl(:,:,:)

     call horizontal_limiter_thuburn( q_a      (:,:,:,:), q_a_pl      (:,:,:),   & ! [INOUT]
                                      q        (:,:,:),   q_pl        (:,:,:),   & ! [IN]
                                      d        (:,:,:),   d_pl        (:,:,:),   & ! [IN]
                                      ch       (:,:,:,:), ch_pl       (:,:,:),   & ! [IN]
                                      cmask    (:,:,:,:), cmask_pl    (:,:,:),   & ! [IN]
                                      Qout_prev(:,:,:,:), Qout_prev_pl(:,:,:,:), & ! [OUT]
                                      Qout_post(:,:,:,:), Qout_post_pl(:,:,:,:)  ) ! [IN]

     call DEBUG_rapend  ('MAIN_dyn_horiz_adv_limiter')

     write(ADM_LOG_FID,*) '### Input ###'
     call DEBUG_valuecheck( 'q_a_prev          ', q_a_prev          (:,:,:,:) )
     call DEBUG_valuecheck( 'q_a_prev_pl       ', q_a_prev_pl       (:,:,:)   )
     call DEBUG_valuecheck( 'check_q_a         ', check_q_a         (:,:,:,:) )
     call DEBUG_valuecheck( 'check_q_a_pl      ', check_q_a_pl      (:,:,:)   )
     call DEBUG_valuecheck( 'q                 ', q                 (:,:,:)   )
     call DEBUG_valuecheck( 'q_pl              ', q_pl              (:,:,:)   )
     call DEBUG_valuecheck( 'd                 ', d                 (:,:,:)   )
     call DEBUG_valuecheck( 'd_pl              ', d_pl              (:,:,:)   )
     call DEBUG_valuecheck( 'ch                ', ch                (:,:,:,:) )
     call DEBUG_valuecheck( 'ch_pl             ', ch_pl             (:,:,:)   )
     call DEBUG_valuecheck( 'cmask             ', cmask             (:,:,:,:) )
     call DEBUG_valuecheck( 'cmask_pl          ', cmask_pl          (:,:,:)   )
     call DEBUG_valuecheck( 'check_Qout_prev   ', check_Qout_prev   (:,:,:,:) )
     call DEBUG_valuecheck( 'check_Qout_prev_pl', check_Qout_prev_pl(:,:,:,:) )
     call DEBUG_valuecheck( 'Qout_post         ', Qout_post         (:,:,:,:) )
     call DEBUG_valuecheck( 'Qout_post_pl      ', Qout_post_pl      (:,:,:,:) )
     write(ADM_LOG_FID,*) '### Output ###'
     call DEBUG_valuecheck( 'q_a               ', q_a               (:,:,:,:) )
     call DEBUG_valuecheck( 'q_a_pl            ', q_a_pl            (:,:,:)   )
     call DEBUG_valuecheck( 'Qout_prev         ', Qout_prev         (:,:,:,:) )
     call DEBUG_valuecheck( 'Qout_prev_pl      ', Qout_prev_pl      (:,:,:,:) )
  enddo

  write(ADM_LOG_FID,*) '### Validation : point-by-point diff ###'
  check_q_a         (:,:,:,:) = check_q_a         (:,:,:,:) - q_a         (:,:,:,:)
  check_q_a_pl      (:,:,:)   = check_q_a_pl      (:,:,:)   - q_a_pl      (:,:,:)
  check_Qout_prev   (:,:,:,:) = check_Qout_prev   (:,:,:,:) - Qout_prev   (:,:,:,:)
  check_Qout_prev_pl(:,:,:,:) = check_Qout_prev_pl(:,:,:,:) - Qout_prev_pl(:,:,:,:)
  call DEBUG_valuecheck( 'check_q_a         ', check_q_a         (:,:,:,:) )
  call DEBUG_valuecheck( 'check_q_a_pl      ', check_q_a_pl      (:,:,:)   )
  call DEBUG_valuecheck( 'check_Qout_prev   ', check_Qout_prev   (:,:,:,:) )
  call DEBUG_valuecheck( 'check_Qout_prev_pl', check_Qout_prev_pl(:,:,:,:) )

  write(*,*) "*** Finish kernel"

  !###############################################################################

  call DEBUG_rapreport

  stop
end program dyn_horiz_adv_limiter
!-------------------------------------------------------------------------------
