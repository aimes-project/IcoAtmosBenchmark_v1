!-------------------------------------------------------------------------------
!
!+  Program dynamics kernel driver (flux calculation for horizontal tracer advection)
!
!-------------------------------------------------------------------------------
program dyn_horiz_adv_flux
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_misc
  use mod_src_tracer, only: &
     horizontal_flux
  !-----------------------------------------------------------------------------
  implicit none

  real(DP), allocatable :: ORG_flx_h       (:,:,:,:)
  real(DP), allocatable :: ORG_flx_h_pl    (:,:,:)
  real(DP), allocatable :: ORG_GRD_xc      (:,:,:,:,:)
  real(DP), allocatable :: ORG_GRD_xc_pl   (:,:,:,:)
  real(DP), allocatable :: ORG_rhog_mean   (:,:,:)
  real(DP), allocatable :: ORG_rhog_mean_pl(:,:,:)
  real(DP), allocatable :: ORG_rhogvx      (:,:,:)
  real(DP), allocatable :: ORG_rhogvx_pl   (:,:,:)
  real(DP), allocatable :: ORG_rhogvy      (:,:,:)
  real(DP), allocatable :: ORG_rhogvy_pl   (:,:,:)
  real(DP), allocatable :: ORG_rhogvz      (:,:,:)
  real(DP), allocatable :: ORG_rhogvz_pl   (:,:,:)
  real(DP), allocatable :: ORG_GRD_xr      (:,:,:,:,:)
  real(DP), allocatable :: ORG_GRD_xr_pl   (:,:,:,:)
  real(DP), allocatable :: ORG_GMTR_p      (:,:,:,:)
  real(DP), allocatable :: ORG_GMTR_p_pl   (:,:,:,:)
  real(DP), allocatable :: ORG_GMTR_t      (:,:,:,:,:)
  real(DP), allocatable :: ORG_GMTR_t_pl   (:,:,:,:)
  real(DP), allocatable :: ORG_GMTR_a      (:,:,:,:,:)
  real(DP), allocatable :: ORG_GMTR_a_pl   (:,:,:,:)

  real(RP), allocatable :: flx_h       (:,:,:,:)
  real(RP), allocatable :: flx_h_pl    (:,:,:)
  real(RP), allocatable :: GRD_xc      (:,:,:,:,:)
  real(RP), allocatable :: GRD_xc_pl   (:,:,:,:)
  real(RP), allocatable :: rhog_mean   (:,:,:)
  real(RP), allocatable :: rhog_mean_pl(:,:,:)
  real(RP), allocatable :: rhogvx      (:,:,:)
  real(RP), allocatable :: rhogvx_pl   (:,:,:)
  real(RP), allocatable :: rhogvy      (:,:,:)
  real(RP), allocatable :: rhogvy_pl   (:,:,:)
  real(RP), allocatable :: rhogvz      (:,:,:)
  real(RP), allocatable :: rhogvz_pl   (:,:,:)

  real(RP), allocatable :: check_flx_h    (:,:,:,:)
  real(RP), allocatable :: check_flx_h_pl (:,:,:)
  real(RP), allocatable :: check_GRD_xc   (:,:,:,:,:)
  real(RP), allocatable :: check_GRD_xc_pl(:,:,:,:)

  integer :: iteration
  !=============================================================================

  write(*,*) "[KERNEL] dyn_horiz_adv_flux"

  !###############################################################################

  write(*,*) "*** Start  initialize"

  allocate( flx_h       (ADM_gall   ,ADM_kall,ADM_lall   ,6  )        )
  allocate( flx_h_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl    )        )
  allocate( GRD_xc      (ADM_gall   ,ADM_kall,ADM_lall   ,3,3)        )
  allocate( GRD_xc_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl  ,3)        )
  allocate( rhog_mean   (ADM_gall   ,ADM_kall,ADM_lall   )            )
  allocate( rhog_mean_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)            )
  allocate( rhogvx      (ADM_gall   ,ADM_kall,ADM_lall   )            )
  allocate( rhogvx_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)            )
  allocate( rhogvy      (ADM_gall   ,ADM_kall,ADM_lall   )            )
  allocate( rhogvy_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)            )
  allocate( rhogvz      (ADM_gall   ,ADM_kall,ADM_lall   )            )
  allocate( rhogvz_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)            )
  flx_h       (:,:,:,:)   = 0.0_RP
  flx_h_pl    (:,:,:)     = 0.0_RP
  GRD_xc      (:,:,:,:,:) = 0.0_RP
  GRD_xc_pl   (:,:,:,:)   = 0.0_RP
  rhog_mean   (:,:,:)     = 0.0_RP
  rhog_mean_pl(:,:,:)     = 0.0_RP
  rhogvx      (:,:,:)     = 0.0_RP
  rhogvx_pl   (:,:,:)     = 0.0_RP
  rhogvy      (:,:,:)     = 0.0_RP
  rhogvy_pl   (:,:,:)     = 0.0_RP
  rhogvz      (:,:,:)     = 0.0_RP
  rhogvz_pl   (:,:,:)     = 0.0_RP

  allocate( GRD_xr      (ADM_gall   ,K0,ADM_lall   ,3,3)              )
  allocate( GRD_xr_pl   (ADM_gall_pl,K0,ADM_lall_pl,  3)              )
  allocate( GMTR_p      (ADM_gall   ,K0,ADM_lall   ,  GMTR_p_nmax   ) )
  allocate( GMTR_p_pl   (ADM_gall_pl,K0,ADM_lall_pl,  GMTR_p_nmax   ) )
  allocate( GMTR_t      (ADM_gall   ,K0,ADM_lall   ,2,GMTR_t_nmax   ) )
  allocate( GMTR_t_pl   (ADM_gall_pl,K0,ADM_lall_pl,  GMTR_t_nmax   ) )
  allocate( GMTR_a      (ADM_gall   ,K0,ADM_lall   ,3,GMTR_a_nmax   ) )
  allocate( GMTR_a_pl   (ADM_gall_pl,K0,ADM_lall_pl,  GMTR_a_nmax_pl) )
  GRD_xr      (:,:,:,:,:) = 0.0_RP
  GRD_xr_pl   (:,:,:,:)   = 0.0_RP
  GMTR_p      (:,:,:,:)   = 0.0_RP
  GMTR_p_pl   (:,:,:,:)   = 0.0_RP
  GMTR_t      (:,:,:,:,:) = 0.0_RP
  GMTR_t_pl   (:,:,:,:)   = 0.0_RP
  GMTR_a      (:,:,:,:,:) = 0.0_RP
  GMTR_a_pl   (:,:,:,:)   = 0.0_RP

  allocate( check_flx_h    (ADM_gall   ,ADM_kall,ADM_lall   ,6  ) )
  allocate( check_flx_h_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl    ) )
  allocate( check_GRD_xc   (ADM_gall   ,ADM_kall,ADM_lall   ,3,3) )
  allocate( check_GRD_xc_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl  ,3) )
  check_flx_h    (:,:,:,:)   = 0.0_RP
  check_flx_h_pl (:,:,:)     = 0.0_RP
  check_GRD_xc   (:,:,:,:,:) = 0.0_RP
  check_GRD_xc_pl(:,:,:,:)   = 0.0_RP

  !###############################################################################

  !---< read input data >---
  allocate( ORG_flx_h       (ADM_gall   ,ADM_kall,ADM_lall   ,6  )        )
  allocate( ORG_flx_h_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl    )        )
  allocate( ORG_GRD_xc      (ADM_gall   ,ADM_kall,ADM_lall   ,3,3)        )
  allocate( ORG_GRD_xc_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl  ,3)        )
  allocate( ORG_rhog_mean   (ADM_gall   ,ADM_kall,ADM_lall   )            )
  allocate( ORG_rhog_mean_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)            )
  allocate( ORG_rhogvx      (ADM_gall   ,ADM_kall,ADM_lall   )            )
  allocate( ORG_rhogvx_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)            )
  allocate( ORG_rhogvy      (ADM_gall   ,ADM_kall,ADM_lall   )            )
  allocate( ORG_rhogvy_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)            )
  allocate( ORG_rhogvz      (ADM_gall   ,ADM_kall,ADM_lall   )            )
  allocate( ORG_rhogvz_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)            )
  allocate( ORG_GRD_xr      (ADM_gall   ,K0,ADM_lall   ,3,3)              )
  allocate( ORG_GRD_xr_pl   (ADM_gall_pl,K0,ADM_lall_pl,  3)              )
  allocate( ORG_GMTR_p      (ADM_gall   ,K0,ADM_lall   ,  GMTR_p_nmax   ) )
  allocate( ORG_GMTR_p_pl   (ADM_gall_pl,K0,ADM_lall_pl,  GMTR_p_nmax   ) )
  allocate( ORG_GMTR_t      (ADM_gall   ,K0,ADM_lall   ,2,GMTR_t_nmax   ) )
  allocate( ORG_GMTR_t_pl   (ADM_gall_pl,K0,ADM_lall_pl,  GMTR_t_nmax   ) )
  allocate( ORG_GMTR_a      (ADM_gall   ,K0,ADM_lall   ,3,GMTR_a_nmax   ) )
  allocate( ORG_GMTR_a_pl   (ADM_gall_pl,K0,ADM_lall_pl,  GMTR_a_nmax_pl) )

  call dumpio_syscheck
  call dumpio_mk_fname(EX_fname,'snapshot.dyn_horiz_adv_flux','pe',SET_prc_me-1,6)
  call dumpio_fopen(EX_fid,EX_fname,IO_FREAD)

  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*ADM_lall   *6         , ORG_flx_h       (:,:,:,:)   )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*ADM_lall_pl           , ORG_flx_h_pl    (:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*ADM_lall   *3*3       , ORG_GRD_xc      (:,:,:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*ADM_lall_pl*  3       , ORG_GRD_xc_pl   (:,:,:,:)   )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*ADM_lall              , ORG_rhog_mean   (:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*ADM_lall_pl           , ORG_rhog_mean_pl(:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*ADM_lall              , ORG_rhogvx      (:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*ADM_lall_pl           , ORG_rhogvx_pl   (:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*ADM_lall              , ORG_rhogvy      (:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*ADM_lall_pl           , ORG_rhogvy_pl   (:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*ADM_lall              , ORG_rhogvz      (:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*ADM_lall_pl           , ORG_rhogvz_pl   (:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_gall   *K0*ADM_lall   *3*3             , ORG_GRD_xr      (:,:,:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall_pl*K0*ADM_lall_pl*  3             , ORG_GRD_xr_pl   (:,:,:,:)   )
  call dumpio_read_data( EX_fid, ADM_gall   *K0*ADM_lall   *  GMTR_p_nmax   , ORG_GMTR_p      (:,:,:,:)   )
  call dumpio_read_data( EX_fid, ADM_gall_pl*K0*ADM_lall_pl*  GMTR_p_nmax   , ORG_GMTR_p_pl   (:,:,:,:)   )
  call dumpio_read_data( EX_fid, ADM_gall   *K0*ADM_lall   *2*GMTR_t_nmax   , ORG_GMTR_t      (:,:,:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall_pl*K0*ADM_lall_pl*  GMTR_t_nmax   , ORG_GMTR_t_pl   (:,:,:,:)   )
  call dumpio_read_data( EX_fid, ADM_gall   *K0*ADM_lall   *3*GMTR_a_nmax   , ORG_GMTR_a      (:,:,:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall_pl*K0*ADM_lall_pl*  GMTR_a_nmax_pl, ORG_GMTR_a_pl   (:,:,:,:)   )

  call dumpio_fclose(EX_fid)

  flx_h       (:,:,:,:)   = real( ORG_flx_h       (:,:,:,:)  , kind=RP )
  flx_h_pl    (:,:,:)     = real( ORG_flx_h_pl    (:,:,:)    , kind=RP )
  GRD_xc      (:,:,:,:,:) = real( ORG_GRD_xc      (:,:,:,:,:), kind=RP )
  GRD_xc_pl   (:,:,:,:)   = real( ORG_GRD_xc_pl   (:,:,:,:)  , kind=RP )
  rhog_mean   (:,:,:)     = real( ORG_rhog_mean   (:,:,:)    , kind=RP )
  rhog_mean_pl(:,:,:)     = real( ORG_rhog_mean_pl(:,:,:)    , kind=RP )
  rhogvx      (:,:,:)     = real( ORG_rhogvx      (:,:,:)    , kind=RP )
  rhogvx_pl   (:,:,:)     = real( ORG_rhogvx_pl   (:,:,:)    , kind=RP )
  rhogvy      (:,:,:)     = real( ORG_rhogvy      (:,:,:)    , kind=RP )
  rhogvy_pl   (:,:,:)     = real( ORG_rhogvy_pl   (:,:,:)    , kind=RP )
  rhogvz      (:,:,:)     = real( ORG_rhogvz      (:,:,:)    , kind=RP )
  rhogvz_pl   (:,:,:)     = real( ORG_rhogvz_pl   (:,:,:)    , kind=RP )
  GRD_xr      (:,:,:,:,:) = real( ORG_GRD_xr      (:,:,:,:,:), kind=RP )
  GRD_xr_pl   (:,:,:,:)   = real( ORG_GRD_xr_pl   (:,:,:,:)  , kind=RP )
  GMTR_p      (:,:,:,:)   = real( ORG_GMTR_p      (:,:,:,:)  , kind=RP )
  GMTR_p_pl   (:,:,:,:)   = real( ORG_GMTR_p_pl   (:,:,:,:)  , kind=RP )
  GMTR_t      (:,:,:,:,:) = real( ORG_GMTR_t      (:,:,:,:,:), kind=RP )
  GMTR_t_pl   (:,:,:,:)   = real( ORG_GMTR_t_pl   (:,:,:,:)  , kind=RP )
  GMTR_a      (:,:,:,:,:) = real( ORG_GMTR_a      (:,:,:,:,:), kind=RP )
  GMTR_a_pl   (:,:,:,:)   = real( ORG_GMTR_a_pl   (:,:,:,:)  , kind=RP )

  deallocate( ORG_flx_h        )
  deallocate( ORG_flx_h_pl     )
  deallocate( ORG_GRD_xc       )
  deallocate( ORG_GRD_xc_pl    )
  deallocate( ORG_rhog_mean    )
  deallocate( ORG_rhog_mean_pl )
  deallocate( ORG_rhogvx       )
  deallocate( ORG_rhogvx_pl    )
  deallocate( ORG_rhogvy       )
  deallocate( ORG_rhogvy_pl    )
  deallocate( ORG_rhogvz       )
  deallocate( ORG_rhogvz_pl    )
  deallocate( ORG_GRD_xr       )
  deallocate( ORG_GRD_xr_pl    )
  deallocate( ORG_GMTR_p       )
  deallocate( ORG_GMTR_p_pl    )
  deallocate( ORG_GMTR_t       )
  deallocate( ORG_GMTR_t_pl    )
  deallocate( ORG_GMTR_a       )
  deallocate( ORG_GMTR_a_pl    )

  !###############################################################################

  call CONST_setup

  check_flx_h    (:,:,:,:)   = flx_h    (:,:,:,:)
  check_flx_h_pl (:,:,:)     = flx_h_pl (:,:,:)
  check_GRD_xc   (:,:,:,:,:) = GRD_xc   (:,:,:,:,:)
  check_GRD_xc_pl(:,:,:,:)   = GRD_xc_pl(:,:,:,:)

  write(*,*) "*** Finish initialize"

  !###############################################################################

  write(*,*) "*** Start  kernel"

  do iteration = 1, SET_iteration

     call DEBUG_rapstart('MAIN_dyn_horiz_adv_flux')

     call horizontal_flux( flx_h    (:,:,:,:),   flx_h_pl    (:,:,:),   & ! [OUT]
                           GRD_xc   (:,:,:,:,:), GRD_xc_pl   (:,:,:,:), & ! [OUT]
                           rhog_mean(:,:,:),     rhog_mean_pl(:,:,:),   & ! [IN]
                           rhogvx   (:,:,:),     rhogvx_pl   (:,:,:),   & ! [IN]
                           rhogvy   (:,:,:),     rhogvy_pl   (:,:,:),   & ! [IN]
                           rhogvz   (:,:,:),     rhogvz_pl   (:,:,:),   & ! [IN]
                           SET_dt_large                                 ) ! [IN]

     call DEBUG_rapend  ('MAIN_dyn_horiz_adv_flux')

     write(ADM_LOG_FID,*) '### Input ###'
     call DEBUG_valuecheck( 'flx_h       ', check_flx_h    (:,:,:,:)   )
     call DEBUG_valuecheck( 'flx_h_pl    ', check_flx_h_pl (:,:,:)     )
     call DEBUG_valuecheck( 'GRD_xc      ', check_GRD_xc   (:,:,:,:,:) )
     call DEBUG_valuecheck( 'GRD_xc_pl   ', check_GRD_xc_pl(:,:,:,:)   )
     call DEBUG_valuecheck( 'rhog_mean   ', rhog_mean      (:,:,:)     )
     call DEBUG_valuecheck( 'rhog_mean_pl', rhog_mean_pl   (:,:,:)     )
     call DEBUG_valuecheck( 'rhogvx      ', rhogvx         (:,:,:)     )
     call DEBUG_valuecheck( 'rhogvx_pl   ', rhogvx_pl      (:,:,:)     )
     call DEBUG_valuecheck( 'rhogvy      ', rhogvy         (:,:,:)     )
     call DEBUG_valuecheck( 'rhogvy_pl   ', rhogvy_pl      (:,:,:)     )
     call DEBUG_valuecheck( 'rhogvz      ', rhogvz         (:,:,:)     )
     call DEBUG_valuecheck( 'rhogvz_pl   ', rhogvz_pl      (:,:,:)     )
     call DEBUG_valuecheck( 'GRD_xr      ', GRD_xr         (:,:,:,:,:) )
     call DEBUG_valuecheck( 'GRD_xr_pl   ', GRD_xr_pl      (:,:,:,:)   )
     call DEBUG_valuecheck( 'GMTR_p      ', GMTR_p         (:,:,:,:)   )
     call DEBUG_valuecheck( 'GMTR_p_pl   ', GMTR_p_pl      (:,:,:,:)   )
     call DEBUG_valuecheck( 'GMTR_t      ', GMTR_t         (:,:,:,:,:) )
     call DEBUG_valuecheck( 'GMTR_t_pl   ', GMTR_t_pl      (:,:,:,:)   )
     call DEBUG_valuecheck( 'GMTR_a      ', GMTR_a         (:,:,:,:,:) )
     call DEBUG_valuecheck( 'GMTR_a_pl   ', GMTR_a_pl      (:,:,:,:)   )
     write(ADM_LOG_FID,*) '### Output ###'
     call DEBUG_valuecheck( 'flx_h       ', flx_h          (:,:,:,:)   )
     call DEBUG_valuecheck( 'flx_h_pl    ', flx_h_pl       (:,:,:)     )
     call DEBUG_valuecheck( 'GRD_xc      ', GRD_xc         (:,:,:,:,:) )
     call DEBUG_valuecheck( 'GRD_xc_pl   ', GRD_xc_pl      (:,:,:,:)   )
  enddo

  write(ADM_LOG_FID,*) '### Validation : grid-by-grid diff ###'
  check_flx_h    (:,:,:,:)   = check_flx_h    (:,:,:,:)   - flx_h    (:,:,:,:)
  check_flx_h_pl (:,:,:)     = check_flx_h_pl (:,:,:)     - flx_h_pl (:,:,:)
  check_GRD_xc   (:,:,:,:,:) = check_GRD_xc   (:,:,:,:,:) - GRD_xc   (:,:,:,:,:)
  check_GRD_xc_pl(:,:,:,:)   = check_GRD_xc_pl(:,:,:,:)   - GRD_xc_pl(:,:,:,:)
  call DEBUG_valuecheck( 'check_flx_h    ', check_flx_h    (:,:,:,:)   )
  call DEBUG_valuecheck( 'check_flx_h_pl ', check_flx_h_pl (:,:,:)     )
  call DEBUG_valuecheck( 'check_GRD_xc   ', check_GRD_xc   (:,:,:,:,:) )
  call DEBUG_valuecheck( 'check_GRD_xc_pl', check_GRD_xc_pl(:,:,:,:)   )

  write(*,*) "*** Finish kernel"

  !###############################################################################

  call DEBUG_rapreport

  stop
end program dyn_horiz_adv_flux
!-------------------------------------------------------------------------------
