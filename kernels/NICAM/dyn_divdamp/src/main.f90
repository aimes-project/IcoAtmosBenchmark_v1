!-------------------------------------------------------------------------------
!
!+  Program dynamics kernel driver (divergence damping operator)
!
!-------------------------------------------------------------------------------
program dyn_divdamp
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_misc
  use mod_oprt3d, only: &
     OPRT3D_divdamp
  !-----------------------------------------------------------------------------
  implicit none

  real(DP), allocatable :: ORG_ddivdx           (:,:,:)
  real(DP), allocatable :: ORG_ddivdx_pl        (:,:,:)
  real(DP), allocatable :: ORG_ddivdy           (:,:,:)
  real(DP), allocatable :: ORG_ddivdy_pl        (:,:,:)
  real(DP), allocatable :: ORG_ddivdz           (:,:,:)
  real(DP), allocatable :: ORG_ddivdz_pl        (:,:,:)
  real(DP), allocatable :: ORG_rhogvx           (:,:,:)
  real(DP), allocatable :: ORG_rhogvx_pl        (:,:,:)
  real(DP), allocatable :: ORG_rhogvy           (:,:,:)
  real(DP), allocatable :: ORG_rhogvy_pl        (:,:,:)
  real(DP), allocatable :: ORG_rhogvz           (:,:,:)
  real(DP), allocatable :: ORG_rhogvz_pl        (:,:,:)
  real(DP), allocatable :: ORG_rhogw            (:,:,:)
  real(DP), allocatable :: ORG_rhogw_pl         (:,:,:)
  real(DP), allocatable :: ORG_coef_intp        (:,:,:,:,:)
  real(DP), allocatable :: ORG_coef_intp_pl     (:,:,:,:)
  real(DP), allocatable :: ORG_coef_diff        (:,:,:,:)
  real(DP), allocatable :: ORG_coef_diff_pl     (:,:,:)
  real(DP), allocatable :: ORG_VMTR_RGSQRTH     (:,:,:)
  real(DP), allocatable :: ORG_VMTR_RGSQRTH_pl  (:,:,:)
  real(DP), allocatable :: ORG_VMTR_RGAM        (:,:,:)
  real(DP), allocatable :: ORG_VMTR_RGAM_pl     (:,:,:)
  real(DP), allocatable :: ORG_VMTR_RGAMH       (:,:,:)
  real(DP), allocatable :: ORG_VMTR_RGAMH_pl    (:,:,:)
  real(DP), allocatable :: ORG_VMTR_C2WfactGz   (:,:,:,:)
  real(DP), allocatable :: ORG_VMTR_C2WfactGz_pl(:,:,:,:)

  real(RP), allocatable :: ddivdx      (:,:,:)
  real(RP), allocatable :: ddivdx_pl   (:,:,:)
  real(RP), allocatable :: ddivdy      (:,:,:)
  real(RP), allocatable :: ddivdy_pl   (:,:,:)
  real(RP), allocatable :: ddivdz      (:,:,:)
  real(RP), allocatable :: ddivdz_pl   (:,:,:)
  real(RP), allocatable :: rhogvx      (:,:,:)
  real(RP), allocatable :: rhogvx_pl   (:,:,:)
  real(RP), allocatable :: rhogvy      (:,:,:)
  real(RP), allocatable :: rhogvy_pl   (:,:,:)
  real(RP), allocatable :: rhogvz      (:,:,:)
  real(RP), allocatable :: rhogvz_pl   (:,:,:)
  real(RP), allocatable :: rhogw       (:,:,:)
  real(RP), allocatable :: rhogw_pl    (:,:,:)
  real(RP), allocatable :: coef_intp   (:,:,:,:,:)
  real(RP), allocatable :: coef_intp_pl(:,:,:,:)
  real(RP), allocatable :: coef_diff   (:,:,:,:)
  real(RP), allocatable :: coef_diff_pl(:,:,:)

  real(RP), allocatable :: check_ddivdx   (:,:,:)
  real(RP), allocatable :: check_ddivdx_pl(:,:,:)
  real(RP), allocatable :: check_ddivdy   (:,:,:)
  real(RP), allocatable :: check_ddivdy_pl(:,:,:)
  real(RP), allocatable :: check_ddivdz   (:,:,:)
  real(RP), allocatable :: check_ddivdz_pl(:,:,:)

  integer :: iteration
  !=============================================================================

  write(*,*) "[KERNEL] dyn_divdamp"

  !###############################################################################

  write(*,*) "*** Start  initialize"

  allocate( ddivdx      (ADM_gall   ,ADM_kall,ADM_lall   )           )
  allocate( ddivdx_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)           )
  allocate( ddivdy      (ADM_gall   ,ADM_kall,ADM_lall   )           )
  allocate( ddivdy_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)           )
  allocate( ddivdz      (ADM_gall   ,ADM_kall,ADM_lall   )           )
  allocate( ddivdz_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)           )
  allocate( rhogvx      (ADM_gall   ,ADM_kall,ADM_lall   )           )
  allocate( rhogvx_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)           )
  allocate( rhogvy      (ADM_gall   ,ADM_kall,ADM_lall   )           )
  allocate( rhogvy_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)           )
  allocate( rhogvz      (ADM_gall   ,ADM_kall,ADM_lall   )           )
  allocate( rhogvz_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)           )
  allocate( rhogw       (ADM_gall   ,ADM_kall,ADM_lall   )           )
  allocate( rhogw_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)           )
  allocate( coef_intp   (ADM_gall   ,1:3,ADM_nxyz,TI:TJ,ADM_lall   ) )
  allocate( coef_intp_pl(ADM_gall_pl,1:3,ADM_nxyz,      ADM_lall_pl) )
  allocate( coef_diff   (ADM_gall   ,1:6,ADM_nxyz,      ADM_lall   ) )
  allocate( coef_diff_pl(    1:ADM_vlink,ADM_nxyz,      ADM_lall_pl) )
  ddivdx      (:,:,:)     = 0.0_RP
  ddivdx_pl   (:,:,:)     = 0.0_RP
  ddivdy      (:,:,:)     = 0.0_RP
  ddivdy_pl   (:,:,:)     = 0.0_RP
  ddivdz      (:,:,:)     = 0.0_RP
  ddivdz_pl   (:,:,:)     = 0.0_RP
  rhogvx      (:,:,:)     = 0.0_RP
  rhogvx_pl   (:,:,:)     = 0.0_RP
  rhogvy      (:,:,:)     = 0.0_RP
  rhogvy_pl   (:,:,:)     = 0.0_RP
  rhogvz      (:,:,:)     = 0.0_RP
  rhogvz_pl   (:,:,:)     = 0.0_RP
  rhogw       (:,:,:)     = 0.0_RP
  rhogw_pl    (:,:,:)     = 0.0_RP
  coef_intp   (:,:,:,:,:) = 0.0_RP
  coef_intp_pl(:,:,:,:)   = 0.0_RP
  coef_diff   (:,:,:,:)   = 0.0_RP
  coef_diff_pl(:,:,:)     = 0.0_RP

  allocate( VMTR_RGSQRTH     (ADM_gall   ,ADM_kall,ADM_lall   )   )
  allocate( VMTR_RGSQRTH_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)   )
  allocate( VMTR_RGAM        (ADM_gall   ,ADM_kall,ADM_lall   )   )
  allocate( VMTR_RGAM_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)   )
  allocate( VMTR_RGAMH       (ADM_gall   ,ADM_kall,ADM_lall   )   )
  allocate( VMTR_RGAMH_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)   )
  allocate( VMTR_C2WfactGz   (ADM_gall   ,ADM_kall,6,ADM_lall   ) )
  allocate( VMTR_C2WfactGz_pl(ADM_gall_pl,ADM_kall,6,ADM_lall_pl) )
  VMTR_RGSQRTH     (:,:,:)   = 0.0_RP
  VMTR_RGSQRTH_pl  (:,:,:)   = 0.0_RP
  VMTR_RGAM        (:,:,:)   = 0.0_RP
  VMTR_RGAM_pl     (:,:,:)   = 0.0_RP
  VMTR_RGAMH       (:,:,:)   = 0.0_RP
  VMTR_RGAMH_pl    (:,:,:)   = 0.0_RP
  VMTR_C2WfactGz   (:,:,:,:) = 0.0_RP
  VMTR_C2WfactGz_pl(:,:,:,:) = 0.0_RP

  allocate( check_ddivdx   (ADM_gall   ,ADM_kall,ADM_lall   ) )
  allocate( check_ddivdx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
  allocate( check_ddivdy   (ADM_gall   ,ADM_kall,ADM_lall   ) )
  allocate( check_ddivdy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
  allocate( check_ddivdz   (ADM_gall   ,ADM_kall,ADM_lall   ) )
  allocate( check_ddivdz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
  check_ddivdx   (:,:,:) = 0.0_RP
  check_ddivdx_pl(:,:,:) = 0.0_RP
  check_ddivdy   (:,:,:) = 0.0_RP
  check_ddivdy_pl(:,:,:) = 0.0_RP
  check_ddivdz   (:,:,:) = 0.0_RP
  check_ddivdz_pl(:,:,:) = 0.0_RP

  !###############################################################################

  !---< read input data >---
  allocate( ORG_ddivdx           (ADM_gall   ,ADM_kall,ADM_lall   )           )
  allocate( ORG_ddivdx_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl)           )
  allocate( ORG_ddivdy           (ADM_gall   ,ADM_kall,ADM_lall   )           )
  allocate( ORG_ddivdy_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl)           )
  allocate( ORG_ddivdz           (ADM_gall   ,ADM_kall,ADM_lall   )           )
  allocate( ORG_ddivdz_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl)           )
  allocate( ORG_rhogvx           (ADM_gall   ,ADM_kall,ADM_lall   )           )
  allocate( ORG_rhogvx_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl)           )
  allocate( ORG_rhogvy           (ADM_gall   ,ADM_kall,ADM_lall   )           )
  allocate( ORG_rhogvy_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl)           )
  allocate( ORG_rhogvz           (ADM_gall   ,ADM_kall,ADM_lall   )           )
  allocate( ORG_rhogvz_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl)           )
  allocate( ORG_rhogw            (ADM_gall   ,ADM_kall,ADM_lall   )           )
  allocate( ORG_rhogw_pl         (ADM_gall_pl,ADM_kall,ADM_lall_pl)           )
  allocate( ORG_coef_intp        (ADM_gall   ,1:3,ADM_nxyz,TI:TJ,ADM_lall   ) )
  allocate( ORG_coef_intp_pl     (ADM_gall_pl,1:3,ADM_nxyz,      ADM_lall_pl) )
  allocate( ORG_coef_diff        (ADM_gall   ,1:6,ADM_nxyz,      ADM_lall   ) )
  allocate( ORG_coef_diff_pl     (    1:ADM_vlink,ADM_nxyz,      ADM_lall_pl) )
  allocate( ORG_VMTR_RGSQRTH     (ADM_gall   ,ADM_kall,ADM_lall   )           )
  allocate( ORG_VMTR_RGSQRTH_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)           )
  allocate( ORG_VMTR_RGAM        (ADM_gall   ,ADM_kall,ADM_lall   )           )
  allocate( ORG_VMTR_RGAM_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)           )
  allocate( ORG_VMTR_RGAMH       (ADM_gall   ,ADM_kall,ADM_lall   )           )
  allocate( ORG_VMTR_RGAMH_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)           )
  allocate( ORG_VMTR_C2WfactGz   (ADM_gall   ,ADM_kall,6,ADM_lall   )         )
  allocate( ORG_VMTR_C2WfactGz_pl(ADM_gall_pl,ADM_kall,6,ADM_lall_pl)         )

  call dumpio_syscheck
  call dumpio_mk_fname(EX_fname,'snapshot.dyn_divdamp','pe',SET_prc_me-1,6)
  call dumpio_fopen(EX_fid,EX_fname,IO_FREAD)

  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*ADM_lall   ,     ORG_ddivdx           (:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*ADM_lall_pl,     ORG_ddivdx_pl        (:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*ADM_lall   ,     ORG_ddivdy           (:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*ADM_lall_pl,     ORG_ddivdy_pl        (:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*ADM_lall   ,     ORG_ddivdz           (:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*ADM_lall_pl,     ORG_ddivdz_pl        (:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*ADM_lall   ,     ORG_rhogvx           (:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*ADM_lall_pl,     ORG_rhogvx_pl        (:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*ADM_lall   ,     ORG_rhogvy           (:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*ADM_lall_pl,     ORG_rhogvy_pl        (:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*ADM_lall   ,     ORG_rhogvz           (:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*ADM_lall_pl,     ORG_rhogvz_pl        (:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*ADM_lall   ,     ORG_rhogw            (:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*ADM_lall_pl,     ORG_rhogw_pl         (:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_gall   *3*ADM_nxyz*2*ADM_lall   , ORG_coef_intp        (:,:,:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall_pl*3*ADM_nxyz*  ADM_lall_pl, ORG_coef_intp_pl     (:,:,:,:)   )
  call dumpio_read_data( EX_fid, ADM_gall   *6*ADM_nxyz*  ADM_lall   , ORG_coef_diff        (:,:,:,:)   )
  call dumpio_read_data( EX_fid,     ADM_vlink*ADM_nxyz*  ADM_lall_pl, ORG_coef_diff_pl     (:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*ADM_lall   ,     ORG_VMTR_RGSQRTH     (:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*ADM_lall_pl,     ORG_VMTR_RGSQRTH_pl  (:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*ADM_lall   ,     ORG_VMTR_RGAM        (:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*ADM_lall_pl,     ORG_VMTR_RGAM_pl     (:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*ADM_lall   ,     ORG_VMTR_RGAMH       (:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*ADM_lall_pl,     ORG_VMTR_RGAMH_pl    (:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*6*ADM_lall   ,   ORG_VMTR_C2WfactGz   (:,:,:,:)   )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*6*ADM_lall_pl,   ORG_VMTR_C2WfactGz_pl(:,:,:,:)   )

  call dumpio_fclose(EX_fid)

  ddivdx           (:,:,:)     = real( ORG_ddivdx           (:,:,:)    , kind=RP )
  ddivdx_pl        (:,:,:)     = real( ORG_ddivdx_pl        (:,:,:)    , kind=RP )
  ddivdy           (:,:,:)     = real( ORG_ddivdy           (:,:,:)    , kind=RP )
  ddivdy_pl        (:,:,:)     = real( ORG_ddivdy_pl        (:,:,:)    , kind=RP )
  ddivdz           (:,:,:)     = real( ORG_ddivdz           (:,:,:)    , kind=RP )
  ddivdz_pl        (:,:,:)     = real( ORG_ddivdz_pl        (:,:,:)    , kind=RP )
  rhogvx           (:,:,:)     = real( ORG_rhogvx           (:,:,:)    , kind=RP )
  rhogvx_pl        (:,:,:)     = real( ORG_rhogvx_pl        (:,:,:)    , kind=RP )
  rhogvy           (:,:,:)     = real( ORG_rhogvy           (:,:,:)    , kind=RP )
  rhogvy_pl        (:,:,:)     = real( ORG_rhogvy_pl        (:,:,:)    , kind=RP )
  rhogvz           (:,:,:)     = real( ORG_rhogvz           (:,:,:)    , kind=RP )
  rhogvz_pl        (:,:,:)     = real( ORG_rhogvz_pl        (:,:,:)    , kind=RP )
  rhogw            (:,:,:)     = real( ORG_rhogw            (:,:,:)    , kind=RP )
  rhogw_pl         (:,:,:)     = real( ORG_rhogw_pl         (:,:,:)    , kind=RP )
  coef_intp        (:,:,:,:,:) = real( ORG_coef_intp        (:,:,:,:,:), kind=RP )
  coef_intp_pl     (:,:,:,:)   = real( ORG_coef_intp_pl     (:,:,:,:)  , kind=RP )
  coef_diff        (:,:,:,:)   = real( ORG_coef_diff        (:,:,:,:)  , kind=RP )
  coef_diff_pl     (:,:,:)     = real( ORG_coef_diff_pl     (:,:,:)    , kind=RP )
  VMTR_RGSQRTH     (:,:,:)     = real( ORG_VMTR_RGSQRTH     (:,:,:)    , kind=RP )
  VMTR_RGSQRTH_pl  (:,:,:)     = real( ORG_VMTR_RGSQRTH_pl  (:,:,:)    , kind=RP )
  VMTR_RGAM        (:,:,:)     = real( ORG_VMTR_RGAM        (:,:,:)    , kind=RP )
  VMTR_RGAM_pl     (:,:,:)     = real( ORG_VMTR_RGAM_pl     (:,:,:)    , kind=RP )
  VMTR_RGAMH       (:,:,:)     = real( ORG_VMTR_RGAMH       (:,:,:)    , kind=RP )
  VMTR_RGAMH_pl    (:,:,:)     = real( ORG_VMTR_RGAMH_pl    (:,:,:)    , kind=RP )
  VMTR_C2WfactGz   (:,:,:,:)   = real( ORG_VMTR_C2WfactGz   (:,:,:,:)  , kind=RP )
  VMTR_C2WfactGz_pl(:,:,:,:)   = real( ORG_VMTR_C2WfactGz_pl(:,:,:,:)  , kind=RP )

  deallocate( ORG_ddivdx            )
  deallocate( ORG_ddivdx_pl         )
  deallocate( ORG_ddivdy            )
  deallocate( ORG_ddivdy_pl         )
  deallocate( ORG_ddivdz            )
  deallocate( ORG_ddivdz_pl         )
  deallocate( ORG_rhogvx            )
  deallocate( ORG_rhogvx_pl         )
  deallocate( ORG_rhogvy            )
  deallocate( ORG_rhogvy_pl         )
  deallocate( ORG_rhogvz            )
  deallocate( ORG_rhogvz_pl         )
  deallocate( ORG_rhogw             )
  deallocate( ORG_rhogw_pl          )
  deallocate( ORG_coef_intp         )
  deallocate( ORG_coef_intp_pl      )
  deallocate( ORG_coef_diff         )
  deallocate( ORG_coef_diff_pl      )
  deallocate( ORG_VMTR_RGSQRTH      )
  deallocate( ORG_VMTR_RGSQRTH_pl   )
  deallocate( ORG_VMTR_RGAM         )
  deallocate( ORG_VMTR_RGAM_pl      )
  deallocate( ORG_VMTR_RGAMH        )
  deallocate( ORG_VMTR_RGAMH_pl     )
  deallocate( ORG_VMTR_C2WfactGz    )
  deallocate( ORG_VMTR_C2WfactGz_pl )

  !###############################################################################

  call GRD_setup ! allocate GRD_rdgz

  check_ddivdx   (:,:,:) = ddivdx   (:,:,:)
  check_ddivdx_pl(:,:,:) = ddivdx_pl(:,:,:)
  check_ddivdy   (:,:,:) = ddivdy   (:,:,:)
  check_ddivdy_pl(:,:,:) = ddivdy_pl(:,:,:)
  check_ddivdz   (:,:,:) = ddivdz   (:,:,:)
  check_ddivdz_pl(:,:,:) = ddivdz_pl(:,:,:)

  write(*,*) "*** Finish initialize"

  !###############################################################################

  write(*,*) "*** Start  kernel"

  do iteration = 1, SET_iteration

     call DEBUG_rapstart('MAIN_dyn_divdamp')

     call OPRT3D_divdamp( ddivdx   (:,:,:),     ddivdx_pl   (:,:,:),   & ! [OUT]
                          ddivdy   (:,:,:),     ddivdy_pl   (:,:,:),   & ! [OUT]
                          ddivdz   (:,:,:),     ddivdz_pl   (:,:,:),   & ! [OUT]
                          rhogvx   (:,:,:),     rhogvx_pl   (:,:,:),   & ! [IN]
                          rhogvy   (:,:,:),     rhogvy_pl   (:,:,:),   & ! [IN]
                          rhogvz   (:,:,:),     rhogvz_pl   (:,:,:),   & ! [IN]
                          rhogw    (:,:,:),     rhogw_pl    (:,:,:),   & ! [IN]
                          coef_intp(:,:,:,:,:), coef_intp_pl(:,:,:,:), & ! [IN]
                          coef_diff(:,:,:,:),   coef_diff_pl(:,:,:)    ) ! [IN]

     call DEBUG_rapend  ('MAIN_dyn_divdamp')

     write(ADM_LOG_FID,*) '### Input ###'
     call DEBUG_valuecheck( 'check_ddivdx   ', check_ddivdx   (:,:,:) )
     call DEBUG_valuecheck( 'check_ddivdx_pl', check_ddivdx_pl(:,:,:) )
     call DEBUG_valuecheck( 'check_ddivdy   ', check_ddivdy   (:,:,:) )
     call DEBUG_valuecheck( 'check_ddivdy_pl', check_ddivdy_pl(:,:,:) )
     call DEBUG_valuecheck( 'check_ddivdz   ', check_ddivdz   (:,:,:) )
     call DEBUG_valuecheck( 'check_ddivdz_pl', check_ddivdz_pl(:,:,:) )
     call DEBUG_valuecheck( 'rhogvx         ', rhogvx         (:,:,:) )
     call DEBUG_valuecheck( 'rhogvx_pl      ', rhogvx_pl      (:,:,:) )
     call DEBUG_valuecheck( 'rhogvy         ', rhogvy         (:,:,:) )
     call DEBUG_valuecheck( 'rhogvy_pl      ', rhogvy_pl      (:,:,:) )
     call DEBUG_valuecheck( 'rhogvz         ', rhogvz         (:,:,:) )
     call DEBUG_valuecheck( 'rhogvz_pl      ', rhogvz_pl      (:,:,:) )
     call DEBUG_valuecheck( 'rhogw          ', rhogw          (:,:,:) )
     call DEBUG_valuecheck( 'rhogw_pl       ', rhogw_pl       (:,:,:) )

     write(ADM_LOG_FID,*) '### Output ###'
     call DEBUG_valuecheck( 'ddivdx         ', ddivdx         (:,:,:) )
     call DEBUG_valuecheck( 'ddivdx_pl      ', ddivdx_pl      (:,:,:) )
     call DEBUG_valuecheck( 'ddivdy         ', ddivdy         (:,:,:) )
     call DEBUG_valuecheck( 'ddivdy_pl      ', ddivdy_pl      (:,:,:) )
     call DEBUG_valuecheck( 'ddivdz         ', ddivdz         (:,:,:) )
     call DEBUG_valuecheck( 'ddivdz_pl      ', ddivdz_pl      (:,:,:) )
  enddo

  write(ADM_LOG_FID,*) '### Validation : grid-by-grid diff ###'
  check_ddivdx   (:,:,:) = check_ddivdx   (:,:,:) - ddivdx   (:,:,:)
  check_ddivdx_pl(:,:,:) = check_ddivdx_pl(:,:,:) - ddivdx_pl(:,:,:)
  check_ddivdy   (:,:,:) = check_ddivdy   (:,:,:) - ddivdy   (:,:,:)
  check_ddivdy_pl(:,:,:) = check_ddivdy_pl(:,:,:) - ddivdy_pl(:,:,:)
  check_ddivdz   (:,:,:) = check_ddivdz   (:,:,:) - ddivdz   (:,:,:)
  check_ddivdz_pl(:,:,:) = check_ddivdz_pl(:,:,:) - ddivdz_pl(:,:,:)
  call DEBUG_valuecheck( 'check_ddivdx   ', check_ddivdx   (:,:,:) )
  call DEBUG_valuecheck( 'check_ddivdx_pl', check_ddivdx_pl(:,:,:) )
  call DEBUG_valuecheck( 'check_ddivdy   ', check_ddivdy   (:,:,:) )
  call DEBUG_valuecheck( 'check_ddivdy_pl', check_ddivdy_pl(:,:,:) )
  call DEBUG_valuecheck( 'check_ddivdz   ', check_ddivdz   (:,:,:) )
  call DEBUG_valuecheck( 'check_ddivdz_pl', check_ddivdz_pl(:,:,:) )

  write(*,*) "*** Finish kernel"

  !###############################################################################

  call DEBUG_rapreport

  stop
end program dyn_divdamp
!-------------------------------------------------------------------------------
