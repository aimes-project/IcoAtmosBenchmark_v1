!-------------------------------------------------------------------------------
!> Module vertical implicit scheme
!!
!! @par Description
!!          This module is for the caluculation of vertical implicit scheme
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
module mod_vi
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_misc
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: vi_rhow_solver_putcoef
  public :: vi_rhow_solver

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
#ifdef _FIXEDINDEX_
  real(RP), public              :: Mc              (ADM_gall   ,ADM_kall,ADM_lall   )
  real(RP), public              :: Mc_pl           (ADM_gall_pl,ADM_kall,ADM_lall_pl)
  real(RP), public              :: Ml              (ADM_gall   ,ADM_kall,ADM_lall   )
  real(RP), public              :: Ml_pl           (ADM_gall_pl,ADM_kall,ADM_lall_pl)
  real(RP), public              :: Mu              (ADM_gall   ,ADM_kall,ADM_lall   )
  real(RP), public              :: Mu_pl           (ADM_gall_pl,ADM_kall,ADM_lall_pl)
#else
  real(RP), public, allocatable :: Mc              (:,:,:)
  real(RP), public, allocatable :: Mc_pl           (:,:,:)
  real(RP), public, allocatable :: Ml              (:,:,:)
  real(RP), public, allocatable :: Ml_pl           (:,:,:)
  real(RP), public, allocatable :: Mu              (:,:,:)
  real(RP), public, allocatable :: Mu_pl           (:,:,:)
#endif

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine vi_rhow_solver_putcoef( &
       Mc_in, Mc_pl_in, &
       Ml_in, Ml_pl_in, &
       Mu_in, Mu_pl_in  )
    implicit none

    real(RP), intent(in)  :: Mc_in   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: Mc_pl_in(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: Ml_in   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: Ml_pl_in(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: Mu_in   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: Mu_pl_in(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !---------------------------------------------------------------------------

#ifndef _FIXEDINDEX_
    allocate( Mc   (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( Mc_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    allocate( Ml   (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( Ml_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    allocate( Mu   (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( Mu_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
#endif

    Mc   (:,:,:) = Mc_in   (:,:,:)
    Mc_pl(:,:,:) = Mc_pl_in(:,:,:)
    Ml   (:,:,:) = Ml_in   (:,:,:)
    Ml_pl(:,:,:) = Ml_pl_in(:,:,:)
    Mu   (:,:,:) = Mu_in   (:,:,:)
    Mu_pl(:,:,:) = Mu_pl_in(:,:,:)

    return
  end subroutine vi_rhow_solver_putcoef

  !-----------------------------------------------------------------------------
  !> Tridiagonal matrix solver
  subroutine vi_rhow_solver( &
       rhogw,  rhogw_pl,  &
       rhogw0, rhogw0_pl, &
       preg0,  preg0_pl,  &
       rhog0,  rhog0_pl,  &
       Srho,   Srho_pl,   &
       Sw,     Sw_pl,     &
       Spre,   Spre_pl,   &
       dt                 )
!ESC!    use mod_adm, only: &
!ESC!       ADM_have_pl, &
!ESC!       ADM_gall,    &
!ESC!       ADM_gall_pl, &
!ESC!       ADM_lall,    &
!ESC!       ADM_lall_pl, &
!ESC!       ADM_kall,    &
!ESC!       ADM_kmin,    &
!ESC!       ADM_kmax
!ESC!    use mod_const, only: &
!ESC!       CONST_GRAV, &
!ESC!       CONST_Rdry, &
!ESC!       CONST_CVdry
!ESC!    use mod_grd, only: &
!ESC!       GRD_rdgzh, &
!ESC!       GRD_afact, &
!ESC!       GRD_bfact
!ESC!    use mod_vmtr, only: &
!ESC!       VMTR_GSGAM2H,     &
!ESC!       VMTR_GSGAM2H_pl,  &
!ESC!       VMTR_RGAM,        &
!ESC!       VMTR_RGAM_pl,     &
!ESC!       VMTR_RGAMH,       &
!ESC!       VMTR_RGAMH_pl,    &
!ESC!       VMTR_RGSGAM2,     &
!ESC!       VMTR_RGSGAM2_pl,  &
!ESC!       VMTR_RGSGAM2H,    &
!ESC!       VMTR_RGSGAM2H_pl
!ESC!    use mod_runconf, only: &
!ESC!       NON_HYDRO_ALPHA
    !$ use omp_lib
    implicit none

    real(RP), intent(inout) :: rhogw    (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*w          ( G^1/2 x gam2 ), n+1
    real(RP), intent(inout) :: rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP), intent(in)    :: rhogw0   (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*w          ( G^1/2 x gam2 )
    real(RP), intent(in)    :: rhogw0_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: preg0    (ADM_gall   ,ADM_kall,ADM_lall   ) ! pressure prime ( G^1/2 x gam2 )
    real(RP), intent(in)    :: preg0_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: rhog0    (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho            ( G^1/2 x gam2 )
    real(RP), intent(in)    :: rhog0_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: Srho     (ADM_gall   ,ADM_kall,ADM_lall   ) ! source term for rho  at the full level
    real(RP), intent(in)    :: Srho_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: Sw       (ADM_gall   ,ADM_kall,ADM_lall   ) ! source term for rhow at the half level
    real(RP), intent(in)    :: Sw_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: Spre     (ADM_gall   ,ADM_kall,ADM_lall   ) ! source term for pres at the full level
    real(RP), intent(in)    :: Spre_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: dt

    real(RP) :: Sall    (ADM_gall,   ADM_kall)
    real(RP) :: Sall_pl (ADM_gall_pl,ADM_kall)
    real(RP) :: beta    (ADM_gall   )
    real(RP) :: beta_pl (ADM_gall_pl)
    real(RP) :: gamma   (ADM_gall,   ADM_kall)
    real(RP) :: gamma_pl(ADM_gall_pl,ADM_kall)

    integer  :: gall, kmin, kmax, lall
    real(RP) :: grav
    real(RP) :: CVovRt2 ! Cv / R / dt**2
    real(RP) :: alpha

    integer  :: g, k, l
    integer  :: gstr, gend
    !$ integer  :: n_per_thread
    !$ integer  :: n_thread
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('____vi_rhow_solver')

    gall = ADM_gall
    kmin = ADM_kmin
    kmax = ADM_kmax
    lall = ADM_lall

    grav    = CONST_GRAV
    CVovRt2 = CONST_CVdry / CONST_Rdry / (dt*dt)
    alpha   = real(NON_HYDRO_ALPHA,kind=RP)

    !$omp parallel default(none),private(g,k,l), &
    !$omp private(gstr,gend,n_thread,n_per_thread) &
    !$omp shared(gall,kmin,kmax,lall,rhogw,rhogw0,preg0,rhog0,Srho,Sw,Spre,dt,Sall,beta,gamma,Mu,Mc,Ml, &
    !$omp        GRD_afact,GRD_bfact,GRD_rdgzh,VMTR_GSGAM2H,VMTR_RGAM,VMTR_RGAMH,VMTR_RGSGAM2,VMTR_RGSGAM2H,grav,alpha,CVovRt2)
    gstr = 1
    gend = gall
    !$ n_thread     = omp_get_num_threads()
    !$ n_per_thread = gall / n_thread + int( 0.5_RP + sign(0.5_RP,mod(gall,n_thread)-0.5_RP) )
    !$ gstr         = n_per_thread * omp_get_thread_num() + 1
    !$ gend         = min( gstr+n_per_thread-1, gall )

    do l = 1, lall
       ! calc Sall
       do k = kmin+1, kmax
       do g = gstr, gend
          Sall(g,k) = (   ( rhogw0(g,k,  l)*alpha + dt * Sw  (g,k,  l) ) * VMTR_RGAMH  (g,k,  l)**2             &
                      - ( ( preg0 (g,k,  l)       + dt * Spre(g,k,  l) ) * VMTR_RGSGAM2(g,k,  l)                &
                        - ( preg0 (g,k-1,l)       + dt * Spre(g,k-1,l) ) * VMTR_RGSGAM2(g,k-1,l)                &
                        ) * dt * GRD_rdgzh(k)                                                                   &
                      - ( ( rhog0 (g,k,  l)       + dt * Srho(g,k,  l) ) * VMTR_RGAM(g,k,  l)**2 * GRD_afact(k) &
                        + ( rhog0 (g,k-1,l)       + dt * Srho(g,k-1,l) ) * VMTR_RGAM(g,k-1,l)**2 * GRD_bfact(k) &
                        ) * dt * grav                                                                           &
                      ) * CVovRt2
       enddo
       enddo

       ! boundary conditions
       do g = gstr, gend
          rhogw(g,kmin,  l) = rhogw(g,kmin,  l) * VMTR_RGSGAM2H(g,kmin,  l)
          rhogw(g,kmax+1,l) = rhogw(g,kmax+1,l) * VMTR_RGSGAM2H(g,kmax+1,l)
          Sall (g,kmin+1)   = Sall (g,kmin+1) - Ml(g,kmin+1,l) * rhogw(g,kmin,  l)
          Sall (g,kmax  )   = Sall (g,kmax  ) - Mu(g,kmax,  l) * rhogw(g,kmax+1,l)
       enddo

       !---< solve tri-daigonal matrix >

       ! condition at kmin+1
       k = kmin+1
       do g = gstr, gend
          beta (g)     = Mc(g,k,l)
          rhogw(g,k,l) = Sall(g,k) / beta(g)
       enddo

       ! forward
       do k = kmin+2, kmax
       do g = gstr, gend
          gamma(g,k)   = Mu(g,k-1,l) / beta(g)
          beta (g)     = Mc(g,k,l) - Ml(g,k,l) * gamma(g,k) ! update beta
          rhogw(g,k,l) = ( Sall(g,k) - Ml(g,k,l) * rhogw(g,k-1,l) ) / beta(g)
       enddo
       enddo

       ! backward
       do k = kmax-1, kmin+1, -1
       do g = gstr, gend
          rhogw(g,k  ,l) = rhogw(g,k  ,l) - gamma(g,k+1) * rhogw(g,k+1,l)
          rhogw(g,k+1,l) = rhogw(g,k+1,l) * VMTR_GSGAM2H(g,k+1,l) ! return value ( G^1/2 x gam2 )
       enddo
       enddo

       ! boundary treatment
       do g = gstr, gend
          rhogw(g,kmin  ,l) = rhogw(g,kmin  ,l) * VMTR_GSGAM2H(g,kmin  ,l)
          rhogw(g,kmin+1,l) = rhogw(g,kmin+1,l) * VMTR_GSGAM2H(g,kmin+1,l)
          rhogw(g,kmax+1,l) = rhogw(g,kmax+1,l) * VMTR_GSGAM2H(g,kmax+1,l)
       enddo
    enddo
    !$omp end parallel

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
          do k  = ADM_kmin+1, ADM_kmax
          do g = 1, ADM_gall_pl
             Sall_pl(g,k) = (   ( rhogw0_pl(g,k,  l)*alpha + dt * Sw_pl  (g,k,  l) ) * VMTR_RGAMH_pl  (g,k,  l)**2             &
                            - ( ( preg0_pl (g,k,  l)       + dt * Spre_pl(g,k,  l) ) * VMTR_RGSGAM2_pl(g,k,  l)                &
                              - ( preg0_pl (g,k-1,l)       + dt * Spre_pl(g,k-1,l) ) * VMTR_RGSGAM2_pl(g,k-1,l)                &
                              ) * dt * GRD_rdgzh(k)                                                                            &
                            - ( ( rhog0_pl (g,k,  l)       + dt * Srho_pl(g,k,  l) ) * VMTR_RGAM_pl(g,k,  l)**2 * GRD_afact(k) &
                              + ( rhog0_pl (g,k-1,l)       + dt * Srho_pl(g,k-1,l) ) * VMTR_RGAM_pl(g,k-1,l)**2 * GRD_bfact(k) &
                              ) * dt * grav                                                                                    &
                            ) * CVovRt2
          enddo
          enddo

          do g = 1, ADM_gall_pl
             rhogw_pl(g,ADM_kmin,  l) = rhogw_pl(g,ADM_kmin,  l) * VMTR_RGSGAM2H_pl(g,ADM_kmin,  l)
             rhogw_pl(g,ADM_kmax+1,l) = rhogw_pl(g,ADM_kmax+1,l) * VMTR_RGSGAM2H_pl(g,ADM_kmax+1,l)
             Sall_pl (g,ADM_kmin+1)   = Sall_pl (g,ADM_kmin+1) - Ml_pl(g,ADM_kmin+1,l) * rhogw_pl(g,ADM_kmin,  l)
             Sall_pl (g,ADM_kmax  )   = Sall_pl (g,ADM_kmax  ) - Mu_pl(g,ADM_kmax,  l) * rhogw_pl(g,ADM_kmax+1,l)
          enddo

          k = ADM_kmin+1
          do g = 1, ADM_gall_pl
             beta_pl (g)     = Mc_pl(g,k,l)
             rhogw_pl(g,k,l) = Sall_pl(g,k) / beta_pl(g)
          enddo

          do k = ADM_kmin+2, ADM_kmax
          do g = 1, ADM_gall_pl
             gamma_pl(g,k)   = Mu_pl(g,k-1,l) / beta_pl(g)
             beta_pl (g)     = Mc_pl(g,k,l) - Ml_pl(g,k,l) * gamma_pl(g,k) ! update beta
             rhogw_pl(g,k,l) = ( Sall_pl(g,k) - Ml_pl(g,k,l) * rhogw_pl(g,k-1,l) ) / beta_pl(g)
          enddo
          enddo

          ! backward
          do k = ADM_kmax-1, ADM_kmin+1, -1
          do g = 1, ADM_gall_pl
             rhogw_pl(g,k  ,l) = rhogw_pl(g,k  ,l) - gamma_pl(g,k+1) * rhogw_pl(g,k+1,l)
             rhogw_pl(g,k+1,l) = rhogw_pl(g,k+1,l) * VMTR_GSGAM2H_pl(g,k+1,l) ! return value ( G^1/2 x gam2 )
          enddo
          enddo

          ! boundary treatment
          do g = 1, ADM_gall_pl
             rhogw_pl(g,ADM_kmin  ,l) = rhogw_pl(g,ADM_kmin  ,l) * VMTR_GSGAM2H_pl(g,ADM_kmin  ,l)
             rhogw_pl(g,ADM_kmin+1,l) = rhogw_pl(g,ADM_kmin+1,l) * VMTR_GSGAM2H_pl(g,ADM_kmin+1,l)
             rhogw_pl(g,ADM_kmax+1,l) = rhogw_pl(g,ADM_kmax+1,l) * VMTR_GSGAM2H_pl(g,ADM_kmax+1,l)
          enddo
       enddo
    endif

    call DEBUG_rapend('____vi_rhow_solver')

    return
  end subroutine vi_rhow_solver

end module mod_vi
