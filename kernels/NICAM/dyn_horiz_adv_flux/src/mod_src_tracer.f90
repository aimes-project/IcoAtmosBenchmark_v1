!-------------------------------------------------------------------------------
!> Module tracer advection
!!
!! @par Description
!!         This module contains subroutines for tracer advection
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
module mod_src_tracer
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
  !++ Public procedures
  !
  public :: horizontal_flux

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
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Prepare horizontal advection term: mass flux, GRD_xc
  subroutine horizontal_flux( &
       flx_h,  flx_h_pl,  &
       GRD_xc, GRD_xc_pl, &
       rho,    rho_pl,    &
       rhovx,  rhovx_pl,  &
       rhovy,  rhovy_pl,  &
       rhovz,  rhovz_pl,  &
       dt                 )
!ESC!    use mod_const, only: &
!ESC!       CONST_EPS
!ESC!    use mod_adm, only: &
!ESC!       ADM_have_pl,    &
!ESC!       ADM_have_sgp,   &
!ESC!       ADM_lall,       &
!ESC!       ADM_lall_pl,    &
!ESC!       ADM_gall,       &
!ESC!       ADM_gall_pl,    &
!ESC!       ADM_kall,       &
!ESC!       ADM_gall_1d,    &
!ESC!       ADM_gmin,       &
!ESC!       ADM_gmax,       &
!ESC!       ADM_gslf_pl,    &
!ESC!       ADM_gmin_pl,    &
!ESC!       ADM_gmax_pl
!ESC!    use mod_grd, only: &
!ESC!       GRD_xr,   &
!ESC!       GRD_xr_pl
!ESC!    use mod_gmtr, only: &
!ESC!       GMTR_p,    &
!ESC!       GMTR_p_pl, &
!ESC!       GMTR_t,    &
!ESC!       GMTR_t_pl, &
!ESC!       GMTR_a,    &
!ESC!       GMTR_a_pl
    implicit none

    real(RP), intent(out) :: flx_h    (ADM_gall   ,ADM_kall,ADM_lall   ,6)               ! horizontal mass flux
    real(RP), intent(out) :: flx_h_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl  )
    real(RP), intent(out) :: GRD_xc   (ADM_gall   ,ADM_kall,ADM_lall   ,AI:AJ,XDIR:ZDIR) ! mass centroid position
    real(RP), intent(out) :: GRD_xc_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,      XDIR:ZDIR)
    real(RP), intent(in)  :: rho      (ADM_gall   ,ADM_kall,ADM_lall   )                 ! rho at cell center
    real(RP), intent(in)  :: rho_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhovx    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhovx_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhovy    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhovy_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhovz    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhovz_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: dt

    real(RP) :: rhot_TI  (ADM_gall   ) ! rho at cell vertex
    real(RP) :: rhot_TJ  (ADM_gall   ) ! rho at cell vertex
    real(RP) :: rhot_pl  (ADM_gall_pl)
    real(RP) :: rhovxt_TI(ADM_gall   )
    real(RP) :: rhovxt_TJ(ADM_gall   )
    real(RP) :: rhovxt_pl(ADM_gall_pl)
    real(RP) :: rhovyt_TI(ADM_gall   )
    real(RP) :: rhovyt_TJ(ADM_gall   )
    real(RP) :: rhovyt_pl(ADM_gall_pl)
    real(RP) :: rhovzt_TI(ADM_gall   )
    real(RP) :: rhovzt_TJ(ADM_gall   )
    real(RP) :: rhovzt_pl(ADM_gall_pl)

    real(RP) :: rhovxt2
    real(RP) :: rhovyt2
    real(RP) :: rhovzt2
    real(RP) :: flux
    real(RP) :: rrhoa2

    integer  :: gmin, gmax, kall, iall
    real(RP) :: EPS

    integer  :: ij
    integer  :: ip1j, ijp1, ip1jp1
    integer  :: im1j, ijm1

    integer  :: i, j, k, l, n, v
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('____horizontal_adv_flux')

    gmin = ADM_gmin
    gmax = ADM_gmax
    kall = ADM_kall
    iall = ADM_gall_1d

    EPS  = CONST_EPS

    do l = 1, ADM_lall
       !$omp parallel default(none), &
       !$omp private(i,j,k,ij,ip1j,ip1jp1,ijp1,im1j,ijm1,                                        &
       !$omp         rrhoa2,rhovxt2,rhovyt2,rhovzt2,flux),                                       &
       !$omp shared(l,ADM_have_sgp,gmin,gmax,kall,iall,rho,rhovx,rhovy,rhovz,flx_h,dt,           &
       !$omp        rhot_TI,rhovxt_TI,rhovyt_TI,rhovzt_TI,rhot_TJ,rhovxt_TJ,rhovyt_TJ,rhovzt_TJ, &
       !$omp        GRD_xc,GRD_xr,GMTR_p,GMTR_t,GMTR_a,EPS)
       do k = 1, kall

          ! (i,j),(i+1,j)
          !$omp do
          do j = gmin-1, gmax
          do i = gmin-1, gmax
             ij     = (j-1)*iall + i
             ip1j   = ij + 1

             rhot_TI  (ij) = rho  (ij  ,k,l) * GMTR_t(ij,K0,l,TI,W1) &
                           + rho  (ip1j,k,l) * GMTR_t(ij,K0,l,TI,W2)
             rhovxt_TI(ij) = rhovx(ij  ,k,l) * GMTR_t(ij,K0,l,TI,W1) &
                           + rhovx(ip1j,k,l) * GMTR_t(ij,K0,l,TI,W2)
             rhovyt_TI(ij) = rhovy(ij  ,k,l) * GMTR_t(ij,K0,l,TI,W1) &
                           + rhovy(ip1j,k,l) * GMTR_t(ij,K0,l,TI,W2)
             rhovzt_TI(ij) = rhovz(ij  ,k,l) * GMTR_t(ij,K0,l,TI,W1) &
                           + rhovz(ip1j,k,l) * GMTR_t(ij,K0,l,TI,W2)

             rhot_TJ  (ij) = rho  (ij  ,k,l) * GMTR_t(ij,K0,l,TJ,W1)
             rhovxt_TJ(ij) = rhovx(ij  ,k,l) * GMTR_t(ij,K0,l,TJ,W1)
             rhovyt_TJ(ij) = rhovy(ij  ,k,l) * GMTR_t(ij,K0,l,TJ,W1)
             rhovzt_TJ(ij) = rhovz(ij  ,k,l) * GMTR_t(ij,K0,l,TJ,W1)
          enddo
          enddo
          !$omp end do

          ! (i,j+1),(i+1,j+1)
          !$omp do
          do j = gmin-1, gmax
          do i = gmin-1, gmax
             ij     = (j-1)*iall + i
             ijp1   = ij + iall
             ip1jp1 = ij + iall + 1

             rhot_TI  (ij) = rhot_TI  (ij) + rho  (ip1jp1,k,l) * GMTR_t(ij,K0,l,TI,W3)
             rhovxt_TI(ij) = rhovxt_TI(ij) + rhovx(ip1jp1,k,l) * GMTR_t(ij,K0,l,TI,W3)
             rhovyt_TI(ij) = rhovyt_TI(ij) + rhovy(ip1jp1,k,l) * GMTR_t(ij,K0,l,TI,W3)
             rhovzt_TI(ij) = rhovzt_TI(ij) + rhovz(ip1jp1,k,l) * GMTR_t(ij,K0,l,TI,W3)

             rhot_TJ  (ij) = rhot_TJ  (ij) + rho  (ip1jp1,k,l) * GMTR_t(ij,K0,l,TJ,W2) &
                                           + rho  (ijp1  ,k,l) * GMTR_t(ij,K0,l,TJ,W3)
             rhovxt_TJ(ij) = rhovxt_TJ(ij) + rhovx(ip1jp1,k,l) * GMTR_t(ij,K0,l,TJ,W2) &
                                           + rhovx(ijp1  ,k,l) * GMTR_t(ij,K0,l,TJ,W3)
             rhovyt_TJ(ij) = rhovyt_TJ(ij) + rhovy(ip1jp1,k,l) * GMTR_t(ij,K0,l,TJ,W2) &
                                           + rhovy(ijp1  ,k,l) * GMTR_t(ij,K0,l,TJ,W3)
             rhovzt_TJ(ij) = rhovzt_TJ(ij) + rhovz(ip1jp1,k,l) * GMTR_t(ij,K0,l,TJ,W2) &
                                           + rhovz(ijp1  ,k,l) * GMTR_t(ij,K0,l,TJ,W3)
          enddo
          enddo
          !$omp end do

          if ( ADM_have_sgp(l) ) then
             !$omp master
             j = gmin-1
             i = gmin-1

             ij   = (j-1)*iall + i
             ip1j = ij + 1

             rhot_TI  (ij) = rhot_TJ  (ip1j)
             rhovxt_TI(ij) = rhovxt_TJ(ip1j)
             rhovyt_TI(ij) = rhovyt_TJ(ip1j)
             rhovzt_TI(ij) = rhovzt_TJ(ip1j)
             !$omp end master
          endif

          !--- calculate flux and mass centroid position

!OCL XFILL
          !$omp do
          do j = 1, iall
          do i = 1, iall
             if (      i < gmin .OR. i > gmax &
                  .OR. j < gmin .OR. j > gmax ) then
                ij = (j-1)*iall + i

                flx_h(ij,k,l,1) = 0.0_RP
                flx_h(ij,k,l,2) = 0.0_RP
                flx_h(ij,k,l,3) = 0.0_RP
                flx_h(ij,k,l,4) = 0.0_RP
                flx_h(ij,k,l,5) = 0.0_RP
                flx_h(ij,k,l,6) = 0.0_RP

                GRD_xc(ij,k,l,AI ,XDIR) = 0.0_RP
                GRD_xc(ij,k,l,AI ,YDIR) = 0.0_RP
                GRD_xc(ij,k,l,AI ,ZDIR) = 0.0_RP
                GRD_xc(ij,k,l,AIJ,XDIR) = 0.0_RP
                GRD_xc(ij,k,l,AIJ,YDIR) = 0.0_RP
                GRD_xc(ij,k,l,AIJ,ZDIR) = 0.0_RP
                GRD_xc(ij,k,l,AJ ,XDIR) = 0.0_RP
                GRD_xc(ij,k,l,AJ ,YDIR) = 0.0_RP
                GRD_xc(ij,k,l,AJ ,ZDIR) = 0.0_RP
             endif
          enddo
          enddo
          !$omp end do

          !$omp do
          do j = gmin  , gmax
          do i = gmin-1, gmax
             ij     = (j-1)*iall + i
             ip1j   = ij + 1
             ijm1   = ij - iall

             rrhoa2  = 1.0_RP / max( rhot_TJ(ijm1) + rhot_TI(ij), EPS ) ! doubled
             rhovxt2 = rhovxt_TJ(ijm1) + rhovxt_TI(ij)
             rhovyt2 = rhovyt_TJ(ijm1) + rhovyt_TI(ij)
             rhovzt2 = rhovzt_TJ(ijm1) + rhovzt_TI(ij)

             flux = 0.5_RP * ( rhovxt2 * GMTR_a(ij,k0,l,AI ,HNX) &
                             + rhovyt2 * GMTR_a(ij,k0,l,AI ,HNY) &
                             + rhovzt2 * GMTR_a(ij,k0,l,AI ,HNZ) )

             flx_h(ij  ,k,l,1) =  flux * GMTR_p(ij  ,k0,l,P_RAREA) * dt
             flx_h(ip1j,k,l,4) = -flux * GMTR_p(ip1j,k0,l,P_RAREA) * dt

             GRD_xc(ij,k,l,AI,XDIR) = GRD_xr(ij,K0,l,AI,XDIR) - rhovxt2 * rrhoa2 * dt * 0.5_RP
             GRD_xc(ij,k,l,AI,YDIR) = GRD_xr(ij,K0,l,AI,YDIR) - rhovyt2 * rrhoa2 * dt * 0.5_RP
             GRD_xc(ij,k,l,AI,ZDIR) = GRD_xr(ij,K0,l,AI,ZDIR) - rhovzt2 * rrhoa2 * dt * 0.5_RP
          enddo
          enddo
          !$omp end do

          !$omp do
          do j = gmin-1, gmax
          do i = gmin-1, gmax
             ij     = (j-1)*iall + i
             ip1jp1 = ij + iall + 1

             rrhoa2  = 1.0_RP / max( rhot_TI(ij) + rhot_TJ(ij), EPS ) ! doubled
             rhovxt2 = rhovxt_TI(ij) + rhovxt_TJ(ij)
             rhovyt2 = rhovyt_TI(ij) + rhovyt_TJ(ij)
             rhovzt2 = rhovzt_TI(ij) + rhovzt_TJ(ij)

             flux = 0.5_RP * ( rhovxt2 * GMTR_a(ij,k0,l,AIJ,HNX) &
                             + rhovyt2 * GMTR_a(ij,k0,l,AIJ,HNY) &
                             + rhovzt2 * GMTR_a(ij,k0,l,AIJ,HNZ) )

             flx_h(ij    ,k,l,2) =  flux * GMTR_p(ij    ,k0,l,P_RAREA) * dt
             flx_h(ip1jp1,k,l,5) = -flux * GMTR_p(ip1jp1,k0,l,P_RAREA) * dt

             GRD_xc(ij,k,l,AIJ,XDIR) = GRD_xr(ij,K0,l,AIJ,XDIR) - rhovxt2 * rrhoa2 * dt * 0.5_RP
             GRD_xc(ij,k,l,AIJ,YDIR) = GRD_xr(ij,K0,l,AIJ,YDIR) - rhovyt2 * rrhoa2 * dt * 0.5_RP
             GRD_xc(ij,k,l,AIJ,ZDIR) = GRD_xr(ij,K0,l,AIJ,ZDIR) - rhovzt2 * rrhoa2 * dt * 0.5_RP
          enddo
          enddo
          !$omp end do

          !$omp do
          do j = gmin-1, gmax
          do i = gmin  , gmax
             ij     = (j-1)*iall + i
             ijp1   = ij + iall
             im1j   = ij - 1

             rrhoa2  = 1.0_RP / max( rhot_TJ(ij) + rhot_TI(im1j), EPS ) ! doubled
             rhovxt2 = rhovxt_TJ(ij) + rhovxt_TI(im1j)
             rhovyt2 = rhovyt_TJ(ij) + rhovyt_TI(im1j)
             rhovzt2 = rhovzt_TJ(ij) + rhovzt_TI(im1j)

             flux = 0.5_RP * ( rhovxt2 * GMTR_a(ij,k0,l,AJ ,HNX) &
                             + rhovyt2 * GMTR_a(ij,k0,l,AJ ,HNY) &
                             + rhovzt2 * GMTR_a(ij,k0,l,AJ ,HNZ) )

             flx_h(ij  ,k,l,3) =  flux * GMTR_p(ij  ,k0,l,P_RAREA) * dt
             flx_h(ijp1,k,l,6) = -flux * GMTR_p(ijp1,k0,l,P_RAREA) * dt

             GRD_xc(ij,k,l,AJ,XDIR) = GRD_xr(ij,K0,l,AJ,XDIR) - rhovxt2 * rrhoa2 * dt * 0.5_RP
             GRD_xc(ij,k,l,AJ,YDIR) = GRD_xr(ij,K0,l,AJ,YDIR) - rhovyt2 * rrhoa2 * dt * 0.5_RP
             GRD_xc(ij,k,l,AJ,ZDIR) = GRD_xr(ij,K0,l,AJ,ZDIR) - rhovzt2 * rrhoa2 * dt * 0.5_RP
          enddo
          enddo
          !$omp end do

          if ( ADM_have_sgp(l) ) then
             !$omp master
             j = gmin
             i = gmin

             ij = (j-1)*iall + i

             flx_h(ij,k,l,6) = 0.0_RP
             !$omp end master
          endif

       enddo
       !$omp end parallel
    enddo

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             if( ijp1 == ADM_gmax_pl + 1 ) ijp1 = ADM_gmin_pl

             rhot_pl  (v) = rho_pl  (n   ,k,l) * GMTR_t_pl(ij,K0,l,W1) &
                          + rho_pl  (ij  ,k,l) * GMTR_t_pl(ij,K0,l,W2) &
                          + rho_pl  (ijp1,k,l) * GMTR_t_pl(ij,K0,l,W3)
             rhovxt_pl(v) = rhovx_pl(n   ,k,l) * GMTR_t_pl(ij,K0,l,W1) &
                          + rhovx_pl(ij  ,k,l) * GMTR_t_pl(ij,K0,l,W2) &
                          + rhovx_pl(ijp1,k,l) * GMTR_t_pl(ij,K0,l,W3)
             rhovyt_pl(v) = rhovy_pl(n   ,k,l) * GMTR_t_pl(ij,K0,l,W1) &
                          + rhovy_pl(ij  ,k,l) * GMTR_t_pl(ij,K0,l,W2) &
                          + rhovy_pl(ijp1,k,l) * GMTR_t_pl(ij,K0,l,W3)
             rhovzt_pl(v) = rhovz_pl(n   ,k,l) * GMTR_t_pl(ij,K0,l,W1) &
                          + rhovz_pl(ij  ,k,l) * GMTR_t_pl(ij,K0,l,W2) &
                          + rhovz_pl(ijp1,k,l) * GMTR_t_pl(ij,K0,l,W3)
          enddo

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijm1 = v - 1
             if( ijm1 == ADM_gmin_pl - 1 ) ijm1 = ADM_gmax_pl

             rrhoa2  = 1.0_RP / max( rhot_pl(ijm1) + rhot_pl(ij), EPS ) ! doubled
             rhovxt2 = rhovxt_pl(ijm1) + rhovxt_pl(ij)
             rhovyt2 = rhovyt_pl(ijm1) + rhovyt_pl(ij)
             rhovzt2 = rhovzt_pl(ijm1) + rhovzt_pl(ij)

             flux = 0.5_RP * ( rhovxt2 * GMTR_a_pl(ij,K0,l,HNX) &
                             + rhovyt2 * GMTR_a_pl(ij,K0,l,HNY) &
                             + rhovzt2 * GMTR_a_pl(ij,K0,l,HNZ) )

             flx_h_pl(v,k,l) = flux * GMTR_p_pl(n,K0,l,P_RAREA) * dt

             GRD_xc_pl(v,k,l,XDIR) = GRD_xr_pl(v,K0,l,XDIR) - rhovxt2 * rrhoa2 * dt * 0.5_RP
             GRD_xc_pl(v,k,l,YDIR) = GRD_xr_pl(v,K0,l,YDIR) - rhovyt2 * rrhoa2 * dt * 0.5_RP
             GRD_xc_pl(v,k,l,ZDIR) = GRD_xr_pl(v,K0,l,ZDIR) - rhovzt2 * rrhoa2 * dt * 0.5_RP
          enddo

       enddo
       enddo
    endif

    call DEBUG_rapend  ('____horizontal_adv_flux')

    return
  end subroutine horizontal_flux

end module mod_src_tracer
