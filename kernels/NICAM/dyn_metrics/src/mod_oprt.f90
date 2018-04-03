!-------------------------------------------------------------------------------
!> Module operator
!!
!! @par Description
!!          This module contains the subroutines for differential operators
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
module mod_oprt
  !-----------------------------------------------------------------------------
  !
  !++ used modules
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
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  public :: OPRT_divergence_setup
  public :: OPRT_rotation_setup
  public :: OPRT_gradient_setup
  public :: OPRT_laplacian_setup
  public :: OPRT_diffusion_setup

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine OPRT_divergence_setup( &
       GMTR_p,   GMTR_p_pl,  &
       GMTR_t,   GMTR_t_pl,  &
       GMTR_a,   GMTR_a_pl,  &
       coef_div, coef_div_pl )
!ESC!    use mod_adm, only: &
!ESC!       ADM_have_pl,  &
!ESC!       ADM_have_sgp, &
!ESC!       ADM_gall_1d,  &
!ESC!       ADM_gmin,     &
!ESC!       ADM_gmax,     &
!ESC!       ADM_gslf_pl,  &
!ESC!       ADM_gmin_pl,  &
!ESC!       ADM_gmax_pl
!ESC!    use mod_gmtr, only: &
!ESC!       P_RAREA => GMTR_p_RAREA, &
!ESC!       W1      => GMTR_t_W1,    &
!ESC!       W2      => GMTR_t_W2,    &
!ESC!       W3      => GMTR_t_W3,    &
!ESC!       HNX     => GMTR_a_HNX,   &
!ESC!       GMTR_p_nmax,             &
!ESC!       GMTR_t_nmax,             &
!ESC!       GMTR_a_nmax,             &
!ESC!       GMTR_a_nmax_pl
    implicit none

    real(RP), intent(in)  :: GMTR_p     (ADM_gall   ,K0,ADM_lall   ,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_p_pl  (ADM_gall_pl,K0,ADM_lall_pl,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_t     (ADM_gall   ,K0,ADM_lall   ,TI:TJ,GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_t_pl  (ADM_gall_pl,K0,ADM_lall_pl,      GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_a     (ADM_gall   ,K0,ADM_lall   ,AI:AJ,GMTR_a_nmax   )
    real(RP), intent(in)  :: GMTR_a_pl  (ADM_gall_pl,K0,ADM_lall_pl,      GMTR_a_nmax_pl)
    real(RP), intent(out) :: coef_div   (ADM_gall,0:6        ,ADM_nxyz,ADM_lall   )
    real(RP), intent(out) :: coef_div_pl(         0:ADM_vlink,ADM_nxyz,ADM_lall_pl)

    integer  :: gmin, gmax, iall, gall, nxyz, lall

    integer  :: ij
    integer  :: ip1j, ijp1, ip1jp1
    integer  :: im1j, ijm1, im1jm1

    real(RP) :: coef
    integer  :: g, l, d, n, v, hn
    !---------------------------------------------------------------------------

    !if( IO_L ) write(IO_FID_LOG,*) '*** setup coefficient of divergence operator'

    gmin = (ADM_gmin-1)*ADM_gall_1d + ADM_gmin
    gmax = (ADM_gmax-1)*ADM_gall_1d + ADM_gmax
    iall = ADM_gall_1d
    gall = ADM_gall
    nxyz = ADM_nxyz
    lall = ADM_lall

    !$omp parallel workshare
    coef_div   (:,:,:,:) = 0.0_RP
    !$omp end parallel workshare
    coef_div_pl(  :,:,:) = 0.0_RP

    !$omp parallel default(none),private(g,d,l,hn,ij,ip1j,ip1jp1,ijp1,im1j,ijm1,im1jm1), &
    !$omp shared(ADM_have_sgp,gmin,gmax,iall,gall,nxyz,lall,coef_div,GMTR_p,GMTR_t,GMTR_a)
    do l = 1, lall
    do d = 1, nxyz
       hn = d + HNX - 1

       !$omp do
       do g = gmin, gmax
          ij     = g
          ip1j   = g + 1
          ip1jp1 = g + iall + 1
          ijp1   = g + iall
          im1j   = g - 1
          im1jm1 = g - iall - 1
          ijm1   = g - iall

          ! ij
          coef_div(ij,0,d,l) = ( + GMTR_t(ij    ,k0,l,TI,W1) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q1 * b6
                                 + GMTR_t(ij    ,k0,l,TI,W1) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q1 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W1) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q2 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W1) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q2 * b2
                                 + GMTR_t(im1j  ,k0,l,TI,W2) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q3 * b2
                                 - GMTR_t(im1j  ,k0,l,TI,W2) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q3 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W2) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q4 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W2) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q4 * b4
                                 - GMTR_t(im1jm1,k0,l,TI,W3) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q5 * b4
                                 - GMTR_t(im1jm1,k0,l,TI,W3) * GMTR_a(ijm1  ,k0,l,AJ ,hn) & ! Q5 * b5
                                 - GMTR_t(ijm1  ,k0,l,TJ,W3) * GMTR_a(ijm1  ,k0,l,AJ ,hn) & ! Q6 * b5
                                 + GMTR_t(ijm1  ,k0,l,TJ,W3) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q6 * b6
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ip1j
          coef_div(ij,1,d,l) = ( - GMTR_t(ijm1  ,k0,l,TJ,W2) * GMTR_a(ijm1  ,k0,l,AJ ,hn) & ! Q6 * b5
                                 + GMTR_t(ijm1  ,k0,l,TJ,W2) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q6 * b6
                                 + GMTR_t(ij    ,k0,l,TI,W2) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q1 * b6
                                 + GMTR_t(ij    ,k0,l,TI,W2) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q1 * b1
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ip1jp1
          coef_div(ij,2,d,l) = ( + GMTR_t(ij    ,k0,l,TI,W3) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q1 * b6
                                 + GMTR_t(ij    ,k0,l,TI,W3) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q1 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W2) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q2 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W2) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q2 * b2
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ijp1
          coef_div(ij,3,d,l) = ( + GMTR_t(ij    ,k0,l,TJ,W3) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q2 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W3) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q2 * b2
                                 + GMTR_t(im1j  ,k0,l,TI,W3) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q3 * b2
                                 - GMTR_t(im1j  ,k0,l,TI,W3) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q3 * b3
                               ) * 0.5_RP*GMTR_p(ij,k0,l,P_RAREA)
          ! im1j
          coef_div(ij,4,d,l) = ( + GMTR_t(im1j  ,k0,l,TI,W1) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q3 * b2
                                 - GMTR_t(im1j  ,k0,l,TI,W1) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q3 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W3) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q4 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W3) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q4 * b4
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! im1jm1
          coef_div(ij,5,d,l) = ( - GMTR_t(im1jm1,k0,l,TJ,W1) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q4 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W1) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q4 * b4
                                 - GMTR_t(im1jm1,k0,l,TI,W1) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q5 * b4
                                 - GMTR_t(im1jm1,k0,l,TI,W1) * GMTR_a(ijm1  ,k0,l,AJ ,hn) & ! Q5 * b5
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ijm1
          coef_div(ij,6,d,l) = ( - GMTR_t(im1jm1,k0,l,TI,W2) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q5 * b4
                                 - GMTR_t(im1jm1,k0,l,TI,W2) * GMTR_a(ijm1  ,k0,l,AJ ,hn) & ! Q5 * b5
                                 - GMTR_t(ijm1  ,k0,l,TJ,W1) * GMTR_a(ijm1  ,k0,l,AJ ,hn) & ! Q6 * b5
                                 + GMTR_t(ijm1  ,k0,l,TJ,W1) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q6 * b6
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
       enddo
       !$omp end do

       if ( ADM_have_sgp(l) ) then ! pentagon
          !$omp master
          ij     = gmin
          ip1j   = gmin + 1
          ip1jp1 = gmin + iall + 1
          ijp1   = gmin + iall
          im1j   = gmin - 1
          im1jm1 = gmin - iall - 1
          ijm1   = gmin - iall

          ! ij
          coef_div(ij,0,d,l) = ( + GMTR_t(ij    ,k0,l,TI,W1) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q1 * b6
                                 + GMTR_t(ij    ,k0,l,TI,W1) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q1 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W1) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q2 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W1) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q2 * b2
                                 + GMTR_t(im1j  ,k0,l,TI,W2) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q3 * b2
                                 - GMTR_t(im1j  ,k0,l,TI,W2) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q3 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W2) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q4 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W2) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q4 * b4
                                 - GMTR_t(ijm1  ,k0,l,TJ,W3) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q6 * b4
                                 + GMTR_t(ijm1  ,k0,l,TJ,W3) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q6 * b6
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ip1j
          coef_div(ij,1,d,l) = ( - GMTR_t(ijm1  ,k0,l,TJ,W2) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q6 * b4
                                 + GMTR_t(ijm1  ,k0,l,TJ,W2) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q6 * b6
                                 + GMTR_t(ij    ,k0,l,TI,W2) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q1 * b6
                                 + GMTR_t(ij    ,k0,l,TI,W2) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q1 * b1
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ip1jp1
          coef_div(ij,2,d,l) = ( + GMTR_t(ij    ,k0,l,TI,W3) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q1 * b6
                                 + GMTR_t(ij    ,k0,l,TI,W3) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q1 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W2) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q2 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W2) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q2 * b2
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ijp1
          coef_div(ij,3,d,l) = ( + GMTR_t(ij    ,k0,l,TJ,W3) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q2 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W3) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q2 * b2
                                 + GMTR_t(im1j  ,k0,l,TI,W3) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q3 * b2
                                 - GMTR_t(im1j  ,k0,l,TI,W3) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q3 * b3
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! im1j
          coef_div(ij,4,d,l) = ( + GMTR_t(im1j  ,k0,l,TI,W1) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q3 * b2
                                 - GMTR_t(im1j  ,k0,l,TI,W1) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q3 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W3) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q4 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W3) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q4 * b4
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! im1jm1
          coef_div(ij,5,d,l) = ( - GMTR_t(im1jm1,k0,l,TJ,W1) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q4 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W1) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q4 * b4
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ijm1
          coef_div(ij,6,d,l) = ( - GMTR_t(ijm1  ,k0,l,TJ,W1) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q6 * b4
                                 + GMTR_t(ijm1  ,k0,l,TJ,W1) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q6 * b6
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          !$omp end master
       endif

    enddo ! loop d
    enddo ! loop l
    !$omp end parallel

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl
       do l = 1, ADM_lall_pl
       do d = 1, ADM_nxyz
          hn = d + HNX - 1

          coef = 0.0_RP
          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             if( ijp1 == ADM_gmax_pl+1 ) ijp1 = ADM_gmin_pl

             coef = coef + ( GMTR_t_pl(ij,k0,l,W1) * GMTR_a_pl(ij  ,k0,l,hn) &
                           + GMTR_t_pl(ij,k0,l,W1) * GMTR_a_pl(ijp1,k0,l,hn) )
          enddo
          coef_div_pl(0,d,l) = coef * 0.5_RP * GMTR_p_pl(n,k0,l,P_RAREA)

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             ijm1 = v - 1
             if( ijp1 == ADM_gmax_pl + 1 ) ijp1 = ADM_gmin_pl
             if( ijm1 == ADM_gmin_pl - 1 ) ijm1 = ADM_gmax_pl

             coef_div_pl(v-1,d,l) = ( + GMTR_t_pl(ijm1,k0,l,W3) * GMTR_a_pl(ijm1,k0,l,hn) &
                                      + GMTR_t_pl(ijm1,k0,l,W3) * GMTR_a_pl(ij  ,k0,l,hn) &
                                      + GMTR_t_pl(ij  ,k0,l,W2) * GMTR_a_pl(ij  ,k0,l,hn) &
                                      + GMTR_t_pl(ij  ,k0,l,W2) * GMTR_a_pl(ijp1,k0,l,hn) &
                                    ) * 0.5_RP * GMTR_p_pl(n,k0,l,P_RAREA)
          enddo
       enddo ! loop d
       enddo ! loop l
    endif

    return
  end subroutine OPRT_divergence_setup

  !-----------------------------------------------------------------------------
  subroutine OPRT_rotation_setup( &
       GMTR_p,   GMTR_p_pl,  &
       GMTR_t,   GMTR_t_pl,  &
       GMTR_a,   GMTR_a_pl,  &
       coef_rot, coef_rot_pl )
!ESC!    use mod_adm, only: &
!ESC!       ADM_have_pl,  &
!ESC!       ADM_have_sgp, &
!ESC!       ADM_gall_1d,  &
!ESC!       ADM_gmin,     &
!ESC!       ADM_gmax,     &
!ESC!       ADM_gslf_pl,  &
!ESC!       ADM_gmin_pl,  &
!ESC!       ADM_gmax_pl
!ESC!    use mod_gmtr, only: &
!ESC!       P_RAREA => GMTR_p_RAREA, &
!ESC!       W1      => GMTR_t_W1,    &
!ESC!       W2      => GMTR_t_W2,    &
!ESC!       W3      => GMTR_t_W3,    &
!ESC!       HTX     => GMTR_a_HTX,   &
!ESC!       GMTR_p_nmax,             &
!ESC!       GMTR_t_nmax,             &
!ESC!       GMTR_a_nmax,             &
!ESC!       GMTR_a_nmax_pl
    implicit none

    real(RP), intent(in)  :: GMTR_p     (ADM_gall   ,K0,ADM_lall   ,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_p_pl  (ADM_gall_pl,K0,ADM_lall_pl,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_t     (ADM_gall   ,K0,ADM_lall   ,TI:TJ,GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_t_pl  (ADM_gall_pl,K0,ADM_lall_pl,      GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_a     (ADM_gall   ,K0,ADM_lall   ,AI:AJ,GMTR_a_nmax   )
    real(RP), intent(in)  :: GMTR_a_pl  (ADM_gall_pl,K0,ADM_lall_pl,      GMTR_a_nmax_pl)
    real(RP), intent(out) :: coef_rot   (ADM_gall,0:6        ,ADM_nxyz,ADM_lall   )
    real(RP), intent(out) :: coef_rot_pl(         0:ADM_vlink,ADM_nxyz,ADM_lall_pl)

    integer  :: gmin, gmax, iall, gall, nxyz, lall

    integer  :: ij
    integer  :: ip1j, ijp1, ip1jp1
    integer  :: im1j, ijm1, im1jm1

    real(RP) :: coef
    integer  :: g, l, d, n, v, ht
    !---------------------------------------------------------------------------

    !if( IO_L ) write(IO_FID_LOG,*) '*** setup coefficient of rotation operator'

    gmin = (ADM_gmin-1)*ADM_gall_1d + ADM_gmin
    gmax = (ADM_gmax-1)*ADM_gall_1d + ADM_gmax
    iall = ADM_gall_1d
    gall = ADM_gall
    nxyz = ADM_nxyz
    lall = ADM_lall

    !$omp parallel workshare
    coef_rot   (:,:,:,:) = 0.0_RP
    !$omp end parallel workshare
    coef_rot_pl(  :,:,:) = 0.0_RP

    !$omp parallel default(none),private(g,d,l,ht,ij,ip1j,ip1jp1,ijp1,im1j,ijm1,im1jm1), &
    !$omp shared(ADM_have_sgp,gmin,gmax,iall,gall,nxyz,lall,coef_rot,GMTR_p,GMTR_t,GMTR_a)
    do l = 1, lall
    do d = 1, nxyz
       ht = d + HTX - 1

       !$omp do
       do g = gmin, gmax
          ij     = g
          ip1j   = g + 1
          ip1jp1 = g + iall + 1
          ijp1   = g + iall
          im1j   = g - 1
          im1jm1 = g - iall - 1
          ijm1   = g - iall

          ! ij
          coef_rot(ij,0,d,l) = ( + GMTR_t(ij    ,k0,l,TI,W1) * GMTR_a(ij    ,k0,l,AI ,ht) & ! Q1 * b6
                                 + GMTR_t(ij    ,k0,l,TI,W1) * GMTR_a(ij    ,k0,l,AIJ,ht) & ! Q1 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W1) * GMTR_a(ij    ,k0,l,AIJ,ht) & ! Q2 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W1) * GMTR_a(ij    ,k0,l,AJ ,ht) & ! Q2 * b2
                                 + GMTR_t(im1j  ,k0,l,TI,W2) * GMTR_a(ij    ,k0,l,AJ ,ht) & ! Q3 * b2
                                 - GMTR_t(im1j  ,k0,l,TI,W2) * GMTR_a(im1j  ,k0,l,AI ,ht) & ! Q3 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W2) * GMTR_a(im1j  ,k0,l,AI ,ht) & ! Q4 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W2) * GMTR_a(im1jm1,k0,l,AIJ,ht) & ! Q4 * b4
                                 - GMTR_t(im1jm1,k0,l,TI,W3) * GMTR_a(im1jm1,k0,l,AIJ,ht) & ! Q5 * b4
                                 - GMTR_t(im1jm1,k0,l,TI,W3) * GMTR_a(ijm1  ,k0,l,AJ ,ht) & ! Q5 * b5
                                 - GMTR_t(ijm1  ,k0,l,TJ,W3) * GMTR_a(ijm1  ,k0,l,AJ ,ht) & ! Q6 * b5
                                 + GMTR_t(ijm1  ,k0,l,TJ,W3) * GMTR_a(ij    ,k0,l,AI ,ht) & ! Q6 * b6
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ip1j
          coef_rot(ij,1,d,l) = ( - GMTR_t(ijm1  ,k0,l,TJ,W2) * GMTR_a(ijm1  ,k0,l,AJ ,ht) & ! Q6 * b5
                                 + GMTR_t(ijm1  ,k0,l,TJ,W2) * GMTR_a(ij    ,k0,l,AI ,ht) & ! Q6 * b6
                                 + GMTR_t(ij    ,k0,l,TI,W2) * GMTR_a(ij    ,k0,l,AI ,ht) & ! Q1 * b6
                                 + GMTR_t(ij    ,k0,l,TI,W2) * GMTR_a(ij    ,k0,l,AIJ,ht) & ! Q1 * b1
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ip1jp1
          coef_rot(ij,2,d,l) = ( + GMTR_t(ij    ,k0,l,TI,W3) * GMTR_a(ij    ,k0,l,AI ,ht) & ! Q1 * b6
                                 + GMTR_t(ij    ,k0,l,TI,W3) * GMTR_a(ij    ,k0,l,AIJ,ht) & ! Q1 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W2) * GMTR_a(ij    ,k0,l,AIJ,ht) & ! Q2 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W2) * GMTR_a(ij    ,k0,l,AJ ,ht) & ! Q2 * b2
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ijp1
          coef_rot(ij,3,d,l) = ( + GMTR_t(ij    ,k0,l,TJ,W3) * GMTR_a(ij    ,k0,l,AIJ,ht) & ! Q2 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W3) * GMTR_a(ij    ,k0,l,AJ ,ht) & ! Q2 * b2
                                 + GMTR_t(im1j  ,k0,l,TI,W3) * GMTR_a(ij    ,k0,l,AJ ,ht) & ! Q3 * b2
                                 - GMTR_t(im1j  ,k0,l,TI,W3) * GMTR_a(im1j  ,k0,l,AI ,ht) & ! Q3 * b3
                               ) * 0.5_RP*GMTR_p(ij,k0,l,P_RAREA)
          ! im1j
          coef_rot(ij,4,d,l) = ( + GMTR_t(im1j  ,k0,l,TI,W1) * GMTR_a(ij    ,k0,l,AJ ,ht) & ! Q3 * b2
                                 - GMTR_t(im1j  ,k0,l,TI,W1) * GMTR_a(im1j  ,k0,l,AI ,ht) & ! Q3 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W3) * GMTR_a(im1j  ,k0,l,AI ,ht) & ! Q4 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W3) * GMTR_a(im1jm1,k0,l,AIJ,ht) & ! Q4 * b4
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! im1jm1
          coef_rot(ij,5,d,l) = ( - GMTR_t(im1jm1,k0,l,TJ,W1) * GMTR_a(im1j  ,k0,l,AI ,ht) & ! Q4 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W1) * GMTR_a(im1jm1,k0,l,AIJ,ht) & ! Q4 * b4
                                 - GMTR_t(im1jm1,k0,l,TI,W1) * GMTR_a(im1jm1,k0,l,AIJ,ht) & ! Q5 * b4
                                 - GMTR_t(im1jm1,k0,l,TI,W1) * GMTR_a(ijm1  ,k0,l,AJ ,ht) & ! Q5 * b5
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ijm1
          coef_rot(ij,6,d,l) = ( - GMTR_t(im1jm1,k0,l,TI,W2) * GMTR_a(im1jm1,k0,l,AIJ,ht) & ! Q5 * b4
                                 - GMTR_t(im1jm1,k0,l,TI,W2) * GMTR_a(ijm1  ,k0,l,AJ ,ht) & ! Q5 * b5
                                 - GMTR_t(ijm1  ,k0,l,TJ,W1) * GMTR_a(ijm1  ,k0,l,AJ ,ht) & ! Q6 * b5
                                 + GMTR_t(ijm1  ,k0,l,TJ,W1) * GMTR_a(ij    ,k0,l,AI ,ht) & ! Q6 * b6
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
       enddo
       !$omp end do

       if ( ADM_have_sgp(l) ) then ! pentagon
          !$omp master
          ij     = gmin
          ip1j   = gmin + 1
          ip1jp1 = gmin + iall + 1
          ijp1   = gmin + iall
          im1j   = gmin - 1
          im1jm1 = gmin - iall - 1
          ijm1   = gmin - iall

          ! ij
          coef_rot(ij,0,d,l) = ( + GMTR_t(ij    ,k0,l,TI,W1) * GMTR_a(ij    ,k0,l,AI ,ht) & ! Q1 * b6
                                 + GMTR_t(ij    ,k0,l,TI,W1) * GMTR_a(ij    ,k0,l,AIJ,ht) & ! Q1 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W1) * GMTR_a(ij    ,k0,l,AIJ,ht) & ! Q2 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W1) * GMTR_a(ij    ,k0,l,AJ ,ht) & ! Q2 * b2
                                 + GMTR_t(im1j  ,k0,l,TI,W2) * GMTR_a(ij    ,k0,l,AJ ,ht) & ! Q3 * b2
                                 - GMTR_t(im1j  ,k0,l,TI,W2) * GMTR_a(im1j  ,k0,l,AI ,ht) & ! Q3 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W2) * GMTR_a(im1j  ,k0,l,AI ,ht) & ! Q4 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W2) * GMTR_a(im1jm1,k0,l,AIJ,ht) & ! Q4 * b4
                                 - GMTR_t(ijm1  ,k0,l,TJ,W3) * GMTR_a(im1jm1,k0,l,AIJ,ht) & ! Q6 * b4
                                 + GMTR_t(ijm1  ,k0,l,TJ,W3) * GMTR_a(ij    ,k0,l,AI ,ht) & ! Q6 * b6
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ip1j
          coef_rot(ij,1,d,l) = ( - GMTR_t(ijm1  ,k0,l,TJ,W2) * GMTR_a(im1jm1,k0,l,AIJ,ht) & ! Q6 * b4
                                 + GMTR_t(ijm1  ,k0,l,TJ,W2) * GMTR_a(ij    ,k0,l,AI ,ht) & ! Q6 * b6
                                 + GMTR_t(ij    ,k0,l,TI,W2) * GMTR_a(ij    ,k0,l,AI ,ht) & ! Q1 * b6
                                 + GMTR_t(ij    ,k0,l,TI,W2) * GMTR_a(ij    ,k0,l,AIJ,ht) & ! Q1 * b1
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ip1jp1
          coef_rot(ij,2,d,l) = ( + GMTR_t(ij    ,k0,l,TI,W3) * GMTR_a(ij    ,k0,l,AI ,ht) & ! Q1 * b6
                                 + GMTR_t(ij    ,k0,l,TI,W3) * GMTR_a(ij    ,k0,l,AIJ,ht) & ! Q1 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W2) * GMTR_a(ij    ,k0,l,AIJ,ht) & ! Q2 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W2) * GMTR_a(ij    ,k0,l,AJ ,ht) & ! Q2 * b2
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ijp1
          coef_rot(ij,3,d,l) = ( + GMTR_t(ij    ,k0,l,TJ,W3) * GMTR_a(ij    ,k0,l,AIJ,ht) & ! Q2 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W3) * GMTR_a(ij    ,k0,l,AJ ,ht) & ! Q2 * b2
                                 + GMTR_t(im1j  ,k0,l,TI,W3) * GMTR_a(ij    ,k0,l,AJ ,ht) & ! Q3 * b2
                                 - GMTR_t(im1j  ,k0,l,TI,W3) * GMTR_a(im1j  ,k0,l,AI ,ht) & ! Q3 * b3
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! im1j
          coef_rot(ij,4,d,l) = ( + GMTR_t(im1j  ,k0,l,TI,W1) * GMTR_a(ij    ,k0,l,AJ ,ht) & ! Q3 * b2
                                 - GMTR_t(im1j  ,k0,l,TI,W1) * GMTR_a(im1j  ,k0,l,AI ,ht) & ! Q3 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W3) * GMTR_a(im1j  ,k0,l,AI ,ht) & ! Q4 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W3) * GMTR_a(im1jm1,k0,l,AIJ,ht) & ! Q4 * b4
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! im1jm1
          coef_rot(ij,5,d,l) = ( - GMTR_t(im1jm1,k0,l,TJ,W1) * GMTR_a(im1j  ,k0,l,AI ,ht) & ! Q4 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W1) * GMTR_a(im1jm1,k0,l,AIJ,ht) & ! Q4 * b4
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ijm1
          coef_rot(ij,6,d,l) = ( - GMTR_t(ijm1  ,k0,l,TJ,W1) * GMTR_a(im1jm1,k0,l,AIJ,ht) & ! Q6 * b4
                                 + GMTR_t(ijm1  ,k0,l,TJ,W1) * GMTR_a(ij    ,k0,l,AI ,ht) & ! Q6 * b6
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          !$omp end master
       endif

    enddo ! loop d
    enddo ! loop l
    !$omp end parallel

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl
       do l = 1, ADM_lall_pl
       do d = 1, ADM_nxyz
          ht = d + HTX - 1

          coef = 0.0_RP
          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             if( ijp1 == ADM_gmax_pl+1 ) ijp1 = ADM_gmin_pl

             coef = coef + ( GMTR_t_pl(ij,k0,l,W1) * GMTR_a_pl(ij  ,k0,l,ht) &
                           + GMTR_t_pl(ij,k0,l,W1) * GMTR_a_pl(ijp1,k0,l,ht) )
          enddo
          coef_rot_pl(0,d,l) = coef * 0.5_RP * GMTR_p_pl(n,k0,l,P_RAREA)

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             ijm1 = v - 1
             if( ijp1 == ADM_gmax_pl + 1 ) ijp1 = ADM_gmin_pl
             if( ijm1 == ADM_gmin_pl - 1 ) ijm1 = ADM_gmax_pl

             coef_rot_pl(v-1,d,l) = ( + GMTR_t_pl(ijm1,k0,l,W3) * GMTR_a_pl(ijm1,k0,l,ht) &
                                      + GMTR_t_pl(ijm1,k0,l,W3) * GMTR_a_pl(ij  ,k0,l,ht) &
                                      + GMTR_t_pl(ij  ,k0,l,W2) * GMTR_a_pl(ij  ,k0,l,ht) &
                                      + GMTR_t_pl(ij  ,k0,l,W2) * GMTR_a_pl(ijp1,k0,l,ht) &
                                    ) * 0.5_RP * GMTR_p_pl(n,k0,l,P_RAREA)
          enddo
       enddo ! loop d
       enddo ! loop l
    endif

    return
  end subroutine OPRT_rotation_setup

  !-----------------------------------------------------------------------------
  subroutine OPRT_gradient_setup( &
       GMTR_p,    GMTR_p_pl,   &
       GMTR_t,    GMTR_t_pl,   &
       GMTR_a,    GMTR_a_pl,   &
       coef_grad, coef_grad_pl )
!ESC!    use mod_adm, only: &
!ESC!       ADM_have_pl,  &
!ESC!       ADM_have_sgp, &
!ESC!       ADM_gall_1d,  &
!ESC!       ADM_gmin,     &
!ESC!       ADM_gmax,     &
!ESC!       ADM_gslf_pl,  &
!ESC!       ADM_gmin_pl,  &
!ESC!       ADM_gmax_pl
!ESC!    use mod_gmtr, only: &
!ESC!       P_RAREA => GMTR_p_RAREA, &
!ESC!       W1      => GMTR_t_W1,    &
!ESC!       W2      => GMTR_t_W2,    &
!ESC!       W3      => GMTR_t_W3,    &
!ESC!       HNX     => GMTR_a_HNX,   &
!ESC!       GMTR_p_nmax,             &
!ESC!       GMTR_t_nmax,             &
!ESC!       GMTR_a_nmax,             &
!ESC!       GMTR_a_nmax_pl
    implicit none

    real(RP), intent(in)  :: GMTR_p      (ADM_gall   ,K0,ADM_lall   ,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_p_pl   (ADM_gall_pl,K0,ADM_lall_pl,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_t      (ADM_gall   ,K0,ADM_lall   ,TI:TJ,GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_t_pl   (ADM_gall_pl,K0,ADM_lall_pl,      GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_a      (ADM_gall   ,K0,ADM_lall   ,AI:AJ,GMTR_a_nmax   )
    real(RP), intent(in)  :: GMTR_a_pl   (ADM_gall_pl,K0,ADM_lall_pl,      GMTR_a_nmax_pl)
    real(RP), intent(out) :: coef_grad   (ADM_gall,0:6        ,ADM_nxyz,ADM_lall   )
    real(RP), intent(out) :: coef_grad_pl(         0:ADM_vlink,ADM_nxyz,ADM_lall_pl)

    integer  :: gmin, gmax, iall, gall, nxyz, lall

    integer  :: ij
    integer  :: ip1j, ijp1, ip1jp1
    integer  :: im1j, ijm1, im1jm1

    real(RP) :: coef
    integer  :: g, l, d, n, v, hn
    !---------------------------------------------------------------------------

    !if( IO_L ) write(IO_FID_LOG,*) '*** setup coefficient of gradient operator'

    gmin = (ADM_gmin-1)*ADM_gall_1d + ADM_gmin
    gmax = (ADM_gmax-1)*ADM_gall_1d + ADM_gmax
    iall = ADM_gall_1d
    gall = ADM_gall
    nxyz = ADM_nxyz
    lall = ADM_lall

    !$omp parallel workshare
    coef_grad   (:,:,:,:) = 0.0_RP
    !$omp end parallel workshare
    coef_grad_pl(  :,:,:) = 0.0_RP

    !$omp parallel default(none),private(g,d,l,hn,ij,ip1j,ip1jp1,ijp1,im1j,ijm1,im1jm1), &
    !$omp shared(ADM_have_sgp,gmin,gmax,iall,gall,nxyz,lall,coef_grad,GMTR_p,GMTR_t,GMTR_a)
    do l = 1, lall
    do d = 1, nxyz
       hn = d + HNX - 1

       !$omp do
       do g = gmin, gmax
          ij     = g
          ip1j   = g + 1
          ip1jp1 = g + iall + 1
          ijp1   = g + iall
          im1j   = g - 1
          im1jm1 = g - iall - 1
          ijm1   = g - iall

          ! ij
          coef_grad(ij,0,d,l) = ( + GMTR_t(ij    ,k0,l,TI,W1) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q1 * b6
                                  + GMTR_t(ij    ,k0,l,TI,W1) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q1 * b1
                                  + GMTR_t(ij    ,k0,l,TJ,W1) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q2 * b1
                                  + GMTR_t(ij    ,k0,l,TJ,W1) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q2 * b2
                                  + GMTR_t(im1j  ,k0,l,TI,W2) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q3 * b2
                                  - GMTR_t(im1j  ,k0,l,TI,W2) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q3 * b3
                                  - GMTR_t(im1jm1,k0,l,TJ,W2) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q4 * b3
                                  - GMTR_t(im1jm1,k0,l,TI,W3) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q4 * b4
                                  - GMTR_t(im1jm1,k0,l,TJ,W2) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q5 * b4
                                  - GMTR_t(im1jm1,k0,l,TI,W3) * GMTR_a(ijm1  ,k0,l,AJ ,hn) & ! Q5 * b5
                                  - GMTR_t(ijm1  ,k0,l,TJ,W3) * GMTR_a(ijm1  ,k0,l,AJ ,hn) & ! Q6 * b5
                                  + GMTR_t(ijm1  ,k0,l,TJ,W3) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q6 * b6
                                  - 2.0_RP * GMTR_a(ij    ,k0,l,AIJ,hn)                    & ! P0 * b1
                                  - 2.0_RP * GMTR_a(ij    ,k0,l,AJ ,hn)                    & ! P0 * b2
                                  + 2.0_RP * GMTR_a(im1j  ,k0,l,AI ,hn)                    & ! P0 * b3
                                  + 2.0_RP * GMTR_a(im1jm1,k0,l,AIJ,hn)                    & ! P0 * b4
                                  + 2.0_RP * GMTR_a(ijm1  ,k0,l,AJ ,hn)                    & ! P0 * b5
                                  - 2.0_RP * GMTR_a(ij    ,k0,l,AI ,hn)                    & ! P0 * b6
                                ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ip1j
          coef_grad(ij,1,d,l) = ( - GMTR_t(ijm1  ,k0,l,TJ,W2) * GMTR_a(ijm1  ,k0,l,AJ ,hn) & ! Q6 * b5
                                  + GMTR_t(ijm1  ,k0,l,TJ,W2) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q6 * b6
                                  + GMTR_t(ij    ,k0,l,TI,W2) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q1 * b6
                                  + GMTR_t(ij    ,k0,l,TI,W2) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q1 * b1
                                ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ip1jp1
          coef_grad(ij,2,d,l) = ( + GMTR_t(ij    ,k0,l,TI,W3) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q1 * b6
                                  + GMTR_t(ij    ,k0,l,TI,W3) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q1 * b1
                                  + GMTR_t(ij    ,k0,l,TJ,W2) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q2 * b1
                                  + GMTR_t(ij    ,k0,l,TJ,W2) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q2 * b2
                                ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ijp1
          coef_grad(ij,3,d,l) = ( + GMTR_t(ij    ,k0,l,TJ,W3) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q2 * b1
                                  + GMTR_t(ij    ,k0,l,TJ,W3) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q2 * b2
                                  + GMTR_t(im1j  ,k0,l,TI,W3) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q3 * b2
                                  - GMTR_t(im1j  ,k0,l,TI,W3) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q3 * b3
                                ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! im1j
          coef_grad(ij,4,d,l) = ( + GMTR_t(im1j  ,k0,l,TI,W1) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q3 * b2
                                  - GMTR_t(im1j  ,k0,l,TI,W1) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q3 * b3
                                  - GMTR_t(im1jm1,k0,l,TJ,W3) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q4 * b3
                                  - GMTR_t(im1jm1,k0,l,TJ,W3) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q4 * b4
                                ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! im1jm1
          coef_grad(ij,5,d,l) = ( - GMTR_t(im1jm1,k0,l,TJ,W1) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q4 * b3
                                  - GMTR_t(im1jm1,k0,l,TJ,W1) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q4 * b4
                                  - GMTR_t(im1jm1,k0,l,TI,W1) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q5 * b4
                                  - GMTR_t(im1jm1,k0,l,TI,W1) * GMTR_a(ijm1  ,k0,l,AJ ,hn) & ! Q5 * b5
                                ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ijm1
          coef_grad(ij,6,d,l) = ( - GMTR_t(im1jm1,k0,l,TI,W2) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q5 * b4
                                  - GMTR_t(im1jm1,k0,l,TI,W2) * GMTR_a(ijm1  ,k0,l,AJ ,hn) & ! Q5 * b5
                                  - GMTR_t(ijm1  ,k0,l,TJ,W1) * GMTR_a(ijm1  ,k0,l,AJ ,hn) & ! Q6 * b5
                                  + GMTR_t(ijm1  ,k0,l,TJ,W1) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q6 * b6
                                ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
       enddo
       !$omp end do

       if ( ADM_have_sgp(l) ) then ! pentagon
          !$omp master
          ij     = gmin
          ip1j   = gmin + 1
          ip1jp1 = gmin + iall + 1
          ijp1   = gmin + iall
          im1j   = gmin - 1
          im1jm1 = gmin - iall - 1
          ijm1   = gmin - iall

          ! ij
          coef_grad(ij,0,d,l) = ( + GMTR_t(ij    ,k0,l,TI,W1) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q1 * b6
                                  + GMTR_t(ij    ,k0,l,TI,W1) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q1 * b1
                                  + GMTR_t(ij    ,k0,l,TJ,W1) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q2 * b1
                                  + GMTR_t(ij    ,k0,l,TJ,W1) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q2 * b2
                                  + GMTR_t(im1j  ,k0,l,TI,W2) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q3 * b2
                                  - GMTR_t(im1j  ,k0,l,TI,W2) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q3 * b3
                                  - GMTR_t(im1jm1,k0,l,TJ,W2) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q4 * b3
                                  - GMTR_t(im1jm1,k0,l,TJ,W2) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q4 * b4
                                  - GMTR_t(ijm1  ,k0,l,TJ,W3) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q6 * b4
                                  + GMTR_t(ijm1  ,k0,l,TJ,W3) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q6 * b6
                                  - 2.0_RP * GMTR_a(ij    ,k0,l,AIJ,hn)                    & ! P0 * b1
                                  - 2.0_RP * GMTR_a(ij    ,k0,l,AJ ,hn)                    & ! P0 * b2
                                  + 2.0_RP * GMTR_a(im1j  ,k0,l,AI ,hn)                    & ! P0 * b3
                                  + 2.0_RP * GMTR_a(im1jm1,k0,l,AIJ,hn)                    & ! P0 * b4
                                  - 2.0_RP * GMTR_a(ij    ,k0,l,AI ,hn)                    & ! P0 * b6
                                ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ip1j
          coef_grad(ij,1,d,l) = ( - GMTR_t(ijm1  ,k0,l,TJ,W2) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q6 * b4
                                  + GMTR_t(ijm1  ,k0,l,TJ,W2) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q6 * b6
                                  + GMTR_t(ij    ,k0,l,TI,W2) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q1 * b6
                                  + GMTR_t(ij    ,k0,l,TI,W2) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q1 * b1
                                ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ip1jp1
          coef_grad(ij,2,d,l) = ( + GMTR_t(ij    ,k0,l,TI,W3) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q1 * b6
                                  + GMTR_t(ij    ,k0,l,TI,W3) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q1 * b1
                                  + GMTR_t(ij    ,k0,l,TJ,W2) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q2 * b1
                                  + GMTR_t(ij    ,k0,l,TJ,W2) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q2 * b2
                                ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ijp1
          coef_grad(ij,3,d,l) = ( + GMTR_t(ij    ,k0,l,TJ,W3) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q2 * b1
                                  + GMTR_t(ij    ,k0,l,TJ,W3) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q2 * b2
                                  + GMTR_t(im1j  ,k0,l,TI,W3) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q3 * b2
                                  - GMTR_t(im1j  ,k0,l,TI,W3) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q3 * b3
                                ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! im1j
          coef_grad(ij,4,d,l) = ( + GMTR_t(im1j  ,k0,l,TI,W1) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q3 * b2
                                  - GMTR_t(im1j  ,k0,l,TI,W1) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q3 * b3
                                  - GMTR_t(im1jm1,k0,l,TJ,W3) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q4 * b3
                                  - GMTR_t(im1jm1,k0,l,TJ,W3) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q4 * b4
                                ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! im1jm1
          coef_grad(ij,5,d,l) = ( - GMTR_t(im1jm1,k0,l,TJ,W1) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q4 * b3
                                  - GMTR_t(im1jm1,k0,l,TJ,W1) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q4 * b4
                                ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ijm1
          coef_grad(ij,6,d,l) = ( - GMTR_t(ijm1  ,k0,l,TJ,W1) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q6 * b4
                                  + GMTR_t(ijm1  ,k0,l,TJ,W1) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q6 * b6
                                ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          !$omp end master
       endif

    enddo ! loop d
    enddo ! loop l
    !$omp end parallel

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
       do d = 1, ADM_nxyz
          hn = d + HNX - 1

          coef = 0.0_RP
          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             if( ijp1 == ADM_gmax_pl + 1 ) ijp1 = ADM_gmin_pl

             coef = coef + 2.0_RP * ( GMTR_t_pl(ij,k0,l,W1) - 1.0_RP ) * GMTR_a_pl(ijp1,k0,l,hn)
          enddo
          coef_grad_pl(0,d,l) = coef * 0.5_RP * GMTR_p_pl(n,k0,l,P_RAREA)

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             ijm1 = v - 1
             if( ijp1 == ADM_gmax_pl + 1 ) ijp1 = ADM_gmin_pl
             if( ijm1 == ADM_gmin_pl - 1 ) ijm1 = ADM_gmax_pl

             coef_grad_pl(v-1,d,l) = ( + GMTR_t_pl(ijm1,k0,l,W3) * GMTR_a_pl(ijm1,k0,l,hn) &
                                       + GMTR_t_pl(ijm1,k0,l,W3) * GMTR_a_pl(ij  ,k0,l,hn) &
                                       + GMTR_t_pl(ij  ,k0,l,W2) * GMTR_a_pl(ij  ,k0,l,hn) &
                                       + GMTR_t_pl(ij  ,k0,l,W2) * GMTR_a_pl(ijp1,k0,l,hn) &
                                     ) * 0.5_RP * GMTR_p_pl(n,k0,l,P_RAREA)
          enddo
       enddo ! loop d
       enddo ! loop l
    endif

    return
  end subroutine OPRT_gradient_setup

  !-----------------------------------------------------------------------------
  subroutine OPRT_laplacian_setup( &
       GMTR_p,   GMTR_p_pl,  &
       GMTR_t,   GMTR_t_pl,  &
       GMTR_a,   GMTR_a_pl,  &
       coef_lap, coef_lap_pl )
!ESC!    use mod_adm, only: &
!ESC!       ADM_have_pl,  &
!ESC!       ADM_have_sgp, &
!ESC!       ADM_gall_1d,  &
!ESC!       ADM_gmin,     &
!ESC!       ADM_gmax,     &
!ESC!       ADM_gslf_pl,  &
!ESC!       ADM_gmin_pl,  &
!ESC!       ADM_gmax_pl
!ESC!    use mod_gmtr, only: &
!ESC!       P_RAREA => GMTR_p_RAREA, &
!ESC!       T_RAREA => GMTR_t_RAREA, &
!ESC!       HNX     => GMTR_a_HNX,   &
!ESC!       TNX     => GMTR_a_TNX,   &
!ESC!       TN2X    => GMTR_a_TN2X,  &
!ESC!       GMTR_p_nmax,             &
!ESC!       GMTR_t_nmax,             &
!ESC!       GMTR_a_nmax,             &
!ESC!       GMTR_a_nmax_pl
    implicit none

    real(RP), intent(in)  :: GMTR_p     (ADM_gall   ,K0,ADM_lall   ,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_p_pl  (ADM_gall_pl,K0,ADM_lall_pl,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_t     (ADM_gall   ,K0,ADM_lall   ,TI:TJ,GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_t_pl  (ADM_gall_pl,K0,ADM_lall_pl,      GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_a     (ADM_gall   ,K0,ADM_lall   ,AI:AJ,GMTR_a_nmax   )
    real(RP), intent(in)  :: GMTR_a_pl  (ADM_gall_pl,K0,ADM_lall_pl,      GMTR_a_nmax_pl)
    real(RP), intent(out) :: coef_lap   (ADM_gall,0:6        ,ADM_lall   )
    real(RP), intent(out) :: coef_lap_pl(         0:ADM_vlink,ADM_lall_pl)

    integer  :: gmin, gmax, iall, gall, nxyz, lall

    integer  :: ij
    integer  :: ip1j, ijp1, ip1jp1
    integer  :: im1j, ijm1, im1jm1

    integer  :: g, l, d, n, v, hn, tn, tn2
    !---------------------------------------------------------------------------

    !if( IO_L ) write(IO_FID_LOG,*) '*** setup coefficient of laplacian operator'

    gmin = (ADM_gmin-1)*ADM_gall_1d + ADM_gmin
    gmax = (ADM_gmax-1)*ADM_gall_1d + ADM_gmax
    iall = ADM_gall_1d
    gall = ADM_gall
    nxyz = ADM_nxyz
    lall = ADM_lall

    !$omp parallel workshare
    coef_lap   (:,:,:) = 0.0_RP
    !$omp end parallel workshare
    coef_lap_pl(  :,:) = 0.0_RP

    !$omp parallel default(none),private(g,d,l,hn,tn,ij,ip1j,ip1jp1,ijp1,im1j,ijm1,im1jm1), &
    !$omp shared(ADM_have_sgp,gmin,gmax,iall,gall,nxyz,lall,coef_lap,GMTR_p,GMTR_t,GMTR_a)
    do l = 1, lall

       do d = 1, nxyz
          hn = d + HNX - 1
          tn = d + TNX - 1

          !$omp do
          do g = gmin, gmax
             ij     = g
             ip1j   = g + 1
             ip1jp1 = g + iall + 1
             ijp1   = g + iall
             im1j   = g - 1
             im1jm1 = g - iall - 1
             ijm1   = g - iall

             ! ij
             coef_lap(ij,0,l) = coef_lap(ij,0,l) &
                              + GMTR_t(ij    ,k0,l,TI,T_RAREA) &
                              * ( - 1.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  + 2.0_RP * GMTR_a(ip1j  ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  - 1.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  + 2.0_RP * GMTR_a(ip1j  ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) )

             coef_lap(ij,0,l) = coef_lap(ij,0,l) &
                              + GMTR_t(ij    ,k0,l,TJ,T_RAREA) &
                              * ( - 1.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  - 2.0_RP * GMTR_a(ijp1  ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  - 1.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  - 2.0_RP * GMTR_a(ijp1  ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) )

             coef_lap(ij,0,l) = coef_lap(ij,0,l) &
                              + GMTR_t(im1j  ,k0,l,TI,T_RAREA) &
                              * ( - 1.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  - 2.0_RP * GMTR_a(im1j  ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  - 1.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  + 2.0_RP * GMTR_a(im1j  ,k0,l,AIJ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  + 1.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) )

             coef_lap(ij,0,l) = coef_lap(ij,0,l) &
                              + GMTR_t(im1jm1,k0,l,TJ,T_RAREA) &
                              * ( - 1.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  + 2.0_RP * GMTR_a(im1jm1,k0,l,AJ ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  + 1.0_RP * GMTR_a(im1jm1,k0,l,AIJ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  - 1.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  + 2.0_RP * GMTR_a(im1jm1,k0,l,AJ ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  + 1.0_RP * GMTR_a(im1jm1,k0,l,AIJ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) )

             coef_lap(ij,0,l) = coef_lap(ij,0,l) &
                              + GMTR_t(im1jm1,k0,l,TI,T_RAREA) &
                              * ( - 1.0_RP * GMTR_a(im1jm1,k0,l,AIJ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  - 2.0_RP * GMTR_a(im1jm1,k0,l,AI ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  + 1.0_RP * GMTR_a(ijm1  ,k0,l,AJ ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  - 1.0_RP * GMTR_a(im1jm1,k0,l,AIJ,tn) * GMTR_a(ijm1  ,k0,l,AJ ,hn) &
                                  - 2.0_RP * GMTR_a(im1jm1,k0,l,AI ,tn) * GMTR_a(ijm1  ,k0,l,AJ ,hn) &
                                  + 1.0_RP * GMTR_a(ijm1  ,k0,l,AJ ,tn) * GMTR_a(ijm1  ,k0,l,AJ ,hn) )

             coef_lap(ij,0,l) = coef_lap(ij,0,l) &
                              + GMTR_t(ijm1  ,k0,l,TJ,T_RAREA) &
                              * ( - 1.0_RP * GMTR_a(ijm1  ,k0,l,AJ ,tn) * GMTR_a(ijm1  ,k0,l,AJ ,hn) &
                                  - 2.0_RP * GMTR_a(ijm1  ,k0,l,AIJ,tn) * GMTR_a(ijm1  ,k0,l,AJ ,hn) &
                                  - 1.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ijm1  ,k0,l,AJ ,hn) &
                                  + 1.0_RP * GMTR_a(ijm1  ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  + 2.0_RP * GMTR_a(ijm1  ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) )

             ! ip1j
             coef_lap(ij,1,l) = coef_lap(ij,1,l) &
                              + GMTR_t(ijm1  ,k0,l,TJ,T_RAREA) &
                              * ( - 1.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ijm1  ,k0,l,AJ ,hn) &
                                  + 2.0_RP * GMTR_a(ijm1  ,k0,l,AJ ,tn) * GMTR_a(ijm1  ,k0,l,AJ ,hn) &
                                  + 1.0_RP * GMTR_a(ijm1  ,k0,l,AIJ,tn) * GMTR_a(ijm1  ,k0,l,AJ ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  - 2.0_RP * GMTR_a(ijm1  ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  - 1.0_RP * GMTR_a(ijm1  ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) )

             coef_lap(ij,1,l) = coef_lap(ij,1,l) &
                              + GMTR_t(ij    ,k0,l,TI,T_RAREA) &
                              * ( - 1.0_RP * GMTR_a(ip1j  ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  - 2.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  - 1.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  - 1.0_RP * GMTR_a(ip1j  ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  - 2.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  - 1.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) )

             ! ip1jp1
             coef_lap(ij,2,l) = coef_lap(ij,2,l) &
                              + GMTR_t(ij    ,k0,l,TI,T_RAREA) &
                              * ( + 1.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  + 2.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  - 1.0_RP * GMTR_a(ip1j  ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  + 2.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  - 1.0_RP * GMTR_a(ip1j  ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) )

             coef_lap(ij,2,l) = coef_lap(ij,2,l) &
                              + GMTR_t(ij    ,k0,l,TJ,T_RAREA) &
                              * ( + 1.0_RP * GMTR_a(ijp1  ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  - 2.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  - 1.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  + 1.0_RP * GMTR_a(ijp1  ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  - 1.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  - 2.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) )

             ! ijp1
             coef_lap(ij,3,l) = coef_lap(ij,3,l) &
                              + GMTR_t(ij    ,k0,l,TJ,T_RAREA) &
                              * ( + 1.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  + 2.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  + 1.0_RP * GMTR_a(ijp1  ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  + 2.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  + 1.0_RP * GMTR_a(ijp1  ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) )

             coef_lap(ij,3,l) = coef_lap(ij,3,l) &
                              + GMTR_t(im1j  ,k0,l,TI,T_RAREA) &
                              * ( + 1.0_RP * GMTR_a(im1j  ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  + 2.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  - 1.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  - 1.0_RP * GMTR_a(im1j  ,k0,l,AIJ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  - 2.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) )

             ! im1j
             coef_lap(ij,4,l) = coef_lap(ij,4,l) &
                              + GMTR_t(im1j  ,k0,l,TI,T_RAREA) &
                              * ( - 1.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  + 2.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  + 1.0_RP * GMTR_a(im1j  ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  + 1.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  - 2.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  - 1.0_RP * GMTR_a(im1j  ,k0,l,AIJ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) )

             coef_lap(ij,4,l) = coef_lap(ij,4,l) &
                              + GMTR_t(im1jm1,k0,l,TJ,T_RAREA) &
                              * ( - 1.0_RP * GMTR_a(im1jm1,k0,l,AJ ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  - 2.0_RP * GMTR_a(im1jm1,k0,l,AIJ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  - 1.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  - 1.0_RP * GMTR_a(im1jm1,k0,l,AJ ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  - 2.0_RP * GMTR_a(im1jm1,k0,l,AIJ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  - 1.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) )

             ! im1jm1
             coef_lap(ij,5,l) = coef_lap(ij,5,l) &
                              + GMTR_t(im1jm1,k0,l,TJ,T_RAREA) &
                              * ( + 1.0_RP * GMTR_a(im1jm1,k0,l,AIJ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  + 2.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  - 1.0_RP * GMTR_a(im1jm1,k0,l,AJ ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  + 1.0_RP * GMTR_a(im1jm1,k0,l,AIJ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  + 2.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  - 1.0_RP * GMTR_a(im1jm1,k0,l,AJ ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) )

             coef_lap(ij,5,l) = coef_lap(ij,5,l) &
                              + GMTR_t(im1jm1,k0,l,TI,T_RAREA) &
                              * ( + 1.0_RP * GMTR_a(im1jm1,k0,l,AI ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  - 2.0_RP * GMTR_a(ijm1  ,k0,l,AJ ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  - 1.0_RP * GMTR_a(im1jm1,k0,l,AIJ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  + 1.0_RP * GMTR_a(im1jm1,k0,l,AI ,tn) * GMTR_a(ijm1  ,k0,l,AJ ,hn) &
                                  - 2.0_RP * GMTR_a(ijm1  ,k0,l,AJ ,tn) * GMTR_a(ijm1  ,k0,l,AJ ,hn) &
                                  - 1.0_RP * GMTR_a(im1jm1,k0,l,AIJ,tn) * GMTR_a(ijm1  ,k0,l,AJ ,hn) )

             ! ijm1
             coef_lap(ij,6,l) = coef_lap(ij,6,l) &
                              + GMTR_t(im1jm1,k0,l,TI,T_RAREA) &
                              * ( + 1.0_RP * GMTR_a(ijm1  ,k0,l,AJ ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  + 2.0_RP * GMTR_a(im1jm1,k0,l,AIJ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  + 1.0_RP * GMTR_a(im1jm1,k0,l,AI ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  + 1.0_RP * GMTR_a(ijm1  ,k0,l,AJ ,tn) * GMTR_a(ijm1  ,k0,l,AJ ,hn) &
                                  + 2.0_RP * GMTR_a(im1jm1,k0,l,AIJ,tn) * GMTR_a(ijm1  ,k0,l,AJ ,hn) &
                                  + 1.0_RP * GMTR_a(im1jm1,k0,l,AI ,tn) * GMTR_a(ijm1  ,k0,l,AJ ,hn) )

             coef_lap(ij,6,l) = coef_lap(ij,6,l) &
                              + GMTR_t(ijm1  ,k0,l,TJ,T_RAREA) &
                              * ( + 1.0_RP * GMTR_a(ijm1  ,k0,l,AIJ,tn) * GMTR_a(ijm1  ,k0,l,AJ ,hn) &
                                  + 2.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ijm1  ,k0,l,AJ ,hn) &
                                  - 1.0_RP * GMTR_a(ijm1  ,k0,l,AJ ,tn) * GMTR_a(ijm1  ,k0,l,AJ ,hn) &
                                  - 1.0_RP * GMTR_a(ijm1  ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  - 2.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  + 1.0_RP * GMTR_a(ijm1  ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) )
          enddo
          !$omp end do
       enddo ! loop d

       if ( ADM_have_sgp(l) ) then ! pentagon
          !$omp master
          coef_lap(gmin,0,l) = 0.0_RP
          coef_lap(gmin,1,l) = 0.0_RP
          coef_lap(gmin,2,l) = 0.0_RP
          coef_lap(gmin,3,l) = 0.0_RP
          coef_lap(gmin,4,l) = 0.0_RP
          coef_lap(gmin,5,l) = 0.0_RP
          coef_lap(gmin,6,l) = 0.0_RP

          do d = 1, ADM_nxyz
             hn = d + HNX - 1
             tn = d + TNX - 1

             ij     = gmin
             ip1j   = gmin + 1
             ip1jp1 = gmin + iall + 1
             ijp1   = gmin + iall
             im1j   = gmin - 1
             im1jm1 = gmin - iall - 1
             ijm1   = gmin - iall

             ! ij
             coef_lap(ij,0,l) = coef_lap(ij,0,l) &
                              + GMTR_t(ij    ,k0,l,TI,T_RAREA) &
                              * ( - 1.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  + 2.0_RP * GMTR_a(ip1j  ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  - 1.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  + 2.0_RP * GMTR_a(ip1j  ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) )

             coef_lap(ij,0,l) = coef_lap(ij,0,l) &
                              + GMTR_t(ij    ,k0,l,TJ,T_RAREA) &
                              * ( - 1.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  - 2.0_RP * GMTR_a(ijp1  ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  - 1.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  - 2.0_RP * GMTR_a(ijp1  ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) )

             coef_lap(ij,0,l) = coef_lap(ij,0,l) &
                              + GMTR_t(im1j  ,k0,l,TI,T_RAREA) &
                              * ( - 1.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  - 2.0_RP * GMTR_a(im1j  ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  - 1.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  + 2.0_RP * GMTR_a(im1j  ,k0,l,AIJ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  + 1.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) )

             coef_lap(ij,0,l) = coef_lap(ij,0,l) &
                              + GMTR_t(im1jm1,k0,l,TJ,T_RAREA) &
                              * ( - 1.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  + 2.0_RP * GMTR_a(im1jm1,k0,l,AJ ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  + 1.0_RP * GMTR_a(im1jm1,k0,l,AIJ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  - 1.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  + 2.0_RP * GMTR_a(im1jm1,k0,l,AJ ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  + 1.0_RP * GMTR_a(im1jm1,k0,l,AIJ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) )

             coef_lap(ij,0,l) = coef_lap(ij,0,l) &
                              + GMTR_t(ijm1  ,k0,l,TJ,T_RAREA) &
                              * ( - 1.0_RP * GMTR_a(ijm1  ,k0,l,AJ ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  - 2.0_RP * GMTR_a(ijm1  ,k0,l,AIJ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  - 1.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  + 1.0_RP * GMTR_a(ijm1  ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  + 2.0_RP * GMTR_a(ijm1  ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) )

             ! ip1j
             coef_lap(ij,1,l) = coef_lap(ij,1,l) &
                              + GMTR_t(ijm1  ,k0,l,TJ,T_RAREA) &
                              * ( + 1.0_RP * GMTR_a(ijm1  ,k0,l,AIJ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  + 2.0_RP * GMTR_a(ijm1  ,k0,l,AJ ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  - 1.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  - 1.0_RP * GMTR_a(ijm1  ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  - 2.0_RP * GMTR_a(ijm1  ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) )

             coef_lap(ij,1,l) = coef_lap(ij,1,l) &
                              + GMTR_t(ij    ,k0,l,TI,T_RAREA) &
                              * ( - 1.0_RP * GMTR_a(ip1j  ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  - 2.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  - 1.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  - 1.0_RP * GMTR_a(ip1j  ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  - 2.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  - 1.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) )

             ! ip1jp1
             coef_lap(ij,2,l) = coef_lap(ij,2,l) &
                              + GMTR_t(ij    ,k0,l,TI,T_RAREA) &
                              * ( + 1.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  + 2.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  - 1.0_RP * GMTR_a(ip1j  ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  + 2.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  - 1.0_RP * GMTR_a(ip1j  ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) )

             coef_lap(ij,2,l) = coef_lap(ij,2,l) &
                              + GMTR_t(ij    ,k0,l,TJ,T_RAREA) &
                              * ( + 1.0_RP * GMTR_a(ijp1  ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  - 2.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  - 1.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  + 1.0_RP * GMTR_a(ijp1  ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  - 2.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  - 1.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) )

             ! ijp1
             coef_lap(ij,3,l) = coef_lap(ij,3,l) &
                              + GMTR_t(ij    ,k0,l,TJ,T_RAREA) &
                              * ( + 1.0_RP * GMTR_a(ijp1  ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  + 2.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  + 1.0_RP * GMTR_a(ijp1  ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  + 2.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) )

             coef_lap(ij,3,l) = coef_lap(ij,3,l) &
                              + GMTR_t(im1j  ,k0,l,TI,T_RAREA) &
                              * ( + 1.0_RP * GMTR_a(im1j  ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  + 2.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  - 1.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  - 1.0_RP * GMTR_a(im1j  ,k0,l,AIJ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  - 2.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) )

             ! im1j
             coef_lap(ij,4,l) = coef_lap(ij,4,l) &
                              + GMTR_t(im1j  ,k0,l,TI,T_RAREA) &
                              * ( + 1.0_RP * GMTR_a(im1j  ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  + 2.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  - 1.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  - 1.0_RP * GMTR_a(im1j  ,k0,l,AIJ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  - 2.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  + 1.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) )

             coef_lap(ij,4,l) = coef_lap(ij,4,l) &
                              + GMTR_t(im1jm1,k0,l,TJ,T_RAREA) &
                              * ( - 1.0_RP * GMTR_a(im1jm1,k0,l,AJ ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  - 2.0_RP * GMTR_a(im1jm1,k0,l,AIJ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  - 1.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  - 1.0_RP * GMTR_a(im1jm1,k0,l,AJ ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  - 2.0_RP * GMTR_a(im1jm1,k0,l,AIJ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  - 1.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) )

             ! im1jm1
             coef_lap(ij,5,l) = coef_lap(ij,5,l) &
                              + GMTR_t(im1jm1,k0,l,TJ,T_RAREA) &
                              * ( - 1.0_RP * GMTR_a(im1jm1,k0,l,AJ ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  + 2.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  + 1.0_RP * GMTR_a(im1jm1,k0,l,AIJ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  - 1.0_RP * GMTR_a(im1jm1,k0,l,AJ ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  + 2.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  + 1.0_RP * GMTR_a(im1jm1,k0,l,AIJ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) )

             ! ijm1
             coef_lap(ij,6,l) = coef_lap(ij,6,l) &
                              + GMTR_t(ijm1  ,k0,l,TJ,T_RAREA) &
                              * ( + 1.0_RP * GMTR_a(ijm1  ,k0,l,AIJ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  + 2.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  - 1.0_RP * GMTR_a(ijm1  ,k0,l,AJ ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  - 1.0_RP * GMTR_a(ijm1  ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  - 2.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  + 1.0_RP * GMTR_a(ijm1  ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) )
          enddo ! loop d
          !$omp end master
       endif

       !$omp do
       do g = 1, gall
          coef_lap(g,0,l) = coef_lap(g,0,l) * GMTR_p(g,k0,l,P_RAREA) / 12.0_RP
          coef_lap(g,1,l) = coef_lap(g,1,l) * GMTR_p(g,k0,l,P_RAREA) / 12.0_RP
          coef_lap(g,2,l) = coef_lap(g,2,l) * GMTR_p(g,k0,l,P_RAREA) / 12.0_RP
          coef_lap(g,3,l) = coef_lap(g,3,l) * GMTR_p(g,k0,l,P_RAREA) / 12.0_RP
          coef_lap(g,4,l) = coef_lap(g,4,l) * GMTR_p(g,k0,l,P_RAREA) / 12.0_RP
          coef_lap(g,5,l) = coef_lap(g,5,l) * GMTR_p(g,k0,l,P_RAREA) / 12.0_RP
          coef_lap(g,6,l) = coef_lap(g,6,l) * GMTR_p(g,k0,l,P_RAREA) / 12.0_RP
       enddo
       !$omp end do

    enddo ! loop l
    !$omp end parallel

    if ( ADM_have_pl ) then
       n  = ADM_gslf_pl

       do l = 1,ADM_lall_pl

          do d = 1, ADM_nxyz
             hn  = d + HNX  - 1
             tn  = d + TNX  - 1
             tn2 = d + TN2X - 1

             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijp1 = v + 1
                ijm1 = v - 1
                if( ijp1 == ADM_gmax_pl + 1 ) ijp1 = ADM_gmin_pl
                if( ijm1 == ADM_gmin_pl - 1 ) ijm1 = ADM_gmax_pl

                coef_lap_pl(0,l) = coef_lap_pl(0,l) &
                                 + GMTR_t_pl(ijm1,k0,l,T_RAREA) &
                                 * ( + 1.0_RP * GMTR_a_pl(ijm1,k0,l,tn ) * GMTR_a_pl(ij,k0,l,hn) &
                                     - 2.0_RP * GMTR_a_pl(ijm1,k0,l,tn2) * GMTR_a_pl(ij,k0,l,hn) &
                                     - 1.0_RP * GMTR_a_pl(ij  ,k0,l,tn ) * GMTR_a_pl(ij,k0,l,hn) )

                coef_lap_pl(0,l) = coef_lap_pl(0,l) &
                                 + GMTR_t_pl(ij  ,k0,l,T_RAREA) &
                                 * ( + 1.0_RP * GMTR_a_pl(ij  ,k0,l,tn ) * GMTR_a_pl(ij,k0,l,hn) &
                                     - 2.0_RP * GMTR_a_pl(ij  ,k0,l,tn2) * GMTR_a_pl(ij,k0,l,hn) &
                                     - 1.0_RP * GMTR_a_pl(ijp1,k0,l,tn ) * GMTR_a_pl(ij,k0,l,hn) )
             enddo

             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijp1 = v + 1
                ijm1 = v - 1
                if( ijp1 == ADM_gmax_pl + 1 ) ijp1 = ADM_gmin_pl
                if( ijm1 == ADM_gmin_pl - 1 ) ijm1 = ADM_gmax_pl

                coef_lap_pl(v-1,l) = coef_lap_pl(v-1,l) &
                                   + GMTR_t_pl(ijm1,k0,l,T_RAREA) &
                                   * ( - 2.0_RP * GMTR_a_pl(ijm1,k0,l,tn ) * GMTR_a_pl(ijm1,k0,l,hn) &
                                       + 1.0_RP * GMTR_a_pl(ijm1,k0,l,tn2) * GMTR_a_pl(ijm1,k0,l,hn) &
                                       - 1.0_RP * GMTR_a_pl(ij  ,k0,l,tn ) * GMTR_a_pl(ijm1,k0,l,hn) &
                                       - 2.0_RP * GMTR_a_pl(ijm1,k0,l,tn ) * GMTR_a_pl(ij  ,k0,l,hn) &
                                       + 1.0_RP * GMTR_a_pl(ijm1,k0,l,tn2) * GMTR_a_pl(ij  ,k0,l,hn) &
                                       - 1.0_RP * GMTR_a_pl(ij  ,k0,l,tn ) * GMTR_a_pl(ij  ,k0,l,hn) )

                coef_lap_pl(v-1,l) = coef_lap_pl(v-1,l) &
                                   + GMTR_t_pl(ij  ,k0,l,T_RAREA) &
                                   * ( + 1.0_RP * GMTR_a_pl(ij  ,k0,l,tn ) * GMTR_a_pl(ij  ,k0,l,hn) &
                                       + 1.0_RP * GMTR_a_pl(ij  ,k0,l,tn2) * GMTR_a_pl(ij  ,k0,l,hn) &
                                       + 2.0_RP * GMTR_a_pl(ijp1,k0,l,tn ) * GMTR_a_pl(ij  ,k0,l,hn) &
                                       + 1.0_RP * GMTR_a_pl(ij  ,k0,l,tn ) * GMTR_a_pl(ijp1,k0,l,hn) &
                                       + 1.0_RP * GMTR_a_pl(ij  ,k0,l,tn2) * GMTR_a_pl(ijp1,k0,l,hn) &
                                       + 2.0_RP * GMTR_a_pl(ijp1,k0,l,tn ) * GMTR_a_pl(ijp1,k0,l,hn) )
             enddo
          enddo ! d loop

          do v = ADM_gslf_pl, ADM_gmax_pl
             coef_lap_pl(v-1,l) = coef_lap_pl(v-1,l) * GMTR_p_pl(n,k0,l,P_RAREA) / 12.0_RP
          enddo

       enddo ! l loop
    endif

    return
  end subroutine OPRT_laplacian_setup

  !-----------------------------------------------------------------------------
  subroutine OPRT_diffusion_setup( &
       GMTR_p,    GMTR_p_pl,    &
       GMTR_t,    GMTR_t_pl,    &
       GMTR_a,    GMTR_a_pl,    &
       coef_intp, coef_intp_pl, &
       coef_diff, coef_diff_pl  )
!ESC!    use mod_adm, only: &
!ESC!       ADM_have_pl,  &
!ESC!       ADM_have_sgp, &
!ESC!       ADM_gall_1d,  &
!ESC!       ADM_gmin,     &
!ESC!       ADM_gmax,     &
!ESC!       ADM_gslf_pl,  &
!ESC!       ADM_gmin_pl,  &
!ESC!       ADM_gmax_pl
!ESC!    use mod_gmtr, only: &
!ESC!       P_RAREA => GMTR_p_RAREA, &
!ESC!       T_RAREA => GMTR_t_RAREA, &
!ESC!       HNX     => GMTR_a_HNX,   &
!ESC!       TNX     => GMTR_a_TNX,   &
!ESC!       TN2X    => GMTR_a_TN2X,  &
!ESC!       GMTR_p_nmax,             &
!ESC!       GMTR_t_nmax,             &
!ESC!       GMTR_a_nmax,             &
!ESC!       GMTR_a_nmax_pl
    implicit none

    real(RP), intent(in)  :: GMTR_p      (ADM_gall   ,K0,ADM_lall   ,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_p_pl   (ADM_gall_pl,K0,ADM_lall_pl,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_t      (ADM_gall   ,K0,ADM_lall   ,TI:TJ,GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_t_pl   (ADM_gall_pl,K0,ADM_lall_pl,      GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_a      (ADM_gall   ,K0,ADM_lall   ,AI:AJ,GMTR_a_nmax   )
    real(RP), intent(in)  :: GMTR_a_pl   (ADM_gall_pl,K0,ADM_lall_pl,      GMTR_a_nmax_pl)
    real(RP), intent(out) :: coef_intp   (ADM_gall   ,1:3,ADM_nxyz,TI:TJ,ADM_lall   )
    real(RP), intent(out) :: coef_intp_pl(ADM_gall_pl,1:3,ADM_nxyz,      ADM_lall_pl)
    real(RP), intent(out) :: coef_diff   (ADM_gall,1:6        ,ADM_nxyz,ADM_lall   )
    real(RP), intent(out) :: coef_diff_pl(         1:ADM_vlink,ADM_nxyz,ADM_lall_pl)

    integer  :: gmin, gmax, iall, gall, nxyz, lall, gminm1

    integer  :: ij
    integer  :: ip1j, ijp1
    integer  :: im1j, ijm1, im1jm1

    integer  :: g, l, d, n, v, hn, tn, tn2
    !---------------------------------------------------------------------------

    !if( IO_L ) write(IO_FID_LOG,*) '*** setup coefficient of diffusion operator'

    gmin = (ADM_gmin-1)*ADM_gall_1d + ADM_gmin
    gmax = (ADM_gmax-1)*ADM_gall_1d + ADM_gmax
    iall = ADM_gall_1d
    gall = ADM_gall
    nxyz = ADM_nxyz
    lall = ADM_lall

    !$omp parallel workshare
    coef_intp   (:,:,:,:,:) = 0.0_RP
    coef_diff   (:,:,:,:)   = 0.0_RP
    !$omp end parallel workshare
    coef_intp_pl(:,:,:,  :) = 0.0_RP
    coef_diff_pl(  :,:,:)   = 0.0_RP

    gminm1 = (ADM_gmin-1-1)*ADM_gall_1d + ADM_gmin-1

    !$omp parallel do default(none),private(g,d,l,tn,ij,ip1j,ijp1), &
    !$omp shared(gminm1,gmax,iall,gall,nxyz,lall,coef_intp,GMTR_t,GMTR_a), &
    !$omp collapse(2)
    do l = 1, lall
    do d = 1, nxyz
       tn = d + TNX - 1

       do g = gminm1, gmax
          ij     = g
          ip1j   = g + 1
          ijp1   = g + iall

          coef_intp(ij,1,d,TI,l) = ( + GMTR_a(ij  ,k0,l,AIJ,tn) - GMTR_a(ij  ,k0,l,AI ,tn) ) &
                                 * 0.5_RP * GMTR_t(ij,k0,l,TI,T_RAREA)
          coef_intp(ij,2,d,TI,l) = ( - GMTR_a(ij  ,k0,l,AI ,tn) - GMTR_a(ip1j,k0,l,AJ ,tn) ) &
                                 * 0.5_RP * GMTR_t(ij,k0,l,TI,T_RAREA)
          coef_intp(ij,3,d,TI,l) = ( - GMTR_a(ip1j,k0,l,AJ ,tn) + GMTR_a(ij  ,k0,l,AIJ,tn) ) &
                                 * 0.5_RP * GMTR_t(ij,k0,l,TI,T_RAREA)

          coef_intp(ij,1,d,TJ,l) = ( + GMTR_a(ij  ,k0,l,AJ ,tn) - GMTR_a(ij  ,k0,l,AIJ,tn) ) &
                                 * 0.5_RP * GMTR_t(ij,k0,l,TJ,T_RAREA)
          coef_intp(ij,2,d,TJ,l) = ( - GMTR_a(ij  ,k0,l,AIJ,tn) + GMTR_a(ijp1,k0,l,AI ,tn) ) &
                                 * 0.5_RP * GMTR_t(ij,k0,l,TJ,T_RAREA)
          coef_intp(ij,3,d,TJ,l) = ( + GMTR_a(ijp1,k0,l,AI ,tn) + GMTR_a(ij  ,k0,l,AJ ,tn) ) &
                                 * 0.5_RP * GMTR_t(ij,k0,l,TJ,T_RAREA)
       enddo
    enddo ! loop d
    enddo ! loop l
    !$omp end parallel do

    !$omp parallel default(none),private(g,d,l,hn,ij,im1j,ijm1,im1jm1), &
    !$omp shared(ADM_have_sgp,gmin,gmax,iall,gall,nxyz,lall,coef_diff,GMTR_p,GMTR_a)
    do l = 1, lall
    do d = 1, nxyz
       hn = d + HNX - 1

       !$omp do
       do g = gmin, gmax
          ij     = g
          im1j   = g - 1
          im1jm1 = g - iall - 1
          ijm1   = g - iall

          coef_diff(ij,1,d,l) = + GMTR_a(ij    ,k0,l,AIJ,hn) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          coef_diff(ij,2,d,l) = + GMTR_a(ij    ,k0,l,AJ ,hn) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          coef_diff(ij,3,d,l) = - GMTR_a(im1j  ,k0,l,AI ,hn) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          coef_diff(ij,4,d,l) = - GMTR_a(im1jm1,k0,l,AIJ,hn) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          coef_diff(ij,5,d,l) = - GMTR_a(ijm1  ,k0,l,AJ ,hn) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          coef_diff(ij,6,d,l) = + GMTR_a(ij    ,k0,l,AI ,hn) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
       enddo
       !$omp end do

       if ( ADM_have_sgp(l) ) then ! pentagon
          !$omp master
          coef_diff(gmin,5,d,l) = 0.0_RP
          !$omp end master
       endif
    enddo ! loop d
    enddo ! loop l
    !$omp end parallel

    if ( ADM_have_pl ) then
       n  = ADM_gslf_pl

       do l = 1, ADM_lall_pl

          do d = 1, ADM_nxyz
             hn  = d + HNX  - 1
             tn  = d + TNX  - 1
             tn2 = d + TN2X - 1

             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijp1 = v + 1
                if( ijp1 == ADM_gmax_pl+1 ) ijp1 = ADM_gmin_pl

                coef_intp_pl(v,1,d,l) = - GMTR_a_pl(ijp1,k0,l,tn ) + GMTR_a_pl(ij  ,k0,l,tn )
                coef_intp_pl(v,2,d,l) = + GMTR_a_pl(ij  ,k0,l,tn ) + GMTR_a_pl(ij  ,k0,l,tn2)
                coef_intp_pl(v,3,d,l) = + GMTR_a_pl(ij  ,k0,l,tn2) - GMTR_a_pl(ijp1,k0,l,tn )

                coef_intp_pl(v,:,d,l) = coef_intp_pl(v,:,d,l) * 0.5_RP * GMTR_t_pl(v,k0,l,T_RAREA)

                coef_diff_pl(v-1,d,l) = GMTR_a_pl(v,k0,l,hn) * 0.5_RP * GMTR_p_pl(n,k0,l,P_RAREA)
             enddo
          enddo

       enddo ! l loop
    endif

    return
  end subroutine OPRT_diffusion_setup

end module mod_oprt
