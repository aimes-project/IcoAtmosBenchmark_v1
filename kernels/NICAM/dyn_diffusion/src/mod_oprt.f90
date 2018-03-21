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
  public :: OPRT_diffusion

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine OPRT_diffusion( &
       dscl,      dscl_pl,      &
       scl,       scl_pl,       &
       kh,        kh_pl,        &
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
!ESC!    use mod_grd, only: &
!ESC!       XDIR => GRD_XDIR, &
!ESC!       YDIR => GRD_YDIR, &
!ESC!       ZDIR => GRD_ZDIR
    implicit none

    real(RP), intent(out) :: dscl        (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: dscl_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: scl         (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: scl_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: kh          (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: kh_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: coef_intp   (ADM_gall   ,1:3,    ADM_nxyz,TI:TJ,ADM_lall   )
    real(RP), intent(in)  :: coef_intp_pl(ADM_gall_pl,1:3,    ADM_nxyz,      ADM_lall_pl)
    real(RP), intent(in)  :: coef_diff   (ADM_gall   ,1:6,    ADM_nxyz,      ADM_lall   )
    real(RP), intent(in)  :: coef_diff_pl(              1:ADM_vlink,ADM_nxyz,      ADM_lall_pl)

    real(RP) :: vt   (ADM_gall   ,ADM_nxyz,TI:TJ)
    real(RP) :: vt_pl(ADM_gall_pl,ADM_nxyz)
    real(RP) :: kf   (1:6)

    integer  :: gmin, gmax, iall, gall, kall, lall, nxyz, gminm1

    integer  :: ij
    integer  :: ip1j, ijp1, ip1jp1
    integer  :: im1j, ijm1, im1jm1

    integer  :: g, k, l, d, n, v
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('OPRT_diffusion')

    gmin = (ADM_gmin-1)*ADM_gall_1d + ADM_gmin
    gmax = (ADM_gmax-1)*ADM_gall_1d + ADM_gmax
    iall = ADM_gall_1d
    gall = ADM_gall
    kall = ADM_kall
    lall = ADM_lall
    nxyz = ADM_nxyz

    gminm1 = (ADM_gmin-1-1)*ADM_gall_1d + ADM_gmin-1

    !$omp parallel default(none),private(g,k,l,d,ij,ip1j,ip1jp1,ijp1,im1j,ijm1,im1jm1), &
    !$omp shared(ADM_have_sgp,gminm1,gmin,gmax,iall,gall,kall,lall,nxyz,dscl,scl,kh,kf,vt,coef_intp,coef_diff)
    do l = 1, lall
    do k = 1, kall

       do d = 1, nxyz
          !$omp do
          do g = gminm1, gmax
             ij     = g
             ip1j   = g + 1
             ip1jp1 = g + iall + 1
             ijp1   = g + iall

             vt(g,d,TI) = ( ( + 2.0_RP * coef_intp(g,1,d,TI,l) &
                              - 1.0_RP * coef_intp(g,2,d,TI,l) &
                              - 1.0_RP * coef_intp(g,3,d,TI,l) ) * scl(ij    ,k,l) &
                          + ( - 1.0_RP * coef_intp(g,1,d,TI,l) &
                              + 2.0_RP * coef_intp(g,2,d,TI,l) &
                              - 1.0_RP * coef_intp(g,3,d,TI,l) ) * scl(ip1j  ,k,l) &
                          + ( - 1.0_RP * coef_intp(g,1,d,TI,l) &
                              - 1.0_RP * coef_intp(g,2,d,TI,l) &
                              + 2.0_RP * coef_intp(g,3,d,TI,l) ) * scl(ip1jp1,k,l) &
                          ) / 3.0_RP
          enddo
          !$omp end do nowait

          !$omp do
          do g = gminm1, gmax
             ij     = g
             ip1j   = g + 1
             ip1jp1 = g + iall + 1
             ijp1   = g + iall

             vt(g,d,TJ) = ( ( + 2.0_RP * coef_intp(g,1,d,TJ,l) &
                              - 1.0_RP * coef_intp(g,2,d,TJ,l) &
                              - 1.0_RP * coef_intp(g,3,d,TJ,l) ) * scl(ij    ,k,l) &
                          + ( - 1.0_RP * coef_intp(g,1,d,TJ,l) &
                              + 2.0_RP * coef_intp(g,2,d,TJ,l) &
                              - 1.0_RP * coef_intp(g,3,d,TJ,l) ) * scl(ip1jp1,k,l) &
                          + ( - 1.0_RP * coef_intp(g,1,d,TJ,l) &
                              - 1.0_RP * coef_intp(g,2,d,TJ,l) &
                              + 2.0_RP * coef_intp(g,3,d,TJ,l) ) * scl(ijp1  ,k,l) &
                          ) / 3.0_RP
          enddo
          !$omp end do
       enddo

       if ( ADM_have_sgp(l) ) then ! pentagon
          !$omp master
          vt(gminm1,XDIR,TI) = vt(gminm1+1,XDIR,TJ)
          vt(gminm1,YDIR,TI) = vt(gminm1+1,YDIR,TJ)
          vt(gminm1,ZDIR,TI) = vt(gminm1+1,ZDIR,TJ)
          !$omp end master
       endif

!OCL XFILL
       !$omp do
       do g = 1, gmin-1
          dscl(g,k,l) = 0.0_RP
       enddo
       !$omp end do nowait

       !$omp do
       do g = gmin, gmax
          ij     = g
          ip1j   = g + 1
          ip1jp1 = g + iall + 1
          ijp1   = g + iall
          im1j   = g - 1
          im1jm1 = g - iall - 1
          ijm1   = g - iall

          kf(1) = 0.5_RP * ( kh(ij    ,k,l) + kh(ip1jp1,k,l) )
          kf(2) = 0.5_RP * ( kh(ij    ,k,l) + kh(ijp1  ,k,l) )
          kf(3) = 0.5_RP * ( kh(im1j  ,k,l) + kh(ij    ,k,l) )
          kf(4) = 0.5_RP * ( kh(im1jm1,k,l) + kh(ij    ,k,l) )
          kf(5) = 0.5_RP * ( kh(ijm1  ,k,l) + kh(ij    ,k,l) )
          kf(6) = 0.5_RP * ( kh(ij    ,k,l) + kh(ip1j  ,k,l) )

          dscl(g,k,l) = ( kf(1) * coef_diff(g,1,XDIR,l) * ( vt(ij    ,XDIR,TI) + vt(ij    ,XDIR,TJ) ) &
                        + kf(2) * coef_diff(g,2,XDIR,l) * ( vt(ij    ,XDIR,TJ) + vt(im1j  ,XDIR,TI) ) &
                        + kf(3) * coef_diff(g,3,XDIR,l) * ( vt(im1j  ,XDIR,TI) + vt(im1jm1,XDIR,TJ) ) &
                        + kf(4) * coef_diff(g,4,XDIR,l) * ( vt(im1jm1,XDIR,TJ) + vt(im1jm1,XDIR,TI) ) &
                        + kf(5) * coef_diff(g,5,XDIR,l) * ( vt(im1jm1,XDIR,TI) + vt(ijm1  ,XDIR,TJ) ) &
                        + kf(6) * coef_diff(g,6,XDIR,l) * ( vt(ijm1  ,XDIR,TJ) + vt(ij    ,XDIR,TI) ) )
       enddo
       !$omp end do

       !$omp do
       do g = gmin, gmax
          ij     = g
          ip1j   = g + 1
          ip1jp1 = g + iall + 1
          ijp1   = g + iall
          im1j   = g - 1
          im1jm1 = g - iall - 1
          ijm1   = g - iall

          kf(1) = 0.5_RP * ( kh(ij    ,k,l) + kh(ip1jp1,k,l) )
          kf(2) = 0.5_RP * ( kh(ij    ,k,l) + kh(ijp1  ,k,l) )
          kf(3) = 0.5_RP * ( kh(im1j  ,k,l) + kh(ij    ,k,l) )
          kf(4) = 0.5_RP * ( kh(im1jm1,k,l) + kh(ij    ,k,l) )
          kf(5) = 0.5_RP * ( kh(ijm1  ,k,l) + kh(ij    ,k,l) )
          kf(6) = 0.5_RP * ( kh(ij    ,k,l) + kh(ip1j  ,k,l) )

          dscl(g,k,l) = dscl(g,k,l) + ( kf(1) * coef_diff(g,1,YDIR,l) * ( vt(ij    ,XDIR,TI) + vt(ij    ,YDIR,TJ) ) &
                                      + kf(2) * coef_diff(g,2,YDIR,l) * ( vt(ij    ,XDIR,TJ) + vt(im1j  ,YDIR,TI) ) &
                                      + kf(3) * coef_diff(g,3,YDIR,l) * ( vt(im1j  ,XDIR,TI) + vt(im1jm1,YDIR,TJ) ) &
                                      + kf(4) * coef_diff(g,4,YDIR,l) * ( vt(im1jm1,XDIR,TJ) + vt(im1jm1,YDIR,TI) ) &
                                      + kf(5) * coef_diff(g,5,YDIR,l) * ( vt(im1jm1,XDIR,TI) + vt(ijm1  ,YDIR,TJ) ) &
                                      + kf(6) * coef_diff(g,6,YDIR,l) * ( vt(ijm1  ,XDIR,TJ) + vt(ij    ,YDIR,TI) ) )
       enddo
       !$omp end do

       !$omp do
       do g = gmin, gmax
          ij     = g
          ip1j   = g + 1
          ip1jp1 = g + iall + 1
          ijp1   = g + iall
          im1j   = g - 1
          im1jm1 = g - iall - 1
          ijm1   = g - iall

          kf(1) = 0.5_RP * ( kh(ij    ,k,l) + kh(ip1jp1,k,l) )
          kf(2) = 0.5_RP * ( kh(ij    ,k,l) + kh(ijp1  ,k,l) )
          kf(3) = 0.5_RP * ( kh(im1j  ,k,l) + kh(ij    ,k,l) )
          kf(4) = 0.5_RP * ( kh(im1jm1,k,l) + kh(ij    ,k,l) )
          kf(5) = 0.5_RP * ( kh(ijm1  ,k,l) + kh(ij    ,k,l) )
          kf(6) = 0.5_RP * ( kh(ij    ,k,l) + kh(ip1j  ,k,l) )

          dscl(g,k,l) = dscl(g,k,l) + ( kf(1) * coef_diff(g,1,ZDIR,l) * ( vt(ij    ,XDIR,TI) + vt(ij    ,ZDIR,TJ) ) &
                                      + kf(2) * coef_diff(g,2,ZDIR,l) * ( vt(ij    ,XDIR,TJ) + vt(im1j  ,ZDIR,TI) ) &
                                      + kf(3) * coef_diff(g,3,ZDIR,l) * ( vt(im1j  ,XDIR,TI) + vt(im1jm1,ZDIR,TJ) ) &
                                      + kf(4) * coef_diff(g,4,ZDIR,l) * ( vt(im1jm1,XDIR,TJ) + vt(im1jm1,ZDIR,TI) ) &
                                      + kf(5) * coef_diff(g,5,ZDIR,l) * ( vt(im1jm1,XDIR,TI) + vt(ijm1  ,ZDIR,TJ) ) &
                                      + kf(6) * coef_diff(g,6,ZDIR,l) * ( vt(ijm1  ,XDIR,TJ) + vt(ij    ,ZDIR,TI) ) )
       enddo
       !$omp end do nowait

!OCL XFILL
       !$omp do
       do g = gmax+1, gall
          dscl(g,k,l) = 0.0_RP
       enddo
       !$omp end do
    enddo ! loop k
    enddo ! loop l
    !$omp end parallel

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall

          do d = 1, ADM_nxyz
             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijp1 = v + 1
                if( ijp1 == ADM_gmax_pl+1 ) ijp1 = ADM_gmin_pl

                vt_pl(ij,d) = ( ( + 2.0_RP * coef_intp_pl(v,1,d,l) &
                                  - 1.0_RP * coef_intp_pl(v,2,d,l) &
                                  - 1.0_RP * coef_intp_pl(v,3,d,l) ) * scl_pl(n   ,k,l) &
                              + ( - 1.0_RP * coef_intp_pl(v,1,d,l) &
                                  + 2.0_RP * coef_intp_pl(v,2,d,l) &
                                  - 1.0_RP * coef_intp_pl(v,3,d,l) ) * scl_pl(ij  ,k,l) &
                              + ( - 1.0_RP * coef_intp_pl(v,1,d,l) &
                                  - 1.0_RP * coef_intp_pl(v,2,d,l) &
                                  + 2.0_RP * coef_intp_pl(v,3,d,l) ) * scl_pl(ijp1,k,l) &
                              ) / 3.0_RP
             enddo
          enddo

          dscl_pl(:,k,l) = 0.0_RP

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijm1 = v - 1
             if( ijm1 == ADM_gmin_pl-1 ) ijm1 = ADM_gmax_pl ! cyclic condition

             dscl_pl(n,k,l) = dscl_pl(n,k,l) &
                            + ( coef_diff_pl(v-1,XDIR,l) * ( vt_pl(ijm1,XDIR) + vt_pl(ij,XDIR) ) &
                              + coef_diff_pl(v-1,YDIR,l) * ( vt_pl(ijm1,YDIR) + vt_pl(ij,YDIR) ) &
                              + coef_diff_pl(v-1,ZDIR,l) * ( vt_pl(ijm1,ZDIR) + vt_pl(ij,ZDIR) ) &
                              ) * 0.5_RP * ( kh_pl(n,k,l) + kh_pl(ij,k,l) )
          enddo

       enddo
       enddo
    else
       dscl_pl(:,:,:) = 0.0_RP
    endif

    call DEBUG_rapend('OPRT_diffusion')

    return
  end subroutine OPRT_diffusion

end module mod_oprt
