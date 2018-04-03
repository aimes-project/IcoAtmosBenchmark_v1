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
  public :: vertical_limiter_thuburn

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
  subroutine vertical_limiter_thuburn( &
       q_h, q_h_pl, &
       q,   q_pl,   &
       d,   d_pl,   &
       ck,  ck_pl   )
!ESC!    use mod_const, only: &
!ESC!       CONST_HUGE, &
!ESC!       CONST_EPS
!ESC!    use mod_adm, only: &
!ESC!       ADM_have_pl, &
!ESC!       ADM_gall,    &
!ESC!       ADM_gall_pl, &
!ESC!       ADM_lall,    &
!ESC!       ADM_lall_pl, &
!ESC!       ADM_kall,    &
!ESC!       ADM_kmin,    &
!ESC!       ADM_kmax
    implicit none

    real(RP), intent(inout) :: q_h   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: q_h_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: q     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: q_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: d     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: d_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: ck    (ADM_gall   ,ADM_kall,ADM_lall   ,2)
    real(RP), intent(in)    :: ck_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,2)

    real(RP) :: Qout_min_k
    real(RP) :: Qout_max_k
    real(RP) :: Qout_min_km1(ADM_gall)
    real(RP) :: Qout_max_km1(ADM_gall)
    real(RP) :: Qout_min_pl(ADM_gall_pl,ADM_kall)
    real(RP) :: Qout_max_pl(ADM_gall_pl,ADM_kall)

    real(RP) :: Qin_minL, Qin_maxL
    real(RP) :: Qin_minU, Qin_maxU
    real(RP) :: qnext_min, qnext_max
    real(RP) :: Cin, Cout
    real(RP) :: CQin_min, CQin_max
    real(RP) :: inflagL, inflagU
    real(RP) :: zerosw

    integer  :: gall, kmin, kmax
    real(RP) :: EPS, BIG

    integer  :: g, k, l
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('____vertical_adv_limiter')

    gall = ADM_gall
    kmin = ADM_kmin
    kmax = ADM_kmax

    EPS  = CONST_EPS
    BIG  = CONST_HUGE

    do l = 1, ADM_lall
       !$omp parallel default(none), &
       !$omp private(g,k,zerosw,inflagL,inflagU,Qin_minL,Qin_minU,Qin_maxL,Qin_maxU, &
       !$omp         qnext_min,qnext_max,Cin,Cout,CQin_min,CQin_max,Qout_min_k,Qout_max_k),  &
       !$omp shared(l,gall,kmin,kmax,q_h,ck,q,d,Qout_min_km1,Qout_max_km1,EPS,BIG)

!OCL XFILL
       !$omp do
       do g = 1, gall
          k = kmin ! peeling

          inflagL = 0.5_RP - sign(0.5_RP,ck(g,k  ,l,1)) ! incoming flux: flag=1
          inflagU = 0.5_RP + sign(0.5_RP,ck(g,k+1,l,1)) ! incoming flux: flag=1

          Qin_minL = min( q(g,k,l), q(g,k-1,l) ) + ( 1.0_RP-inflagL ) * BIG
          Qin_minU = min( q(g,k,l), q(g,k+1,l) ) + ( 1.0_RP-inflagU ) * BIG
          Qin_maxL = max( q(g,k,l), q(g,k-1,l) ) - ( 1.0_RP-inflagL ) * BIG
          Qin_maxU = max( q(g,k,l), q(g,k+1,l) ) - ( 1.0_RP-inflagU ) * BIG

          qnext_min = min( Qin_minL, Qin_minU, q(g,k,l) )
          qnext_max = max( Qin_maxL, Qin_maxU, q(g,k,l) )

          Cin      = (        inflagL ) * ck(g,k,l,1) &
                   + (        inflagU ) * ck(g,k,l,2)
          Cout     = ( 1.0_RP-inflagL ) * ck(g,k,l,1) &
                   + ( 1.0_RP-inflagU ) * ck(g,k,l,2)

          CQin_min = (        inflagL ) * ck(g,k,l,1) * Qin_minL &
                   + (        inflagU ) * ck(g,k,l,2) * Qin_minU
          CQin_max = (        inflagL ) * ck(g,k,l,1) * Qin_maxL &
                   + (        inflagU ) * ck(g,k,l,2) * Qin_maxU

          zerosw = 0.5_RP - sign(0.5_RP,abs(Cout)-EPS) ! if Cout = 0, sw = 1

          Qout_min_k = ( ( q(g,k,l) - qnext_max ) + qnext_max*(Cin+Cout-d(g,k,l)) - CQin_max ) &
                     / ( Cout + zerosw ) * ( 1.0_RP - zerosw )                                 &
                     + q(g,k,l) * zerosw
          Qout_max_k = ( ( q(g,k,l) - qnext_min ) + qnext_min*(Cin+Cout-d(g,k,l)) - CQin_min ) &
                     / ( Cout + zerosw ) * ( 1.0_RP - zerosw )                                 &
                     + q(g,k,l) * zerosw

          Qout_min_km1(g) = Qout_min_k
          Qout_max_km1(g) = Qout_max_k
       enddo
       !$omp end do

       do k = kmin+1, kmax
!OCL XFILL
          !$omp do
          do g = 1, gall
             inflagL = 0.5_RP - sign(0.5_RP,ck(g,k  ,l,1)) ! incoming flux: flag=1
             inflagU = 0.5_RP + sign(0.5_RP,ck(g,k+1,l,1)) ! incoming flux: flag=1

             Qin_minL = min( q(g,k,l), q(g,k-1,l) ) + ( 1.0_RP-inflagL ) * BIG
             Qin_minU = min( q(g,k,l), q(g,k+1,l) ) + ( 1.0_RP-inflagU ) * BIG
             Qin_maxL = max( q(g,k,l), q(g,k-1,l) ) - ( 1.0_RP-inflagL ) * BIG
             Qin_maxU = max( q(g,k,l), q(g,k+1,l) ) - ( 1.0_RP-inflagU ) * BIG

             qnext_min = min( Qin_minL, Qin_minU, q(g,k,l) )
             qnext_max = max( Qin_maxL, Qin_maxU, q(g,k,l) )

             Cin      = (        inflagL ) * ck(g,k,l,1) &
                      + (        inflagU ) * ck(g,k,l,2)
             Cout     = ( 1.0_RP-inflagL ) * ck(g,k,l,1) &
                      + ( 1.0_RP-inflagU ) * ck(g,k,l,2)

             CQin_min = (        inflagL ) * ck(g,k,l,1) * Qin_minL &
                      + (        inflagU ) * ck(g,k,l,2) * Qin_minU
             CQin_max = (        inflagL ) * ck(g,k,l,1) * Qin_maxL &
                      + (        inflagU ) * ck(g,k,l,2) * Qin_maxU

             zerosw = 0.5_RP - sign(0.5_RP,abs(Cout)-EPS) ! if Cout = 0, sw = 1

             Qout_min_k = ( ( q(g,k,l) - qnext_max ) + qnext_max*(Cin+Cout-d(g,k,l)) - CQin_max ) &
                        / ( Cout + zerosw ) * ( 1.0_RP - zerosw )                                 &
                        + q(g,k,l) * zerosw
             Qout_max_k = ( ( q(g,k,l) - qnext_min ) + qnext_min*(Cin+Cout-d(g,k,l)) - CQin_min ) &
                        / ( Cout + zerosw ) * ( 1.0_RP - zerosw )                                 &
                        + q(g,k,l) * zerosw

             q_h(g,k,l) = (        inflagL ) * max( min( q_h(g,k,l), Qout_max_km1(g) ), Qout_min_km1(g) ) &
                        + ( 1.0_RP-inflagL ) * max( min( q_h(g,k,l), Qout_max_k      ), Qout_min_k      )

             Qout_min_km1(g) = Qout_min_k
             Qout_max_km1(g) = Qout_max_k
          enddo
          !$omp end do
       enddo

       !$omp end parallel
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl

          do k = ADM_kmin, ADM_kmax
          do g = 1, ADM_gall_pl
             inflagL = 0.5_RP - sign(0.5_RP,ck_pl(g,k  ,l,1)) ! incoming flux: flag=1
             inflagU = 0.5_RP + sign(0.5_RP,ck_pl(g,k+1,l,1)) ! incoming flux: flag=1

             Qin_minL = min( q_pl(g,k,l), q_pl(g,k-1,l) ) + ( 1.0_RP-inflagL ) * BIG
             Qin_minU = min( q_pl(g,k,l), q_pl(g,k+1,l) ) + ( 1.0_RP-inflagU ) * BIG
             Qin_maxL = max( q_pl(g,k,l), q_pl(g,k-1,l) ) - ( 1.0_RP-inflagL ) * BIG
             Qin_maxU = max( q_pl(g,k,l), q_pl(g,k+1,l) ) - ( 1.0_RP-inflagU ) * BIG

             qnext_min = min( Qin_minL, Qin_minU, q_pl(g,k,l) )
             qnext_max = max( Qin_maxL, Qin_maxU, q_pl(g,k,l) )

             Cin      = (        inflagL ) * ( ck_pl(g,k,l,1) ) &
                      + (        inflagU ) * ( ck_pl(g,k,l,2) )
             Cout     = ( 1.0_RP-inflagL ) * ( ck_pl(g,k,l,1) ) &
                      + ( 1.0_RP-inflagU ) * ( ck_pl(g,k,l,2) )

             CQin_max = (        inflagL ) * ( ck_pl(g,k,l,1) * Qin_maxL ) &
                      + (        inflagU ) * ( ck_pl(g,k,l,2) * Qin_maxU )
             CQin_min = (        inflagL ) * ( ck_pl(g,k,l,1) * Qin_minL ) &
                      + (        inflagU ) * ( ck_pl(g,k,l,2) * Qin_minU )

             zerosw = 0.5_RP - sign(0.5_RP,abs(Cout)-EPS) ! if Cout = 0, sw = 1

             Qout_min_pl(g,k) = ( ( q_pl(g,k,l) - qnext_max ) + qnext_max*(Cin+Cout-d_pl(g,k,l)) - CQin_max ) &
                              / ( Cout + zerosw ) * ( 1.0_RP - zerosw )                                       &
                              + q_pl(g,k,l) * zerosw
             Qout_max_pl(g,k) = ( ( q_pl(g,k,l) - qnext_min ) + qnext_min*(Cin+Cout-d_pl(g,k,l)) - CQin_min ) &
                              / ( Cout + zerosw ) * ( 1.0_RP - zerosw )                                       &
                              + q_pl(g,k,l) * zerosw
          enddo
          enddo

          do k = ADM_kmin+1, ADM_kmax
          do g = 1, ADM_gall_pl
             inflagL = 0.5_RP - sign(0.5_RP,ck_pl(g,k,l,1)) ! incoming flux: flag=1

             q_h_pl(g,k,l) = (        inflagL ) * max( min( q_h_pl(g,k,l), Qout_max_pl(g,k-1) ), Qout_min_pl(g,k-1) ) &
                           + ( 1.0_RP-inflagL ) * max( min( q_h_pl(g,k,l), Qout_max_pl(g,k  ) ), Qout_min_pl(g,k  ) )
          enddo
          enddo

       enddo
    endif

    call DEBUG_rapend  ('____vertical_adv_limiter')

    return
  end subroutine vertical_limiter_thuburn

end module mod_src_tracer
