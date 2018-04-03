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
  public :: horizontal_limiter_thuburn

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
  !> Miura(2007)'s scheme with Thuburn(1996) limiter
  subroutine horizontal_limiter_thuburn( &
       q_a,    q_a_pl,  &
       q,      q_pl,    &
       d,      d_pl,    &
       ch,     ch_pl,   &
       cmask,  cmask_pl, &
       Qout_prev, Qout_prev_pl, & ! KERNEL
       Qout_post, Qout_post_pl  ) ! KERNEL
!ESC!    use mod_const, only: &
!ESC!       CONST_HUGE, &
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
!ESC!    use mod_comm, only: &
!ESC!       COMM_data_transfer
!ESC!    implicit none

    real(RP), intent(inout) :: q_a     (ADM_gall   ,ADM_kall,ADM_lall   ,6)
    real(RP), intent(inout) :: q_a_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl  )
    real(RP), intent(in)    :: q       (ADM_gall   ,ADM_kall,ADM_lall     )
    real(RP), intent(in)    :: q_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl  )
    real(RP), intent(in)    :: d       (ADM_gall   ,ADM_kall,ADM_lall     )
    real(RP), intent(in)    :: d_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl  )
    real(RP), intent(in)    :: ch      (ADM_gall   ,ADM_kall,ADM_lall   ,6)
    real(RP), intent(in)    :: ch_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl  )
    real(RP), intent(in)    :: cmask   (ADM_gall   ,ADM_kall,ADM_lall   ,6)
    real(RP), intent(in)    :: cmask_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl  )
    real(RP), intent(out)   :: Qout_prev   (ADM_gall   ,ADM_kall,ADM_lall   ,2 ) ! before communication (for check)
    real(RP), intent(out)   :: Qout_prev_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,2 ) !
    real(RP), intent(in)    :: Qout_post   (ADM_gall   ,ADM_kall,ADM_lall   ,2 ) ! after communication (additional input)
    real(RP), intent(in)    :: Qout_post_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,2 ) !

    real(RP) :: q_min_AI, q_min_AIJ, q_min_AJ, q_min_pl
    real(RP) :: q_max_AI, q_max_AIJ, q_max_AJ, q_max_pl

    real(RP) :: qnext_min   , qnext_min_pl
    real(RP) :: qnext_max   , qnext_max_pl
    real(RP) :: Cin_sum     , Cin_sum_pl
    real(RP) :: Cout_sum    , Cout_sum_pl
    real(RP) :: CQin_max_sum, CQin_max_sum_pl
    real(RP) :: CQin_min_sum, CQin_min_sum_pl

    integer, parameter :: I_min = 1
    integer, parameter :: I_max = 2
    real(RP) :: Qin    (ADM_gall   ,ADM_kall,ADM_lall   ,2,6)
    real(RP) :: Qin_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,2,2)
    real(RP) :: Qout   (ADM_gall   ,ADM_kall,ADM_lall   ,2  )
    real(RP) :: Qout_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,2  )

    real(RP) :: ch_masked1
    real(RP) :: ch_masked2
    real(RP) :: ch_masked3
    real(RP) :: ch_masked4
    real(RP) :: ch_masked5
    real(RP) :: ch_masked6
    real(RP) :: ch_masked
    real(RP) :: zerosw

    integer  :: gmin, gmax, kall, iall
    real(RP) :: EPS, BIG

    integer  :: ij
    integer  :: ip1j, ijp1, ip1jp1, ip2jp1
    integer  :: im1j, ijm1

    integer  :: i, j, k, l, n, v
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('____horizontal_adv_limiter')

    gmin = ADM_gmin
    gmax = ADM_gmax
    kall = ADM_kall
    iall = ADM_gall_1d

    EPS  = CONST_EPS
    BIG  = CONST_HUGE

    do l = 1, ADM_lall
       !$omp parallel default(none), &
       !$omp private(i,j,k,ij,ip1j,ip1jp1,ijp1,im1j,ijm1,ip2jp1,                        &
       !$omp         q_min_AI,q_min_AIJ,q_min_AJ,q_max_AI,q_max_AIJ,q_max_AJ,zerosw,    &
       !$omp         ch_masked1,ch_masked2,ch_masked3,ch_masked4,ch_masked5,ch_masked6, &
       !$omp         qnext_min,qnext_max,Cin_sum,Cout_sum,CQin_min_sum,CQin_max_sum),   &
       !$omp shared(l,ADM_have_sgp,gmin,gmax,kall,iall,q,cmask,d,ch,Qin,Qout,EPS,BIG)
       do k = 1, kall
          !---< (i) define inflow bounds, eq.(32)&(33) >---
!OCL XFILL
          !$omp do
          do j = gmin-1, gmax
          do i = gmin-1, gmax
             ij     = (j-1)*iall + i
             ip1j   = ij + 1
             ip1jp1 = ij + iall + 1
             ijp1   = ij + iall
             im1j   = ij - 1
             ijm1   = ij - iall

             im1j   = max( im1j  , 1 )
             ijm1   = max( ijm1  , 1 )

             q_min_AI  = min( q(ij,k,l), q(ijm1,k,l), q(ip1j,k,l), q(ip1jp1,k,l) )
             q_max_AI  = max( q(ij,k,l), q(ijm1,k,l), q(ip1j,k,l), q(ip1jp1,k,l) )
             q_min_AIJ = min( q(ij,k,l), q(ip1j,k,l), q(ip1jp1,k,l), q(ijp1,k,l) )
             q_max_AIJ = max( q(ij,k,l), q(ip1j,k,l), q(ip1jp1,k,l), q(ijp1,k,l) )
             q_min_AJ  = min( q(ij,k,l), q(ip1jp1,k,l), q(ijp1,k,l), q(im1j,k,l) )
             q_max_AJ  = max( q(ij,k,l), q(ip1jp1,k,l), q(ijp1,k,l), q(im1j,k,l) )

             Qin(ij,    k,l,I_min,1) = (        cmask(ij,k,l,1) ) * q_min_AI &
                                     + ( 1.0_RP-cmask(ij,k,l,1) ) * BIG
             Qin(ip1j,  k,l,I_min,4) = (        cmask(ij,k,l,1) ) * BIG      &
                                     + ( 1.0_RP-cmask(ij,k,l,1) ) * q_min_AI
             Qin(ij,    k,l,I_max,1) = (        cmask(ij,k,l,1) ) * q_max_AI &
                                     + ( 1.0_RP-cmask(ij,k,l,1) ) * (-BIG)
             Qin(ip1j,  k,l,I_max,4) = (        cmask(ij,k,l,1) ) * (-BIG)   &
                                     + ( 1.0_RP-cmask(ij,k,l,1) ) * q_max_AI

             Qin(ij,    k,l,I_min,2) = (        cmask(ij,k,l,2) ) * q_min_AIJ &
                                     + ( 1.0_RP-cmask(ij,k,l,2) ) * BIG
             Qin(ip1jp1,k,l,I_min,5) = (        cmask(ij,k,l,2) ) * BIG       &
                                     + ( 1.0_RP-cmask(ij,k,l,2) ) * q_min_AIJ
             Qin(ij,    k,l,I_max,2) = (        cmask(ij,k,l,2) ) * q_max_AIJ &
                                     + ( 1.0_RP-cmask(ij,k,l,2) ) * (-BIG)
             Qin(ip1jp1,k,l,I_max,5) = (        cmask(ij,k,l,2) ) * (-BIG)    &
                                     + ( 1.0_RP-cmask(ij,k,l,2) ) * q_max_AIJ

             Qin(ij,    k,l,I_min,3) = (        cmask(ij,k,l,3) ) * q_min_AJ &
                                     + ( 1.0_RP-cmask(ij,k,l,3) ) * BIG
             Qin(ijp1,  k,l,I_min,6) = (        cmask(ij,k,l,3) ) * BIG      &
                                     + ( 1.0_RP-cmask(ij,k,l,3) ) * q_min_AJ
             Qin(ij,    k,l,I_max,3) = (        cmask(ij,k,l,3) ) * q_max_AJ &
                                     + ( 1.0_RP-cmask(ij,k,l,3) ) * (-BIG)
             Qin(ijp1,  k,l,I_max,6) = (        cmask(ij,k,l,3) ) * (-BIG)   &
                                     + ( 1.0_RP-cmask(ij,k,l,3) ) * q_max_AJ
          enddo
          enddo
          !$omp end do

          if ( ADM_have_sgp(l) ) then
             !$omp master
             j = gmin-1
             i = gmin-1

             ij     = (j-1)*iall + i
             ijp1   = ij + iall
             ip1jp1 = ij + iall + 1
             ip2jp1 = ij + iall + 2

             q_min_AIJ = min( q(ij,k,l), q(ip1jp1,k,l), q(ip2jp1,k,l), q(ijp1,k,l) )
             q_max_AIJ = max( q(ij,k,l), q(ip1jp1,k,l), q(ip2jp1,k,l), q(ijp1,k,l) )

             Qin(ij,    k,l,I_min,2) = (        cmask(ij,k,l,2) ) * q_min_AIJ &
                                     + ( 1.0_RP-cmask(ij,k,l,2) ) * BIG
             Qin(ip1jp1,k,l,I_min,5) = (        cmask(ij,k,l,2) ) * BIG       &
                                     + ( 1.0_RP-cmask(ij,k,l,2) ) * q_min_AIJ
             Qin(ij,    k,l,I_max,2) = (        cmask(ij,k,l,2) ) * q_max_AIJ &
                                     + ( 1.0_RP-cmask(ij,k,l,2) ) * (-BIG)
             Qin(ip1jp1,k,l,I_max,5) = (        cmask(ij,k,l,2) ) * (-BIG)    &
                                     + ( 1.0_RP-cmask(ij,k,l,2) ) * q_max_AIJ
             !$omp end master
          endif

          !---< (iii) define allowable range of q at next step, eq.(42)&(43) >---
!OCL XFILL
          !$omp do
          do j = gmin, gmax
          do i = gmin, gmax
             ij = (j-1)*iall + i

             qnext_min = min( q(ij,k,l),           &
                              Qin(ij,k,l,I_min,1), &
                              Qin(ij,k,l,I_min,2), &
                              Qin(ij,k,l,I_min,3), &
                              Qin(ij,k,l,I_min,4), &
                              Qin(ij,k,l,I_min,5), &
                              Qin(ij,k,l,I_min,6)  )

             qnext_max = max( q(ij,k,l),           &
                              Qin(ij,k,l,I_max,1), &
                              Qin(ij,k,l,I_max,2), &
                              Qin(ij,k,l,I_max,3), &
                              Qin(ij,k,l,I_max,4), &
                              Qin(ij,k,l,I_max,5), &
                              Qin(ij,k,l,I_max,6)  )

             ch_masked1 = min( ch(ij,k,l,1), 0.0_RP )
             ch_masked2 = min( ch(ij,k,l,2), 0.0_RP )
             ch_masked3 = min( ch(ij,k,l,3), 0.0_RP )
             ch_masked4 = min( ch(ij,k,l,4), 0.0_RP )
             ch_masked5 = min( ch(ij,k,l,5), 0.0_RP )
             ch_masked6 = min( ch(ij,k,l,6), 0.0_RP )

             Cin_sum      = ch_masked1 &
                          + ch_masked2 &
                          + ch_masked3 &
                          + ch_masked4 &
                          + ch_masked5 &
                          + ch_masked6

             Cout_sum     = ch(ij,k,l,1) - ch_masked1 &
                          + ch(ij,k,l,2) - ch_masked2 &
                          + ch(ij,k,l,3) - ch_masked3 &
                          + ch(ij,k,l,4) - ch_masked4 &
                          + ch(ij,k,l,5) - ch_masked5 &
                          + ch(ij,k,l,6) - ch_masked6

             CQin_min_sum = ch_masked1 * Qin(ij,k,l,I_min,1) &
                          + ch_masked2 * Qin(ij,k,l,I_min,2) &
                          + ch_masked3 * Qin(ij,k,l,I_min,3) &
                          + ch_masked4 * Qin(ij,k,l,I_min,4) &
                          + ch_masked5 * Qin(ij,k,l,I_min,5) &
                          + ch_masked6 * Qin(ij,k,l,I_min,6)

             CQin_max_sum = ch_masked1 * Qin(ij,k,l,I_max,1) &
                          + ch_masked2 * Qin(ij,k,l,I_max,2) &
                          + ch_masked3 * Qin(ij,k,l,I_max,3) &
                          + ch_masked4 * Qin(ij,k,l,I_max,4) &
                          + ch_masked5 * Qin(ij,k,l,I_max,5) &
                          + ch_masked6 * Qin(ij,k,l,I_max,6)

             zerosw = 0.5_RP - sign(0.5_RP,abs(Cout_sum)-EPS) ! if Cout_sum = 0, sw = 1

             Qout(ij,k,l,I_min) = ( q(ij,k,l) - CQin_max_sum - qnext_max*(1.0_RP-Cin_sum-Cout_sum+d(ij,k,l)) ) &
                                / ( Cout_sum + zerosw ) * ( 1.0_RP - zerosw )                                         &
                                + q(ij,k,l) * zerosw
             Qout(ij,k,l,I_max) = ( q(ij,k,l) - CQin_min_sum - qnext_min*(1.0_RP-Cin_sum-Cout_sum+d(ij,k,l)) ) &
                                / ( Cout_sum + zerosw ) * ( 1.0_RP - zerosw )                                         &
                                + q(ij,k,l) * zerosw
          enddo
          enddo
          !$omp end do

!OCL XFILL
          !$omp do
          do j = 1, iall
          do i = 1, iall
             if (      i < gmin .OR. i > gmax &
                  .OR. j < gmin .OR. j > gmax ) then
                ij = (j-1)*iall + i

                Qout(ij,k,l,I_min) = q(ij,k,l)
                Qout(ij,k,l,I_min) = q(ij,k,l)
                Qout(ij,k,l,I_max) = q(ij,k,l)
                Qout(ij,k,l,I_max) = q(ij,k,l)
             endif
          enddo
          enddo
          !$omp end do

       enddo ! k loop
       !$omp end parallel
    enddo ! l loop

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijp1 = v + 1
                ijm1 = v - 1
                if( ijp1 == ADM_gmax_pl+1 ) ijp1 = ADM_gmin_pl
                if( ijm1 == ADM_gmin_pl-1 ) ijm1 = ADM_gmax_pl

                q_min_pl = min( q_pl(n,k,l), q_pl(ij,k,l), q_pl(ijm1,k,l), q_pl(ijp1,k,l) )
                q_max_pl = max( q_pl(n,k,l), q_pl(ij,k,l), q_pl(ijm1,k,l), q_pl(ijp1,k,l) )

                Qin_pl(ij,k,l,I_min,1) = (        cmask_pl(ij,k,l) ) * q_min_pl &
                                       + ( 1.0_RP-cmask_pl(ij,k,l) ) * BIG
                Qin_pl(ij,k,l,I_min,2) = (        cmask_pl(ij,k,l) ) * BIG      &
                                       + ( 1.0_RP-cmask_pl(ij,k,l) ) * q_min_pl
                Qin_pl(ij,k,l,I_max,1) = (        cmask_pl(ij,k,l) ) * q_max_pl &
                                       + ( 1.0_RP-cmask_pl(ij,k,l) ) * (-BIG)
                Qin_pl(ij,k,l,I_max,2) = (        cmask_pl(ij,k,l) ) * (-BIG)   &
                                       + ( 1.0_RP-cmask_pl(ij,k,l) ) * q_max_pl
             enddo

             qnext_min_pl = q_pl(n,k,l)
             qnext_max_pl = q_pl(n,k,l)
             do v = ADM_gmin_pl, ADM_gmax_pl
                qnext_min_pl = min( qnext_min_pl, Qin_pl(v,k,l,I_min,1) )
                qnext_max_pl = max( qnext_max_pl, Qin_pl(v,k,l,I_max,1) )
             enddo

             Cin_sum_pl      = 0.0_RP
             Cout_sum_pl     = 0.0_RP
             CQin_max_sum_pl = 0.0_RP
             CQin_min_sum_pl = 0.0_RP
             do v = ADM_gmin_pl, ADM_gmax_pl
                ch_masked = cmask_pl(v,k,l) * ch_pl(v,k,l)

                Cin_sum_pl      = Cin_sum_pl      + ch_masked
                Cout_sum_pl     = Cout_sum_pl     - ch_masked + ch_pl (v,k,l)
                CQin_min_sum_pl = CQin_min_sum_pl + ch_masked * Qin_pl(v,k,l,I_min,1)
                CQin_max_sum_pl = CQin_max_sum_pl + ch_masked * Qin_pl(v,k,l,I_max,1)
             enddo

             zerosw = 0.5_RP - sign(0.5_RP,abs(Cout_sum_pl)-EPS) ! if Cout_sum_pl = 0, sw = 1

             Qout_pl(n,k,l,I_min) = ( q_pl(n,k,l) - CQin_max_sum_pl - qnext_max_pl*(1.0_RP-Cin_sum_pl-Cout_sum_pl+d_pl(n,k,l)) ) &
                                  / ( Cout_sum_pl + zerosw ) * ( 1.0_RP - zerosw )                                               &
                                  + q_pl(n,k,l) * zerosw
             Qout_pl(n,k,l,I_max) = ( q_pl(n,k,l) - CQin_min_sum_pl - qnext_min_pl*(1.0_RP-Cin_sum_pl-Cout_sum_pl+d_pl(n,k,l)) ) &
                                  / ( Cout_sum_pl + zerosw ) * ( 1.0_RP - zerosw )                                               &
                                  + q_pl(n,k,l) * zerosw

          enddo
       enddo
    endif

    !#####################################################################KERNEL
    call DEBUG_rapend  ('____horizontal_adv_limiter')
    Qout_pl(ADM_gmin_pl:ADM_gmax_pl,:,:,:) = 0.0_RP

    Qout_prev   (:,:,:,:) = Qout        (:,:,:,:)
    Qout_prev_pl(:,:,:,:) = Qout_pl     (:,:,:,:)
    !call COMM_data_transfer( Qout(:,:,:,:), Qout_pl(:,:,:,:) )
    Qout        (:,:,:,:) = Qout_post   (:,:,:,:)
    Qout_pl     (:,:,:,:) = Qout_post_pl(:,:,:,:)
    call DEBUG_rapstart('____horizontal_adv_limiter')
    !#####################################################################KERNEL

    !---- apply inflow/outflow limiter
    do l = 1, ADM_lall
       !$omp parallel do default(none),private(i,j,k,ij,ip1j,ip1jp1,ijp1), &
       !$omp shared(l,gmin,gmax,kall,iall,q_a,cmask,Qin,Qout)
       do k = 1, kall
          do j = gmin-1, gmax
          do i = gmin-1, gmax
             ij     = (j-1)*iall + i
             ip1j   = ij + 1
             ip1jp1 = ij + iall + 1
             ijp1   = ij + iall

             q_a(ij,k,l,1) = (        cmask(ij,k,l,1) ) * min( max( q_a(ij,k,l,1), Qin (ij    ,k,l,I_min,1) ), Qin (ij    ,k,l,I_max,1) ) &
                           + ( 1.0_RP-cmask(ij,k,l,1) ) * min( max( q_a(ij,k,l,1), Qin (ip1j  ,k,l,I_min,4) ), Qin (ip1j  ,k,l,I_max,4) )
             q_a(ij,k,l,1) = (        cmask(ij,k,l,1) ) * max( min( q_a(ij,k,l,1), Qout(ip1j  ,k,l,I_max  ) ), Qout(ip1j  ,k,l,I_min  ) ) &
                           + ( 1.0_RP-cmask(ij,k,l,1) ) * max( min( q_a(ij,k,l,1), Qout(ij    ,k,l,I_max  ) ), Qout(ij    ,k,l,I_min  ) )
             q_a(ip1j,k,l,4) = q_a(ij,k,l,1)

             q_a(ij,k,l,2) = (        cmask(ij,k,l,2) ) * min( max( q_a(ij,k,l,2), Qin (ij    ,k,l,I_min,2) ), Qin (ij    ,k,l,I_max,2) ) &
                           + ( 1.0_RP-cmask(ij,k,l,2) ) * min( max( q_a(ij,k,l,2), Qin (ip1jp1,k,l,I_min,5) ), Qin (ip1jp1,k,l,I_max,5) )
             q_a(ij,k,l,2) = (        cmask(ij,k,l,2) ) * max( min( q_a(ij,k,l,2), Qout(ip1jp1,k,l,I_max  ) ), Qout(ip1jp1,k,l,I_min  ) ) &
                           + ( 1.0_RP-cmask(ij,k,l,2) ) * max( min( q_a(ij,k,l,2), Qout(ij    ,k,l,I_max  ) ), Qout(ij    ,k,l,I_min  ) )
             q_a(ip1jp1,k,l,5) = q_a(ij,k,l,2)

             q_a(ij,k,l,3) = (        cmask(ij,k,l,3) ) * min( max( q_a(ij,k,l,3), Qin (ij    ,k,l,I_min,3) ), Qin (ij    ,k,l,I_max,3) ) &
                           + ( 1.0_RP-cmask(ij,k,l,3) ) * min( max( q_a(ij,k,l,3), Qin (ijp1  ,k,l,I_min,6) ), Qin (ijp1  ,k,l,I_max,6) )
             q_a(ij,k,l,3) = (        cmask(ij,k,l,3) ) * max( min( q_a(ij,k,l,3), Qout(ijp1  ,k,l,I_max  ) ), Qout(ijp1  ,k,l,I_min  ) ) &
                           + ( 1.0_RP-cmask(ij,k,l,3) ) * max( min( q_a(ij,k,l,3), Qout(ij    ,k,l,I_max  ) ), Qout(ij    ,k,l,I_min  ) )
             q_a(ijp1,k,l,6) = q_a(ij,k,l,3)
          enddo
          enddo
       enddo
       !$omp end parallel do
    enddo

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
       do v = ADM_gmin_pl, ADM_gmax_pl
          q_a_pl(v,k,l) = (        cmask_pl(v,k,l) ) * min(max(q_a_pl(v,k,l), Qin_pl (v,k,l,I_min,1)), Qin_pl (v,k,l,I_max,1)) &
                        + ( 1.0_RP-cmask_pl(v,k,l) ) * min(max(q_a_pl(v,k,l), Qin_pl (v,k,l,I_min,2)), Qin_pl (v,k,l,I_max,2))
          q_a_pl(v,k,l) = (        cmask_pl(v,k,l) ) * max(min(q_a_pl(v,k,l), Qout_pl(v,k,l,I_max  )), Qout_pl(v,k,l,I_min  )) &
                        + ( 1.0_RP-cmask_pl(v,k,l) ) * max(min(q_a_pl(v,k,l), Qout_pl(n,k,l,I_max  )), Qout_pl(n,k,l,I_min  ))
       enddo
       enddo
       enddo
    endif

    call DEBUG_rapend  ('____horizontal_adv_limiter')

    return
  end subroutine horizontal_limiter_thuburn

end module mod_src_tracer
