!-------------------------------------------------------------------------------
!
!+  
!
!-------------------------------------------------------------------------------
program comp_pvort
  use prec
  use mod_misc
  use caldyn_gcm_mod, only : compute_pvort
  implicit none


  real(rstd), allocatable :: ORG_u           (:,:)
  real(rstd), allocatable :: ORG_ps          (:)  
  real(rstd), allocatable :: ORG_theta_rhodz (:,:)
  real(rstd), allocatable :: ORG_theta       (:,:)
  real(rstd), allocatable :: ORG_rhodz       (:,:)
  real(rstd), allocatable :: ORG_qu          (:,:)
  real(rstd), allocatable :: ORG_qv          (:,:)

  real(rstd), allocatable :: ORG_theta_prev  (:,:)
  real(rstd), allocatable :: ORG_rhodz_prev  (:,:)
  real(rstd), allocatable :: ORG_qu_prev     (:,:)
  real(rstd), allocatable :: ORG_qv_prev     (:,:)

  real(rstd), allocatable :: u           (:,:)
  real(rstd), allocatable :: ps          (:)  
  real(rstd), allocatable :: theta_rhodz (:,:)
  real(rstd), allocatable :: theta       (:,:)
  real(rstd), allocatable :: rhodz       (:,:)
  real(rstd), allocatable :: qu          (:,:)
  real(rstd), allocatable :: qv          (:,:)

  real(rstd), allocatable :: theta_prev  (:,:)
  real(rstd), allocatable :: rhodz_prev  (:,:)
  real(rstd), allocatable :: qu_prev     (:,:)
  real(rstd), allocatable :: qv_prev     (:,:)


  integer :: iteration
!!$       indices(:) = (/&
!!$            & iim,jjm,llm,&
!!$            & ij_begin,ij_end,&
!!$            & ij_begin_ext,ij_end_ext,&
!!$            & ll_begin,ll_end,&
!!$            & t_right, t_rup, t_lup, t_left, t_ldown, t_rdown,&
!!$            & u_right, u_rup, u_lup, u_left, u_ldown, u_rdown, &
!!$            & caldyn_eta
!!$            & /)
  integer,parameter :: indices_size = 28
  integer(i8) :: indices(indices_size)


  !=============================================================================

  write(*,*) "[KERNEL] comp_pvort"

  !###############################################################################

  write(*,*) "*** Start  initialize"



  !###############################################################################

  !---< read input data >---

  call dumpio_syscheck
  call dumpio_mk_fname(EX_fname,'snapshot.comp_pvort','pe',SET_prc_me-1,6)
  call dumpio_fopen(EX_fid,EX_fname,IO_FREAD)

  call dumpio_read_data( EX_fid, size(indices), indices )

  iim = indices(1)
  jjm = indices(2)
  llm = indices(3)
!!$   write(*,'(A30,3I6)') 'dbg: iim,jjm,llm:',iim,jjm,llm

  allocate( ORG_u           (iim*3*jjm,llm )     )
  allocate( ORG_ps          (iim*jjm       )     )
  allocate( ORG_theta_rhodz (iim*jjm,llm   )     )
  allocate( ORG_theta       (iim*jjm,llm   )     )
  allocate( ORG_rhodz       (iim*jjm,llm   )     )
  allocate( ORG_qu          (iim*3*jjm,llm )     )
  allocate( ORG_qv          (iim*2*jjm,llm )     )

  allocate( ORG_theta_prev  (iim*jjm,llm   )     )
  allocate( ORG_rhodz_prev  (iim*jjm,llm   )     )
  allocate( ORG_qu_prev     (iim*3*jjm,llm )     )
  allocate( ORG_qv_prev     (iim*2*jjm,llm )     )

  allocate( Av          (iim*2*jjm) ) !! field_type is field_Z
  allocate( de          (iim*3*jjm) ) !! field_type is field_U
  allocate( Riv2        (iim*jjm,6) ) !! field_type is field_t,dim1=6
  allocate( fv          (iim*2*jjm) ) !! field_type is field_Z

  allocate( mass_dak(llm) )
  allocate( mass_dbk(llm) ) 

  call dumpio_read_data( EX_fid, 1 , g )

  call dumpio_read_data( EX_fid, iim*3*jjm*llm , ORG_u           (:,:)     )
  call dumpio_read_data( EX_fid, iim*jjm       , ORG_ps          (:)       )
  call dumpio_read_data( EX_fid, iim*jjm*llm   , ORG_theta_rhodz (:,:)     )
  call dumpio_read_data( EX_fid, iim*jjm*llm   , ORG_theta       (:,:)     )
  call dumpio_read_data( EX_fid, iim*jjm*llm   , ORG_rhodz       (:,:)     )
  call dumpio_read_data( EX_fid, iim*3*jjm*llm , ORG_qu          (:,:)     )
  call dumpio_read_data( EX_fid, iim*2*jjm*llm , ORG_qv          (:,:)     )

  call dumpio_read_data( EX_fid, iim*jjm*llm   , ORG_theta_prev  (:,:)     )
  call dumpio_read_data( EX_fid, iim*jjm*llm   , ORG_rhodz_prev  (:,:)     )
  call dumpio_read_data( EX_fid, iim*3*jjm*llm , ORG_qu_prev     (:,:)     )
  call dumpio_read_data( EX_fid, iim*2*jjm*llm , ORG_qv_prev     (:,:)     )

  call dumpio_read_data( EX_fid, iim*2*jjm , Av(:) )
  call dumpio_read_data( EX_fid, iim*3*jjm , de(:) )
  call dumpio_read_data( EX_fid, iim*jjm*6 , Riv2(:,:) )
  call dumpio_read_data( EX_fid, iim*2*jjm , fv(:) )

  call dumpio_read_data( EX_fid, llm , mass_dak(:) )
  call dumpio_read_data( EX_fid, llm , mass_dbk(:) )

  call dumpio_fclose(EX_fid)


  !! I know these are implicit cast.
  ij_begin  = indices(4)
  ij_end    = indices(5)
  ij_begin_ext  = indices(6)
  ij_end_ext    = indices(7)
  ll_begin  = indices(8)
  ll_end    = indices(9)

  t_right = indices(10)
  t_rup   = indices(11)
  t_lup   = indices(12)
  t_left  = indices(13)
  t_ldown = indices(14)
  t_rdown = indices(15)

  u_right = indices(16)
  u_rup   = indices(17)
  u_lup   = indices(18)
  u_left  = indices(19)
  u_ldown = indices(20)
  u_rdown = indices(21)

  z_rup   = indices(22)
  z_up    = indices(23)
  z_lup   = indices(24)
  z_ldown = indices(25)
  z_down  = indices(26)
  z_rdown = indices(27)

  caldyn_eta = indices(28)

  allocate( u           (iim*3*jjm,llm )     )
  allocate( ps          (iim*jjm       )     )
  allocate( theta_rhodz (iim*jjm,llm   )     )
  allocate( theta       (iim*jjm,llm   )     )
  allocate( rhodz       (iim*jjm,llm   )     )
  allocate( qu          (iim*3*jjm,llm )     )
  allocate( qv          (iim*2*jjm,llm )     )

  allocate( theta_prev  (iim*jjm,llm   )     )
  allocate( rhodz_prev  (iim*jjm,llm   )     )
  allocate( qu_prev     (iim*3*jjm,llm )     )
  allocate( qv_prev     (iim*2*jjm,llm )     )


  u           (:,:) =   ORG_u           (:,:)
  ps          (:)   =   ORG_ps          (:)  
  theta_rhodz (:,:) =   ORG_theta_rhodz (:,:)
!!$   theta       (:,:) =   ORG_theta       (:,:)
!!$   rhodz       (:,:) =   ORG_rhodz       (:,:)
!!$   qu          (:,:) =   ORG_qu          (:,:)
!!$   qv          (:,:) =   ORG_qv          (:,:)

  theta_prev  (:,:) =   ORG_theta_prev  (:,:)
  rhodz_prev  (:,:) =   ORG_rhodz_prev  (:,:)
  qu_prev     (:,:) =   ORG_qu_prev     (:,:)
  qv_prev     (:,:) =   ORG_qv_prev     (:,:)

  write(*,'(A30,3I6)') 'iim, jjm, llm:',iim,jjm,llm
  write(*,'(A30,2I6)') 'ij_begin, ij_end:',ij_begin,ij_end
  write(*,'(A30,2I6)') 'ij_begin_ext, ij_end_ext:',ij_begin_ext,ij_end_ext
  write(*,'(A30,2I6)') 'll_begin, ll_end:',ll_begin,ll_end

  write(*,'(A30,3I6)') 't_right, t_rup, t_lup:',t_right,t_rup,t_lup
  write(*,'(A30,3I6)') 't_left, t_ldown, t_rdown:',t_left,t_ldown,t_rdown

  write(*,'(A30,3I6)') 'u_right, u_rup, u_lup:',u_right,u_rup,u_lup
  write(*,'(A30,3I6)') 'u_left, u_ldown, u_rdown:',u_left,u_ldown,u_rdown

  write(*,'(A30,3I6)') 'z_rup, z_up, z_lup:',z_rup,z_up,z_lup
  write(*,'(A30,3I6)') 'z_ldown, z_down, z_rdown:',z_ldown,z_down,z_rdown

  write(*,'(A30, I6)') 'caldyn_eta:',caldyn_eta
  write(*,'(A30,F15.8)') 'g:',g


  call DEBUG_valuecheck(  'Av   ', Av   (:) )
  call DEBUG_valuecheck(  'de   ', de   (:) )
  call DEBUG_valuecheck(  'Riv2 ', Riv2 (:,:) )
  call DEBUG_valuecheck(  'fv   ', fv   (:) )
  call DEBUG_valuecheck(  'mass_dak   ', mass_dak   (:) )
  call DEBUG_valuecheck(  'mass_dbk   ', mass_dbk   (:) )

  write(*,*) "*** Finish initialize"

  !###############################################################################

  write(*,*) "*** Start kernel"

  do iteration = 1, SET_iteration

    theta(:,:) = theta_prev(:,:)
    rhodz(:,:) = rhodz_prev(:,:)
    qu   (:,:) = qu_prev   (:,:)
    qv   (:,:) = qv_prev   (:,:)

    call DEBUG_rapstart('MAIN_comp_pvort')

    call compute_pvort( & ! [IN]
         & ps,          & ! [IN]
         & u,           & ! [IN]
         & theta_rhodz, & ! [IN]
         & rhodz,       & ! [INOUT]
         & theta,       & ! [INOUT]
         & qu,          & ! [INOUT]
         & qv           & ! [INOUT]
         )


    call DEBUG_rapend('MAIN_comp_pvort')

    if ( SET_iteration<5 .or. mod(iteration, SET_iteration/1)==0 ) then
      write(ADM_LOG_FID,*) '### check point iteration:',iteration
      write(ADM_LOG_FID,*) '### Input ###'

      call DEBUG_valuecheck(  'u          ', u           (:,:) )
      call DEBUG_valuecheck(  'ps         ', ps          (:) )
      call DEBUG_valuecheck(  'theta_rhodz', theta_rhodz (:,:) )
      call DEBUG_valuecheck(  'theta_prev ', theta_prev  (:,:) )
      call DEBUG_valuecheck(  'rhodz_prev ', rhodz_prev  (:,:) )
      call DEBUG_valuecheck(  'qu_prev    ', qu_prev     (:,:) )
      call DEBUG_valuecheck(  'qv_prev    ', qv_prev     (:,:) )


      write(ADM_LOG_FID,*) '### Output ###'
      call DEBUG_valuecheck(  'theta      ', theta       (:,:) )
      call DEBUG_valuecheck(  'rhodz      ', rhodz       (:,:) )
      call DEBUG_valuecheck(  'qu         ', qu          (:,:) )
      call DEBUG_valuecheck(  'qv         ', qv          (:,:) )
    end if
  end do

  write(ADM_LOG_FID,*) '### final iteration:',iteration-1
  write(ADM_LOG_FID,*) '### Validation : grid-by-grid diff ###'

  theta(:,:) = ORG_theta(:,:) - theta(:,:)
  rhodz(:,:) = ORG_rhodz(:,:) - rhodz(:,:)
  qu   (:,:) = ORG_qu   (:,:) - qu   (:,:) 
  qv   (:,:) = ORG_qv   (:,:) - qv   (:,:) 

  call DEBUG_valuecheck( 'theta', theta(:,:) )
  call DEBUG_valuecheck( 'rhodz', rhodz(:,:) )
  call DEBUG_valuecheck( 'qu'   , qu   (:,:) )
  call DEBUG_valuecheck( 'qv'   , qv   (:,:) )

  write(*,*) "*** Finish kernel"

  !###############################################################################

  call DEBUG_rapreport

end program comp_pvort




!!$ Local Variables:
!!$ mode: f90
!!$ f90-do-indent: 2
!!$ f90-if-indent: 2
!!$ f90-program-indent: 2
!!$ f90-type-indent: 2
!!$ f90-directive-comment-re: "!\$omp "
!!$ end:

