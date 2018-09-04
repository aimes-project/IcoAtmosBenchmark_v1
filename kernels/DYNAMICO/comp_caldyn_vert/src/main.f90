!-------------------------------------------------------------------------------
!
!+  
!
!-------------------------------------------------------------------------------
program comp_caldyn_vert
  use prec
  use mod_misc
  use caldyn_gcm_mod, only : compute_caldyn_vert
  implicit none


  REAL(rstd),ALLOCATABLE :: ORG_u      (:,:) ! prognostic "velocity"
  REAL(rstd),ALLOCATABLE :: ORG_theta  (:,:) ! potential temperature
  REAL(rstd),ALLOCATABLE :: ORG_rhodz  (:,:)

  REAL(rstd),ALLOCATABLE :: ORG_convm        (:,:) ! mass flux convergence
  REAL(rstd),ALLOCATABLE :: ORG_wflux        (:,:) ! vertical mass flux (kg/m2/s)
  REAL(rstd),allocatable :: ORG_wwuu         (:,:) !(iim*3*jjm,llm+1)
  REAL(rstd),ALLOCATABLE :: ORG_du           (:,:)
  REAL(rstd),ALLOCATABLE :: ORG_dtheta_rhodz (:,:)
  REAL(rstd),allocatable :: ORG_dps          (:)   !(iim*jjm)

  REAL(rstd),ALLOCATABLE :: ORG_convm_prev        (:,:) ! mass flux convergence
  REAL(rstd),ALLOCATABLE :: ORG_wflux_prev        (:,:) ! vertical mass flux (kg/m2/s)
  REAL(rstd),allocatable :: ORG_wwuu_prev         (:,:) !(iim*3*jjm,llm+1)
  REAL(rstd),ALLOCATABLE :: ORG_du_prev           (:,:)
  REAL(rstd),ALLOCATABLE :: ORG_dtheta_rhodz_prev (:,:)
  REAL(rstd),allocatable :: ORG_dps_prev          (:)   !(iim*jjm)

  REAL(rstd),ALLOCATABLE :: u      (:,:) ! prognostic "velocity"
  REAL(rstd),ALLOCATABLE :: theta  (:,:) ! potential temperature
  REAL(rstd),ALLOCATABLE :: rhodz  (:,:)

  REAL(rstd),ALLOCATABLE :: convm        (:,:) ! mass flux convergence
  REAL(rstd),ALLOCATABLE :: wflux        (:,:) ! vertical mass flux (kg/m2/s)
  REAL(rstd),allocatable :: wwuu         (:,:) !(iim*3*jjm,llm+1)
  REAL(rstd),ALLOCATABLE :: du           (:,:)
  REAL(rstd),ALLOCATABLE :: dtheta_rhodz (:,:)
  REAL(rstd),allocatable :: dps          (:)   !(iim*jjm)

  REAL(rstd),ALLOCATABLE :: convm_prev        (:,:) ! mass flux convergence
  REAL(rstd),ALLOCATABLE :: wflux_prev        (:,:) ! vertical mass flux (kg/m2/s)
  REAL(rstd),allocatable :: wwuu_prev         (:,:) !(iim*3*jjm,llm+1)
  REAL(rstd),ALLOCATABLE :: du_prev           (:,:)
  REAL(rstd),ALLOCATABLE :: dtheta_rhodz_prev (:,:)
  REAL(rstd),allocatable :: dps_prev          (:)   !(iim*jjm)


  integer :: iteration
!!$       indices(:) = (/&
!!$            & iim,jjm,llm,&
!!$            & ij_begin,ij_end,&
!!$            & ij_begin_ext,ij_end_ext,&
!!$            & ll_begin,ll_end,&
!!$            & t_right, t_rup, t_lup, t_left, t_ldown, t_rdown,&
!!$            & u_right, u_rup, u_lup, u_left, u_ldown, u_rdown, &
!!$            & z_rup, z_up, z_lup, z_ldown, z_down, z_rdown,&
!!$            & ij_omp_begin, ij_omp_end, &
!!$            & ll_beginp1, ll_endm1 &
!!$            & /)
  integer,parameter :: indices_size = 31
  integer(i8) :: indices(indices_size)



  !=============================================================================

  write(*,*) "[KERNEL] comp_caldyn_vert"

  !###############################################################################

  write(*,*) "*** Start  initialize"



  !###############################################################################

  !---< read input data >---

  call dumpio_syscheck
  call dumpio_mk_fname(EX_fname,'snapshot.comp_caldyn_vert','pe',SET_prc_me-1,6)
  call dumpio_fopen(EX_fid,EX_fname,IO_FREAD)

  call dumpio_read_data( EX_fid, size(indices), indices )

  iim = indices(1)
  jjm = indices(2)
  llm = indices(3)
!!$  write(*,'(A30,3I6)') 'dbg: iim,jjm,llm:',iim,jjm,llm

  allocate( ORG_u       (iim*3*jjm, llm)   ) ! prognostic "velocity"
  allocate( ORG_theta   (iim*jjm,   llm)   ) ! potential temperature
  allocate( ORG_rhodz   (iim*jjm,   llm)   )

  allocate( ORG_convm        (iim*jjm,   llm) ) ! mass flux convergence
  allocate( ORG_wflux        (iim*jjm  , llm+1) )
  allocate( ORG_wwuu         (iim*3*jjm, llm+1) )
  allocate( ORG_du           (iim*3*jjm, llm  ) )
  allocate( ORG_dtheta_rhodz (iim*jjm  , llm  ) )
  allocate( ORG_dps          (iim*jjm         ) )

  allocate( ORG_convm_prev        (iim*jjm,   llm) ) ! mass flux convergence
  allocate( ORG_wflux_prev        (iim*jjm  , llm+1) )
  allocate( ORG_wwuu_prev         (iim*3*jjm, llm+1) )
  allocate( ORG_du_prev           (iim*3*jjm, llm  ) )
  allocate( ORG_dtheta_rhodz_prev (iim*jjm  , llm  ) )
  allocate( ORG_dps_prev          (iim*jjm         ) )

  allocate( bp ( llm+1 ) )

  call dumpio_read_data( EX_fid, iim*3*jjm*llm  , ORG_u           (:,:)     )
  call dumpio_read_data( EX_fid, iim*jjm*llm    , ORG_theta       (:,:)     )
  call dumpio_read_data( EX_fid, iim*jjm*llm    , ORG_rhodz       (:,:)     )

  call dumpio_read_data( EX_fid, iim*jjm*llm      , ORG_convm       (:,:)     )
  call dumpio_read_data( EX_fid, iim*jjm*(llm+1)  , ORG_wflux       (:,:)     )
  call dumpio_read_data( EX_fid, iim*3*jjm*(llm+1), ORG_wwuu        (:,:)     )
  call dumpio_read_data( EX_fid, iim*3*jjm*llm    , ORG_du          (:,:)     )
  call dumpio_read_data( EX_fid, iim*jjm*llm      , ORG_dtheta_rhodz(:,:)     )
  call dumpio_read_data( EX_fid, iim*jjm          , ORG_dps         (:)       )

  call dumpio_read_data( EX_fid, iim*jjm*llm      , ORG_convm_prev       (:,:)     )
  call dumpio_read_data( EX_fid, iim*jjm*(llm+1)  , ORG_wflux_prev       (:,:)     )
  call dumpio_read_data( EX_fid, iim*3*jjm*(llm+1), ORG_wwuu_prev        (:,:)     )
  call dumpio_read_data( EX_fid, iim*3*jjm*llm    , ORG_du_prev          (:,:)     )
  call dumpio_read_data( EX_fid, iim*jjm*llm      , ORG_dtheta_rhodz_prev(:,:)     )
  call dumpio_read_data( EX_fid, iim*jjm          , ORG_dps_prev         (:)       )

  call dumpio_read_data( EX_fid, size(bp)          , bp          (:) )
  call dumpio_read_data( EX_fid, 1         , g)

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

  ij_omp_begin = indices(28)
  ij_omp_end   = indices(29)
  ll_beginp1   = indices(30)
  ll_endm1     = indices(31)

  allocate( u       (iim*3*jjm, llm)   ) ! prognostic "velocity"
  allocate( theta   (iim*jjm,   llm)   ) ! potential temperature
  allocate( rhodz   (iim*jjm,   llm)   )

  allocate( convm        (iim*jjm,   llm) ) ! mass flux convergence
  allocate( wflux        (iim*jjm  , llm+1) )
  allocate( wwuu         (iim*3*jjm, llm+1) )
  allocate( du           (iim*3*jjm, llm  ) )
  allocate( dtheta_rhodz (iim*jjm  , llm  ) )
  allocate( dps          (iim*jjm         ) )

  allocate( convm_prev        (iim*jjm,   llm) ) ! mass flux convergence
  allocate( wflux_prev        (iim*jjm  , llm+1) )
  allocate( wwuu_prev         (iim*3*jjm, llm+1) )
  allocate( du_prev           (iim*3*jjm, llm  ) )
  allocate( dtheta_rhodz_prev (iim*jjm  , llm  ) )
  allocate( dps_prev          (iim*jjm         ) )

  u       (:,:) = ORG_u      (:,:) 
  theta   (:,:) = ORG_theta  (:,:) 
  rhodz   (:,:) = ORG_rhodz  (:,:) 

  convm   (:,:) = ORG_convm  (:,:) 
  wflux       (:,:) = ORG_wflux       (:,:) 
  wwuu        (:,:) = ORG_wwuu        (:,:) 
  du          (:,:) = ORG_du          (:,:) 
  dtheta_rhodz(:,:) = ORG_dtheta_rhodz(:,:) 
  dps         (:)   = ORG_dps         (:) 


  convm_prev       (:,:) = ORG_convm_prev       (:,:) 
  wflux_prev       (:,:) = ORG_wflux_prev       (:,:) 
  wwuu_prev        (:,:) = ORG_wwuu_prev        (:,:) 
  du_prev          (:,:) = ORG_du_prev          (:,:) 
  dtheta_rhodz_prev(:,:) = ORG_dtheta_rhodz_prev(:,:) 
  dps_prev         (:)   = ORG_dps_prev         (:) 

  write(*,'(A30,3I6)') 'iim, jjm, llm:',iim,jjm,llm
  write(*,'(A30,2I6)') 'ij_begin, ij_end:',ij_begin,ij_end
  write(*,'(A30,2I6)') 'ij_begin_ext, ij_end_ext:',ij_begin_ext,ij_end_ext
  write(*,'(A30,2I6)') 'll_begin, ll_end:',ll_begin,ll_end
  write(*,'(A30,2I6)') 'll_beginp1, ll_endm1:',ll_beginp1,ll_endm1

  write(*,'(A30,3I6)') 't_right, t_rup, t_lup:',t_right,t_rup,t_lup
  write(*,'(A30,3I6)') 't_left, t_ldown, t_rdown:',t_left,t_ldown,t_rdown

  write(*,'(A30,3I6)') 'u_right, u_rup, u_lup:',u_right,u_rup,u_lup
  write(*,'(A30,3I6)') 'u_left, u_ldown, u_rdown:',u_left,u_ldown,u_rdown

  write(*,'(A30,3I6)') 'z_rup, z_up, z_lup:',z_rup,z_up,z_lup
  write(*,'(A30,3I6)') 'z_ldown, z_down, z_rdown:',z_ldown,z_down,z_rdown

  write(*,'(A30,F15.8)') 'dbg: g:',g

!!$  call DEBUG_valuecheck( 'convm_prev        :', convm_prev        (:,:) )
!!$  call DEBUG_valuecheck( 'wflux_prev        :', wflux_prev        (:,:) )
!!$  call DEBUG_valuecheck( 'wwuu_prev         :', wwuu_prev         (:,:) )
!!$  call DEBUG_valuecheck( 'du_prev           :', du_prev           (:,:) )
!!$  call DEBUG_valuecheck( 'dtheta_rhodz_prev :', dtheta_rhodz_prev (:,:) )
!!$  call DEBUG_valuecheck( 'dps_prev          :', dps_prev          (:) )

  call DEBUG_valuecheck( 'bp   ',bp(:) )

  write(*,*) "*** Finish initialize"

  !###############################################################################

  write(*,*) "*** Start kernel"

  do iteration = 1, SET_iteration

    convm       (:,:) = convm_prev       (:,:) 
    wflux       (:,:) = wflux_prev       (:,:) 
    wwuu        (:,:) = wwuu_prev        (:,:) 
    du          (:,:) = du_prev          (:,:) 
    dtheta_rhodz(:,:) = dtheta_rhodz_prev(:,:) 
    dps         (:)   = dps_prev         (:) 

    call DEBUG_rapstart('MAIN_comp_caldyn_vert')

    call compute_caldyn_vert( & 
         & u,           &    ! [IN]
         & theta,       &    ! [IN]
         & rhodz,       &    ! [IN]
         & convm,       &    ! [INOUT]
         & wflux,       &    ! [INOUT]
         & wwuu,        &    ! [INOUT]
         & dps,         &    ! [INOUT]
         & dtheta_rhodz,&    ! [INOUT]
         & du           &    ! [INOUT]
         &)


    call DEBUG_rapend('MAIN_comp_caldyn_vert')

    if ( SET_iteration<5 .or. mod(iteration, SET_iteration/1)==0 ) then
      write(ADM_LOG_FID,*) '### check point iteration:',iteration
      write(ADM_LOG_FID,*) '### Input ###'

      call DEBUG_valuecheck(  'u       ', u       (:,:) )
      call DEBUG_valuecheck(  'theta   ', theta   (:,:) )
      call DEBUG_valuecheck(  'rhodz   ', rhodz   (:,:) )

      call DEBUG_valuecheck(  'convm_prev        ', convm_prev       (:,:) )
      call DEBUG_valuecheck(  'wflux_prev        ', wflux_prev       (:,:) )
      call DEBUG_valuecheck(  'wwuu_prev         ', wwuu_prev        (:,:) )
      call DEBUG_valuecheck(  'du_prev           ', du_prev          (:,:) )
      call DEBUG_valuecheck(  'dtheta_rhodz_prev ', dtheta_rhodz_prev(:,:) )
      call DEBUG_valuecheck(  'dps_prev          ', dps_prev         (:) )

      write(ADM_LOG_FID,*) '### Output ###'
      call DEBUG_valuecheck(  'convm        ', convm       (:,:) )
      call DEBUG_valuecheck(  'wflux        ', wflux       (:,:) )
      call DEBUG_valuecheck(  'wwuu         ', wwuu        (:,:) )
      call DEBUG_valuecheck(  'du           ', du          (:,:) )
      call DEBUG_valuecheck(  'dtheta_rhodz ', dtheta_rhodz(:,:) )
      call DEBUG_valuecheck(  'dps          ', dps         (:) )
    end if
  end do

  write(ADM_LOG_FID,*) '### final iteration:',iteration-1
  write(ADM_LOG_FID,*) '### Validation : point-by-point diff ###'

  convm       (:,:) = ORG_convm       (:,:) - convm       (:,:)
  wflux       (:,:) = ORG_wflux       (:,:) - wflux       (:,:)
  wwuu        (:,:) = ORG_wwuu        (:,:) - wwuu        (:,:)
  du          (:,:) = ORG_du          (:,:) - du          (:,:)
  dtheta_rhodz(:,:) = ORG_dtheta_rhodz(:,:) - dtheta_rhodz(:,:)
  dps         (:)   = ORG_dps         (:)   - dps         (:) 

  call DEBUG_valuecheck(  'convm        ', convm       (:,:) )
  call DEBUG_valuecheck(  'wflux        ', wflux       (:,:) )
  call DEBUG_valuecheck(  'wwuu         ', wwuu        (:,:) )
  call DEBUG_valuecheck(  'du           ', du          (:,:) )
  call DEBUG_valuecheck(  'dtheta_rhodz ', dtheta_rhodz(:,:) )
  call DEBUG_valuecheck(  'dps          ', dps         (:) )


  write(*,*) "*** Finish kernel"

  !###############################################################################

  call DEBUG_rapreport

end program comp_caldyn_vert




!!$ Local Variables:
!!$ mode: f90
!!$ f90-do-indent: 2
!!$ f90-if-indent: 2
!!$ f90-program-indent: 2
!!$ f90-type-indent: 2
!!$ f90-directive-comment-re: "!\$omp "
!!$ end:

