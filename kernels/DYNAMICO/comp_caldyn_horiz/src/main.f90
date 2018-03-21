!-------------------------------------------------------------------------------
!
!+  
!
!-------------------------------------------------------------------------------
program comp_caldyn_horiz
  use prec
  use mod_misc
  use caldyn_gcm_mod, only : compute_caldyn_horiz
  implicit none

  real(rstd),allocatable :: ORG_u      (:,:) ! prognostic "velocity"
  real(rstd),allocatable :: ORG_rhodz  (:,:)
  real(rstd),allocatable :: ORG_qu     (:,:)
  real(rstd),allocatable :: ORG_theta  (:,:) ! potential temperature
  real(rstd),allocatable :: ORG_pk     (:,:) ! exner function
  real(rstd),allocatable :: ORG_geopot (:,:) ! geopotential

  real(rstd),allocatable :: ORG_hflux        (:,:) ! hflux in kg/s
  real(rstd),allocatable :: ORG_convm        (:,:) ! mass flux convergence
  real(rstd),allocatable :: ORG_dtheta_rhodz (:,:)
  real(rstd),allocatable :: ORG_du           (:,:)

  real(rstd),allocatable :: ORG_pk_prev           (:,:) ! Exner function
  real(rstd),allocatable :: ORG_hflux_prev        (:,:) ! hflux in kg/s
  real(rstd),allocatable :: ORG_convm_prev        (:,:) ! mass flux convergence
  real(rstd),allocatable :: ORG_dtheta_rhodz_prev (:,:)
  real(rstd),allocatable :: ORG_du_prev           (:,:)

  real(rstd),allocatable :: u       (:,:) ! prognostic "velocity"
  real(rstd),allocatable :: rhodz   (:,:)
  real(rstd),allocatable :: qu      (:,:)
  real(rstd),allocatable :: theta   (:,:) ! potential temperature
  real(rstd),allocatable :: pk      (:,:) ! Exner function
  real(rstd),allocatable :: geopot  (:,:) ! geopotential

  real(rstd),allocatable :: hflux        (:,:) ! hflux in kg/s
  real(rstd),allocatable :: convm        (:,:) ! mass flux convergence
  real(rstd),allocatable :: dtheta_rhodz (:,:)
  real(rstd),allocatable :: du           (:,:)

  real(rstd),allocatable :: pk_prev           (:,:) ! Exner function
  real(rstd),allocatable :: hflux_prev        (:,:) ! hflux in kg/s
  real(rstd),allocatable :: convm_prev        (:,:) ! mass flux convergence
  real(rstd),allocatable :: dtheta_rhodz_prev (:,:)
  real(rstd),allocatable :: du_prev           (:,:)


  integer :: iteration
!!$       indices(:) = (/&
!!$            & iim,jjm,llm,&
!!$            & ij_begin,ij_end,&
!!$            & ij_begin_ext,ij_end_ext,&
!!$            & ll_begin,ll_end,&
!!$            & t_right, t_rup, t_lup, t_left, t_ldown, t_rdown,&
!!$            & u_right, u_rup, u_lup, u_left, u_ldown, u_rdown &
!!$            & caldyn_conserv, >0 if(boussinesq)
!!$            & /)
  integer,parameter :: indices_size = 29
  integer(i8) :: indices(indices_size)


  !=============================================================================

  write(*,*) "[KERNEL] comp_caldyn_horiz"

  !###############################################################################

  write(*,*) "*** Start  initialize"



  !###############################################################################

  !---< read input data >---

  call dumpio_syscheck
  call dumpio_mk_fname(EX_fname,'snapshot.comp_caldyn_horiz','pe',SET_prc_me-1,6)
  call dumpio_fopen(EX_fid,EX_fname,IO_FREAD)

  call dumpio_read_data( EX_fid, size(indices), indices )

  iim = indices(1)
  jjm = indices(2)
  llm = indices(3)
!!$   write(*,'(A30,3I6)') 'dbg: iim,jjm,llm:',iim,jjm,llm

  allocate( ORG_u       (iim*3*jjm, llm)   ) ! prognostic "velocity"
  allocate( ORG_rhodz   (iim*jjm,   llm)   )
  allocate( ORG_qu      (iim*3*jjm, llm)   )
  allocate( ORG_theta   (iim*jjm,   llm)   ) ! potential temperature
  allocate( ORG_pk      (iim*jjm,   llm)   ) ! Exner function
  allocate( ORG_geopot  (iim*jjm,   llm+1) ) ! geopotential

  allocate( ORG_hflux        (iim*3*jjm, llm) ) ! hflux in kg/s
  allocate( ORG_convm        (iim*jjm,   llm) ) ! mass flux convergence
  allocate( ORG_dtheta_rhodz (iim*jjm,   llm) )
  allocate( ORG_du           (iim*3*jjm, llm) )

  allocate( ORG_pk_prev           (iim*jjm,   llm) ) ! Exner function
  allocate( ORG_hflux_prev        (iim*3*jjm, llm) ) ! hflux in kg/s
  allocate( ORG_convm_prev        (iim*jjm,   llm) ) ! mass flux convergence
  allocate( ORG_dtheta_rhodz_prev (iim*jjm,   llm) )
  allocate( ORG_du_prev           (iim*3*jjm, llm) )

  allocate( le          (iim*3*jjm) ) !! field_type is field_U
  allocate( Ai          (iim*jjm)   ) !! field_type is field_T

  allocate( de          (iim*3*jjm) ) !! field_type is field_U
  allocate( Av          (iim*2*jjm) ) !! field_type is field_Z
  allocate( Wee         (iim*3*jjm,5,2) ) !! field_type is field_U

!!$   !!KERNEL: Are these initialization necessary ?
!!$   ORG_u      (:,:) = 0.0_rstd
!!$   ORG_rhodz  (:,:) = 0.0_rstd
!!$   ORG_qu     (:,:) = 0.0_rstd
!!$   ORG_theta  (:,:) = 0.0_rstd
!!$   ORG_pk     (:,:) = 0.0_rstd
!!$   ORG_geopot (:,:) = 0.0_rstd
!!$
!!$   ORG_hflux        (:,:) = 0.0_rstd
!!$   ORG_convm        (:,:) = 0.0_rstd
!!$   ORG_dtheta_rhodz (:,:) = 0.0_rstd
!!$   ORG_du           (:,:) = 0.0_rstd
!!$
!!$   ORG_pk_prev           (:,:) = 0.0_rstd
!!$   ORG_hflux_prev        (:,:) = 0.0_rstd
!!$   ORG_convm_prev        (:,:) = 0.0_rstd
!!$   ORG_dtheta_rhodz_prev (:,:) = 0.0_rstd
!!$   ORG_du_prev           (:,:) = 0.0_rstd
!!$
!!$
!!$   le  (:) = 0.0_rstd
!!$   Ai  (:) = 0.0_rstd
!!$   de  (:) = 0.0_rstd
!!$   Av  (:) = 0.0_rstd
!!$   Wee (:,:,:) = 0.0_rstd

  call dumpio_read_data( EX_fid, iim*3*jjm*llm  , ORG_u           (:,:)     )
  call dumpio_read_data( EX_fid, iim*jjm*llm    , ORG_rhodz       (:,:)     )
  call dumpio_read_data( EX_fid, iim*3*jjm*llm  , ORG_qu          (:,:)     )
  call dumpio_read_data( EX_fid, iim*jjm*llm    , ORG_theta       (:,:)     )
  call dumpio_read_data( EX_fid, iim*jjm*llm    , ORG_pk          (:,:)     )
  call dumpio_read_data( EX_fid, iim*jjm*(llm+1), ORG_geopot      (:,:)     )

  call dumpio_read_data( EX_fid, iim*3*jjm*llm  , ORG_hflux       (:,:)     )
  call dumpio_read_data( EX_fid, iim*jjm*llm    , ORG_convm       (:,:)     )
  call dumpio_read_data( EX_fid, iim*jjm*llm    , ORG_dtheta_rhodz(:,:)     )
  call dumpio_read_data( EX_fid, iim*3*jjm*llm  , ORG_du          (:,:)     )

  call dumpio_read_data( EX_fid, iim*jjm*llm    , ORG_pk_prev          (:,:)     )
  call dumpio_read_data( EX_fid, iim*3*jjm*llm  , ORG_hflux_prev       (:,:)     )
  call dumpio_read_data( EX_fid, iim*jjm*llm    , ORG_convm_prev       (:,:)     )
  call dumpio_read_data( EX_fid, iim*jjm*llm    , ORG_dtheta_rhodz_prev(:,:)     )
  call dumpio_read_data( EX_fid, iim*3*jjm*llm  , ORG_du_prev          (:,:)     )

  call dumpio_read_data( EX_fid, size(le)  , le (:)    )
  call dumpio_read_data( EX_fid, size(Ai)  , Ai (:)    )
  call dumpio_read_data( EX_fid, size(de)  , de (:)    )
  call dumpio_read_data( EX_fid, size(Av)  , Av (:)    )
  call dumpio_read_data( EX_fid, size(Wee) , Wee(:,:,:))
  call dumpio_read_data( EX_fid, 1         , g )

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

  caldyn_conserv = indices(28)
  boussinesq = ( indices(29) > 0 )

  allocate( u       (iim*3*jjm, llm)   ) ! prognostic "velocity"
  allocate( rhodz   (iim*jjm,   llm)   )
  allocate( qu      (iim*3*jjm, llm)   )
  allocate( theta   (iim*jjm,   llm)   ) ! potential temperature
  allocate( pk      (iim*jjm,   llm)   ) ! Exner function
  allocate( geopot  (iim*jjm,   llm+1) ) ! geopotential

  allocate( hflux        (iim*3*jjm, llm) ) ! hflux in kg/s
  allocate( convm        (iim*jjm,   llm) ) ! mass flux convergence
  allocate( dtheta_rhodz (iim*jjm,   llm) )
  allocate( du           (iim*3*jjm, llm) )

  allocate( pk_prev           (iim*jjm,   llm) ) ! Exner function
  allocate( hflux_prev        (iim*3*jjm, llm) ) ! hflux in kg/s
  allocate( convm_prev        (iim*jjm,   llm) ) ! mass flux convergence
  allocate( dtheta_rhodz_prev (iim*jjm,   llm) )
  allocate( du_prev           (iim*3*jjm, llm) )

  u       (:,:) = ORG_u      (:,:) 
  rhodz   (:,:) = ORG_rhodz  (:,:) 
  qu      (:,:) = ORG_qu     (:,:) 
  theta   (:,:) = ORG_theta  (:,:) 
  pk      (:,:) = ORG_pk     (:,:) 
  geopot  (:,:) = ORG_geopot (:,:) 

  hflux        (:,:) = ORG_hflux        (:,:) 
  convm        (:,:) = ORG_convm        (:,:) 
  dtheta_rhodz (:,:) = ORG_dtheta_rhodz (:,:) 
  du           (:,:) = ORG_du           (:,:) 

  pk_prev           (:,:) = ORG_pk_prev           (:,:) 
  hflux_prev        (:,:) = ORG_hflux_prev        (:,:) 
  convm_prev        (:,:) = ORG_convm_prev        (:,:) 
  dtheta_rhodz_prev (:,:) = ORG_dtheta_rhodz_prev (:,:) 
  du_prev           (:,:) = ORG_du_prev           (:,:) 

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

  write(*,'(A30, I6)') 'caldyn_conserv:',caldyn_conserv
  write(*,'(A30, L6)') 'boussinesq:',boussinesq
  write(*,'(A30,F15.8)') 'g:',g


  call DEBUG_valuecheck( 'pk_prev          ', pk_prev           (:,:) )
  call DEBUG_valuecheck( 'hflux_prev       ', hflux_prev        (:,:) )
  call DEBUG_valuecheck( 'convm_prev       ', convm_prev        (:,:) )
  call DEBUG_valuecheck( 'dtheta_rhodz_prev', dtheta_rhodz_prev (:,:) )
  call DEBUG_valuecheck( 'du_prev          ', du_prev           (:,:) )

  call DEBUG_valuecheck( 'le   ',le(:) )
  call DEBUG_valuecheck( 'Ai   ',Ai(:) )
  call DEBUG_valuecheck( 'de   ',de(:) )
  call DEBUG_valuecheck( 'Av   ',Av(:) )
  call DEBUG_valuecheck( 'Wee  ',Wee(:,:,:) )

  write(*,*) "*** Finish initialize"

  !###############################################################################

  write(*,*) "*** Start kernel"

  do iteration = 1, SET_iteration

    pk           (:,:) = pk_prev           (:,:) 
    hflux        (:,:) = hflux_prev        (:,:) 
    convm        (:,:) = convm_prev        (:,:) 
    dtheta_rhodz (:,:) = dtheta_rhodz_prev (:,:) 
    du           (:,:) = du_prev           (:,:) 

    call DEBUG_rapstart('MAIN_comp_caldyn_horiz')

    call compute_caldyn_horiz( & 
         & u,            &    ! [IN]
         & rhodz,        &    ! [IN]
         & qu,           &    ! [IN]
         & theta,        &    ! [IN]
         & pk,           &    ! [INOUT]
         & geopot,       &    ! [IN]
         & hflux,        &    ! [INOUT]
         & convm,        &    ! [INOUT]
         & dtheta_rhodz, &    ! [INOUT]
         & du            )    ! [INOUT]


    call DEBUG_rapend('MAIN_comp_caldyn_horiz')

    if ( SET_iteration<5 .or. mod(iteration, SET_iteration/1)==1 ) then
      write(ADM_LOG_FID,*) '### check point iteration:',iteration
      write(ADM_LOG_FID,*) '### Input ###'

      call DEBUG_valuecheck(  'u       ', u       (:,:) )
      call DEBUG_valuecheck(  'rhodz   ', rhodz   (:,:) )
      call DEBUG_valuecheck(  'qu      ', qu      (:,:) )
      call DEBUG_valuecheck(  'theta   ', theta   (:,:) )
      call DEBUG_valuecheck(  'geopot  ', geopot  (:,:) )

      call DEBUG_valuecheck(  'pk_prev          ', pk_prev           (:,:) )
      call DEBUG_valuecheck(  'hflux_prev       ', hflux_prev        (:,:) )
      call DEBUG_valuecheck(  'convm_prev       ', convm_prev        (:,:) )
      call DEBUG_valuecheck(  'dtheta_rhodz_prev', dtheta_rhodz_prev (:,:) )
      call DEBUG_valuecheck(  'du_prev          ', du_prev           (:,:) )

      write(ADM_LOG_FID,*) '### Output ###'
      call DEBUG_valuecheck(  'pk          ', pk           (:,:) )
      call DEBUG_valuecheck(  'hflux       ', hflux        (:,:) )
      call DEBUG_valuecheck(  'convm       ', convm        (:,:) )
      call DEBUG_valuecheck(  'dtheta_rhodz', dtheta_rhodz (:,:) )
      call DEBUG_valuecheck(  'du          ', du           (:,:) )
    end if
  end do

  write(ADM_LOG_FID,*) '### final iteration:',iteration-1
  write(ADM_LOG_FID,*) '### Validation : grid-by-grid diff ###'

  pk           (:,:) = ORG_pk           (:,:) - pk           (:,:) 
  hflux        (:,:) = ORG_hflux        (:,:) - hflux        (:,:) 
  convm        (:,:) = ORG_convm        (:,:) - convm        (:,:) 
  dtheta_rhodz (:,:) = ORG_dtheta_rhodz (:,:) - dtheta_rhodz (:,:) 
  du           (:,:) = ORG_du           (:,:) - du           (:,:) 

  call DEBUG_valuecheck(  'pk          ', pk           (:,:) )
  call DEBUG_valuecheck(  'hflux       ', hflux        (:,:) )
  call DEBUG_valuecheck(  'convm       ', convm        (:,:) )
  call DEBUG_valuecheck(  'dtheta_rhodz', dtheta_rhodz (:,:) )
  call DEBUG_valuecheck(  'du          ', du           (:,:) )

  write(*,*) "*** Finish kernel"

  !###############################################################################

  call DEBUG_rapreport

end program comp_caldyn_horiz




!!$ Local Variables:
!!$ mode: f90
!!$ f90-do-indent: 2
!!$ f90-if-indent: 2
!!$ f90-program-indent: 2
!!$ f90-type-indent: 2
!!$ f90-directive-comment-re: "!\$omp "
!!$ end:

