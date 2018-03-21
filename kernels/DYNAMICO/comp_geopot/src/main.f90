!-------------------------------------------------------------------------------
!
!+  
!
!-------------------------------------------------------------------------------
program comp_geopot
  use prec
  use mod_misc
  use caldyn_gcm_mod, only : compute_geopot
  implicit none


  real(rstd), allocatable :: ORG_ps          (:)  
  real(rstd), allocatable :: ORG_rhodz       (:,:)
  real(rstd), allocatable :: ORG_theta       (:,:)
  real(rstd), allocatable :: ORG_pk          (:,:)
  real(rstd), allocatable :: ORG_geopot      (:,:)

  real(rstd), allocatable :: ORG_ps_prev     (:)
  real(rstd), allocatable :: ORG_pk_prev     (:,:)
  real(rstd), allocatable :: ORG_geopot_prev (:,:)

  real(rstd), allocatable :: ps          (:)  
  real(rstd), allocatable :: rhodz       (:,:)
  real(rstd), allocatable :: theta       (:,:)
  real(rstd), allocatable :: pk          (:,:)
  real(rstd), allocatable :: geopot      (:,:)

  real(rstd), allocatable :: ps_prev     (:)
  real(rstd), allocatable :: pk_prev     (:,:)
  real(rstd), allocatable :: geopot_prev (:,:)


  integer :: iteration
!!$       indices(:) = (/&
!!$            & iim,jjm,llm,&
!!$            & ij_begin,ij_end,&
!!$            & ij_begin_ext,ij_end_ext,&
!!$            & ll_begin,ll_end,&
!!$            & t_right, t_rup, t_lup, t_left, t_ldown, t_rdown,&
!!$            & u_right, u_rup, u_lup, u_left, u_ldown, u_rdown &
!!$            & caldyn_eta, >0 if(boussinesq)
!!$            & /)
  integer,parameter :: indices_size = 29
  integer(i8) :: indices(indices_size)


  !=============================================================================

  write(*,*) "[KERNEL] comp_geopot"

  !###############################################################################

  write(*,*) "*** Start  initialize"



  !###############################################################################

  !---< read input data >---

  call dumpio_syscheck
  call dumpio_mk_fname(EX_fname,'snapshot.comp_geopot','pe',SET_prc_me-1,6)
  call dumpio_fopen(EX_fid,EX_fname,IO_FREAD)

  call dumpio_read_data( EX_fid, size(indices), indices )

  iim = indices(1)
  jjm = indices(2)
  llm = indices(3)
!!$   write(*,'(A30,3I6)') 'dbg: iim,jjm,llm:',iim,jjm,llm

  allocate( ORG_ps          (iim*jjm        )     )
  allocate( ORG_rhodz       (iim*jjm,llm    )     )
  allocate( ORG_theta       (iim*jjm,llm    )     )
  allocate( ORG_pk          (iim*jjm,llm    )     )
  allocate( ORG_geopot      (iim*jjm,(llm+1))     )

  allocate( ORG_ps_prev          (iim*jjm        )     )
  allocate( ORG_pk_prev          (iim*jjm,llm    )     )
  allocate( ORG_geopot_prev      (iim*jjm,(llm+1))     )

  allocate( mass_ak (llm) )
  allocate( mass_bk (llm) )

  call dumpio_read_data( EX_fid, iim*jjm        , ORG_ps          (:)       )
  call dumpio_read_data( EX_fid, iim*jjm*llm    , ORG_rhodz       (:,:)     )
  call dumpio_read_data( EX_fid, iim*jjm*llm    , ORG_theta       (:,:)     )
  call dumpio_read_data( EX_fid, iim*jjm*llm    , ORG_pk          (:,:)     )
  call dumpio_read_data( EX_fid, iim*jjm*(llm+1), ORG_geopot      (:,:)     )

  call dumpio_read_data( EX_fid, iim*jjm        , ORG_ps_prev          (:)       )
  call dumpio_read_data( EX_fid, iim*jjm*llm    , ORG_pk_prev          (:,:)     )
  call dumpio_read_data( EX_fid, iim*jjm*(llm+1), ORG_geopot_prev      (:,:)     )


  call dumpio_read_data( EX_fid, 1  , ptop )
  call dumpio_read_data( EX_fid, llm, mass_ak(:) )
  call dumpio_read_data( EX_fid, llm, mass_bk(:) )
  call dumpio_read_data( EX_fid, 1  , g     )
  call dumpio_read_data( EX_fid, 1  , kappa )
  call dumpio_read_data( EX_fid, 1  , cpp   )
  call dumpio_read_data( EX_fid, 1  , preff )

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
  boussinesq = ( indices(29) > 0 )

  allocate( ps          (iim*jjm        )     )
  allocate( rhodz       (iim*jjm,llm    )     )
  allocate( theta       (iim*jjm,llm    )     )
  allocate( pk          (iim*jjm,llm    )     )
  allocate( geopot      (iim*jjm,(llm+1))     )

  allocate( ps_prev          (iim*jjm        )     )
  allocate( pk_prev          (iim*jjm,llm    )     )
  allocate( geopot_prev      (iim*jjm,(llm+1))     )

  ps     (:)   =   ORG_ps     (:)  
  rhodz  (:,:) =   ORG_rhodz  (:,:)
  theta  (:,:) =   ORG_theta  (:,:)
  pk     (:,:) =   ORG_pk     (:,:)
  geopot (:,:) =   ORG_geopot (:,:)

  ps_prev     (:)   =   ORG_ps_prev     (:)  
  pk_prev     (:,:) =   ORG_pk_prev     (:,:)
  geopot_prev (:,:) =   ORG_geopot_prev (:,:)

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
  write(*,'(A30, L6)') 'boussinesq:',boussinesq 
  write(*,'(A30,F15.8)') 'g:',g

  call DEBUG_valuecheck(  'mass_ak   ', mass_ak   (:) )
  call DEBUG_valuecheck(  'mass_bk   ', mass_bk   (:) )

  write(*,*) "*** Finish initialize"

  !###############################################################################

  write(*,*) "*** Start kernel"

  do iteration = 1, SET_iteration

    ps     (:)   = ps_prev     (:)
    pk     (:,:) = pk_prev     (:,:)
    geopot (:,:) = geopot_prev (:,:)

    call DEBUG_rapstart('MAIN_comp_geopot')

    call compute_geopot( & 
         & ps,          & ! [INOUT]
         & rhodz,       & ! [IN]
         & theta,       & ! [IN]
         & pk,          & ! [INOUT]
         & geopot       & ! [INOUT]
         )


    call DEBUG_rapend('MAIN_comp_geopot')

    if ( SET_iteration<5 .or. mod(iteration, SET_iteration/1)==0 ) then
      write(ADM_LOG_FID,*) '### check point iteration:',iteration
      write(ADM_LOG_FID,*) '### Input ###'

      call DEBUG_valuecheck(  'ps_prev    ', ps_prev     (:)   )
      call DEBUG_valuecheck(  'rhodz      ', rhodz       (:,:) )
      call DEBUG_valuecheck(  'theta      ', theta       (:,:) )
      call DEBUG_valuecheck(  'pk_prev    ', pk_prev     (:,:) )
      call DEBUG_valuecheck(  'geopot_prev', geopot_prev (:,:) )


      write(ADM_LOG_FID,*) '### Output ###'
      call DEBUG_valuecheck(  'ps         ', ps          (:)   )
      call DEBUG_valuecheck(  'pk         ', pk          (:,:) )
      call DEBUG_valuecheck(  'geopot     ', geopot      (:,:) )
    end if
  end do

  write(ADM_LOG_FID,*) '### final iteration:',iteration-1
  write(ADM_LOG_FID,*) '### Validation : grid-by-grid diff ###'

  ps     (:)   = ORG_ps     (:)   - ps     (:)
  pk     (:,:) = ORG_pk     (:,:) - pk     (:,:) 
  geopot (:,:) = ORG_geopot (:,:) - geopot (:,:) 

  call DEBUG_valuecheck( 'ps    ', ps     (:) )
  call DEBUG_valuecheck( 'pk    ', pk     (:,:) )
  call DEBUG_valuecheck( 'geopot', geopot (:,:) )

  write(*,*) "*** Finish kernel"

  !###############################################################################

  call DEBUG_rapreport

end program comp_geopot




!!$ Local Variables:
!!$ mode: f90
!!$ f90-do-indent: 2
!!$ f90-if-indent: 2
!!$ f90-program-indent: 2
!!$ f90-type-indent: 2
!!$ f90-directive-comment-re: "!\$omp "
!!$ end:

