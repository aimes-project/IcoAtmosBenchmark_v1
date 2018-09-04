MODULE caldyn_gcm_mod
   use mod_misc
   PRIVATE
   TYPE(t_message) :: req_ps, req_mass, req_theta_rhodz, req_u, req_qu

   public compute_caldyn_vert

CONTAINS

   SUBROUTINE compute_caldyn_vert(u,theta,rhodz,convm, wflux,wwuu, dps,dtheta_rhodz,du)
      use prec
      use mod_misc
      IMPLICIT NONE
      REAL(rstd),INTENT(IN)  :: u(iim*3*jjm,llm)
      REAL(rstd),INTENT(IN)  :: theta(iim*jjm,llm)
      REAL(rstd),INTENT(IN)  :: rhodz(iim*jjm,llm)
      REAL(rstd),INTENT(INOUT)  :: convm(iim*jjm,llm)  ! mass flux convergence

      REAL(rstd),INTENT(INOUT) :: wflux(iim*jjm,llm+1) ! vertical mass flux (kg/m2/s)
      REAL(rstd),INTENT(INOUT) :: wwuu(iim*3*jjm,llm+1)
      REAL(rstd),INTENT(INOUT) :: du(iim*3*jjm,llm)
      REAL(rstd),INTENT(INOUT) :: dtheta_rhodz(iim*jjm,llm)
      REAL(rstd),INTENT(INOUT) :: dps(iim*jjm)

      ! temporary variable    
      INTEGER :: i,j,ij,l
      REAL(rstd) :: p_ik, exner_ik
!!$!!$KERN    INTEGER,SAVE ::ij_omp_begin, ij_omp_end
!!$!!$KERN!$OMP THREADPRIVATE(ij_omp_begin, ij_omp_end)
!!$!!$KERN    LOGICAL,SAVE :: first=.TRUE.
      LOGICAL,SAVE :: first=.false.
      !$OMP THREADPRIVATE(first)

      CALL trace_start("compute_geopot")

!!$!!##KERNEL: this section is never called when kernelize.
!!$!!$KERN    IF (first) THEN
!!$!!$KERN      first=.FALSE.
!!$!!$KERN      CALL distrib_level(ij_end-ij_begin+1,ij_omp_begin,ij_omp_end)
!!$!!$KERN      ij_omp_begin=ij_omp_begin+ij_begin-1
!!$!!$KERN      ij_omp_end=ij_omp_end+ij_begin-1
!!$!!$KERN    ENDIF

!    REAL(rstd) :: wwuu(iim*3*jjm,llm+1) ! tmp var, don't know why but gain 30% on the whole code in opemp
      ! need to be understood

      !    wwuu=wwuu_out
      CALL trace_start("compute_caldyn_vert")

      !$OMP BARRIER   
!!! cumulate mass flux convergence from top to bottom
      !  IF (is_omp_level_master) THEN
      DO  l = llm-1, 1, -1
         !    IF (caldyn_conserv==energy) CALL test_message(req_qu) 

!!$OMP DO SCHEDULE(STATIC) 
         !DIR$ SIMD
         DO ij=ij_omp_begin,ij_omp_end         
            convm(ij,l) = convm(ij,l) + convm(ij,l+1)
         ENDDO
      ENDDO
      !  ENDIF

      !$OMP BARRIER
      ! FLUSH on convm
!!!!!!!!!!!!!!!!!!!!!!!!!

      ! compute dps
      IF (is_omp_first_level) THEN
         !DIR$ SIMD
         DO ij=ij_begin,ij_end         
            ! dps/dt = -int(div flux)dz
            dps(ij) = convm(ij,1) * g 
         ENDDO
      ENDIF

!!! Compute vertical mass flux (l=1,llm+1 done by caldyn_BC)
      DO l=ll_beginp1,ll_end
         !    IF (caldyn_conserv==energy) CALL test_message(req_qu) 
         !DIR$ SIMD
         DO ij=ij_begin,ij_end         
            ! w = int(z,ztop,div(flux)dz) + B(eta)dps/dt
            ! => w>0 for upward transport
            wflux( ij, l ) = bp(l) * convm( ij, 1 ) - convm( ij, l ) 
         ENDDO
      ENDDO

      !--> flush wflux  
      !$OMP BARRIER

      DO l=ll_begin,ll_endm1
         !DIR$ SIMD
         DO ij=ij_begin,ij_end         
            dtheta_rhodz(ij, l ) = dtheta_rhodz(ij, l )  -   0.5 * (  wflux(ij,l+1) * (theta(ij,l) + theta(ij,l+1)))
         ENDDO
      ENDDO

      DO l=ll_beginp1,ll_end
         !DIR$ SIMD
         DO ij=ij_begin,ij_end         
            dtheta_rhodz(ij, l ) = dtheta_rhodz(ij, l )  +   0.5 * ( wflux(ij,l  ) * (theta(ij,l-1) + theta(ij,l) ) )
         ENDDO
      ENDDO


      ! Compute vertical transport
      DO l=ll_beginp1,ll_end
         !DIR$ SIMD
         DO ij=ij_begin,ij_end         
            wwuu(ij+u_right,l) = 0.5*( wflux(ij,l) + wflux(ij+t_right,l)) * (u(ij+u_right,l) - u(ij+u_right,l-1))
            wwuu(ij+u_lup,l) = 0.5* ( wflux(ij,l) + wflux(ij+t_lup,l)) * (u(ij+u_lup,l) - u(ij+u_lup,l-1))
            wwuu(ij+u_ldown,l) = 0.5*( wflux(ij,l) + wflux(ij+t_ldown,l)) * (u(ij+u_ldown,l) - u(ij+u_ldown,l-1))
         ENDDO
      ENDDO

      !--> flush wwuu  
      !$OMP BARRIER

      ! Add vertical transport to du
      DO l=ll_begin,ll_end
         !DIR$ SIMD
         DO ij=ij_begin,ij_end         
            du(ij+u_right, l )   = du(ij+u_right,l)  - (wwuu(ij+u_right,l+1)+ wwuu(ij+u_right,l)) / (rhodz(ij,l)+rhodz(ij+t_right,l))
            du(ij+u_lup, l )     = du(ij+u_lup,l)    - (wwuu(ij+u_lup,l+1)  + wwuu(ij+u_lup,l))   / (rhodz(ij,l)+rhodz(ij+t_lup,l))
            du(ij+u_ldown, l )   = du(ij+u_ldown,l)  - (wwuu(ij+u_ldown,l+1)+ wwuu(ij+u_ldown,l)) / (rhodz(ij,l)+rhodz(ij+t_ldown,l))
         ENDDO
      ENDDO

      !  DO l=ll_beginp1,ll_end
      !!DIR$ SIMD
      !      DO ij=ij_begin,ij_end         
      !        wwuu_out(ij+u_right,l) = wwuu(ij+u_right,l)
      !        wwuu_out(ij+u_lup,l)   = wwuu(ij+u_lup,l)
      !        wwuu_out(ij+u_ldown,l) = wwuu(ij+u_ldown,l)
      !     ENDDO
      !   ENDDO

      CALL trace_end("compute_caldyn_vert")

   END SUBROUTINE compute_caldyn_vert

END MODULE caldyn_gcm_mod

