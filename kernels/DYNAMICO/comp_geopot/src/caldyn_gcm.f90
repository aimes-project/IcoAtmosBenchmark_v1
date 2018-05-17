MODULE caldyn_gcm_mod
   use mod_misc
   PRIVATE
   TYPE(t_message) :: req_ps, req_mass, req_theta_rhodz, req_u, req_qu

   public compute_geopot

CONTAINS

   SUBROUTINE compute_geopot(ps,rhodz,theta, pk,geopot)
      use prec
      use mod_misc
      IMPLICIT NONE
      REAL(rstd),INTENT(INOUT) :: ps(iim*jjm)
      REAL(rstd),INTENT(IN)    :: rhodz(iim*jjm,llm)
      REAL(rstd),INTENT(IN)    :: theta(iim*jjm,llm)    ! potential temperature
      REAL(rstd),INTENT(INOUT) :: pk(iim*jjm,llm)       ! Exner function
      REAL(rstd),INTENT(INOUT) :: geopot(iim*jjm,llm+1) ! geopotential

      INTEGER :: i,j,ij,l
      REAL(rstd) :: p_ik, exner_ik

      CALL trace_start("compute_geopot")

      IF(caldyn_eta==eta_mass) THEN

!!! Compute exner function and geopotential
         DO l = 1,llm
            !DIR$ SIMD
            DO ij=ij_begin_ext,ij_end_ext
               p_ik = ptop + mass_ak(l) + mass_bk(l)*ps(ij) ! FIXME : leave ps for the moment ; change ps to Ms later
               !         p_ik = ptop + g*(mass_ak(l)+ mass_bk(l)*ps(i,j))
               exner_ik = cpp * (p_ik/preff) ** kappa
               pk(ij,l) = exner_ik
               ! specific volume v = kappa*theta*pi/p = dphi/g/rhodz
               if ( p_ik > 1.D-30 ) then
                  geopot(ij,l+1) = geopot(ij,l) + (g*kappa)*rhodz(ij,l)*theta(ij,l)*exner_ik/p_ik
               endif
            ENDDO
         ENDDO
         !       ENDIF
      ELSE
         ! We are using a Lagrangian vertical coordinate
         ! Pressure must be computed first top-down (temporarily stored in pk)
         ! Then Exner pressure and geopotential are computed bottom-up
         ! Notice that the computation below should work also when caldyn_eta=eta_mass

         IF(boussinesq) THEN ! compute only geopotential : pressure pk will be computed in compute_caldyn_horiz
            ! specific volume 1 = dphi/g/rhodz
            !         IF (is_omp_level_master) THEN ! no openMP on vertical due to dependency
            DO l = 1,llm
               !DIR$ SIMD
               DO ij=ij_begin_ext,ij_end_ext
                  geopot(ij,l+1) = geopot(ij,l) + g*rhodz(ij,l)
               ENDDO
            ENDDO
         ELSE ! non-Boussinesq, compute geopotential and Exner pressure
            ! uppermost layer

            !DIR$ SIMD
            DO ij=ij_begin_ext,ij_end_ext
               pk(ij,llm) = ptop + (.5*g)*rhodz(ij,llm)
            END DO
            ! other layers
            DO l = llm-1, 1, -1
               !DIR$ SIMD
               DO ij=ij_begin_ext,ij_end_ext
                  pk(ij,l) = pk(ij,l+1) + (.5*g)*(rhodz(ij,l)+rhodz(ij,l+1))
               END DO
            END DO
            ! surface pressure (for diagnostics)
            DO ij=ij_begin_ext,ij_end_ext
               ps(ij) = pk(ij,1) + (.5*g)*rhodz(ij,1)
            END DO

            ! specific volume v = kappa*theta*pi/p = dphi/g/rhodz
            DO l = 1,llm
               !DIR$ SIMD
               DO ij=ij_begin_ext,ij_end_ext
                  p_ik = pk(ij,l)
                  exner_ik = cpp * (p_ik/preff) ** kappa
                  geopot(ij,l+1) = geopot(ij,l) + (g*kappa)*rhodz(ij,l)*theta(ij,l)*exner_ik/p_ik
                  pk(ij,l) = exner_ik
               ENDDO
            ENDDO
         END IF

      END IF

      !ym flush geopot
      !$OMP BARRIER

      CALL trace_end("compute_geopot")

   END SUBROUTINE compute_geopot


END MODULE caldyn_gcm_mod
