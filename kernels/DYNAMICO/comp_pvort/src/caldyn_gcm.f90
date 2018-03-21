MODULE caldyn_gcm_mod
   use mod_misc
   PRIVATE
   TYPE(t_message) :: req_ps, req_mass, req_theta_rhodz, req_u, req_qu ! dummy here

   public compute_pvort

CONTAINS

   SUBROUTINE compute_pvort(ps,u,theta_rhodz, rhodz,theta,qu,qv)
      use prec
      use mod_misc
      IMPLICIT NONE
      REAL(rstd),INTENT(IN)  :: u(iim*3*jjm,llm)
      REAL(rstd),INTENT(IN)  :: ps(iim*jjm)
      REAL(rstd),INTENT(IN)  :: theta_rhodz(iim*jjm,llm)
      REAL(rstd),INTENT(INOUT) :: rhodz(iim*jjm,llm)
      REAL(rstd),INTENT(INOUT) :: theta(iim*jjm,llm)
      REAL(rstd),INTENT(INOUT) :: qu(iim*3*jjm,llm)
      REAL(rstd),INTENT(INOUT) :: qv(iim*2*jjm,llm)

      INTEGER :: i,j,ij,l
      REAL(rstd) :: etav,hv, m

      CALL trace_start("compute_pvort")

      IF(caldyn_eta==eta_mass) THEN
         CALL wait_message(req_ps)
      ELSE
         CALL wait_message(req_mass)
      END IF
      CALL wait_message(req_theta_rhodz)

      IF(caldyn_eta==eta_mass) THEN ! Compute mass & theta
         DO l = ll_begin,ll_end
            CALL test_message(req_u)
            !DIR$ SIMD
            DO ij=ij_begin_ext,ij_end_ext
               m = ( mass_dak(l)+ps(ij)*mass_dbk(l) )/g
               rhodz(ij,l) = m
               if ( rhodz(ij,l) > 0.0 ) then
                  theta(ij,l) = theta_rhodz(ij,l)/rhodz(ij,l)
               else
                  theta(ij,l) = 0.0
               endif
            ENDDO
         ENDDO
      ELSE ! Compute only theta
         DO l = ll_begin,ll_end
            CALL test_message(req_u)
            !DIR$ SIMD
            DO ij=ij_begin_ext,ij_end_ext
               theta(ij,l) = theta_rhodz(ij,l)/rhodz(ij,l)
            ENDDO
         ENDDO
      END IF

      CALL wait_message(req_u)

!!! Compute shallow-water potential vorticity
      DO l = ll_begin,ll_end
         !DIR$ SIMD
         DO ij=ij_begin_ext,ij_end_ext
            etav= 1./Av(ij+z_up)*(  ne_rup        * u(ij+u_rup,l)        * de(ij+u_rup)         &
                 + ne_left * u(ij+t_rup+u_left,l) * de(ij+t_rup+u_left)  &
                 - ne_lup        * u(ij+u_lup,l)        * de(ij+u_lup) )

            hv =  Riv2(ij,vup)          * rhodz(ij,l)            &
                 + Riv2(ij+t_rup,vldown) * rhodz(ij+t_rup,l)     &
                 + Riv2(ij+t_lup,vrdown) * rhodz(ij+t_lup,l)

            if ( hv > 0.0 ) then
               qv(ij+z_up,l) = ( etav+fv(ij+z_up) )/hv
            else
               qv(ij+z_up,l) = 0.0
            endif

            etav = 1./Av(ij+z_down)*(  ne_ldown         * u(ij+u_ldown,l)          * de(ij+u_ldown)          &
                 + ne_right * u(ij+t_ldown+u_right,l)  * de(ij+t_ldown+u_right)  &
                 - ne_rdown         * u(ij+u_rdown,l)          * de(ij+u_rdown) )

            hv =  Riv2(ij,vdown)        * rhodz(ij,l)          &
                 + Riv2(ij+t_ldown,vrup) * rhodz(ij+t_ldown,l)  &
                 + Riv2(ij+t_rdown,vlup) * rhodz(ij+t_rdown,l)

            if ( hv > 0.0 ) then
               qv(ij+z_down,l) =( etav+fv(ij+z_down) )/hv
            else
               qv(ij+z_down,l) = 0.0
            endif

         ENDDO

         !DIR$ SIMD
         DO ij=ij_begin,ij_end
            qu(ij+u_right,l) = 0.5*(qv(ij+z_rdown,l)+qv(ij+z_rup,l))
            qu(ij+u_lup,l) = 0.5*(qv(ij+z_up,l)+qv(ij+z_lup,l))
            qu(ij+u_ldown,l) = 0.5*(qv(ij+z_ldown,l)+qv(ij+z_down,l))
         END DO

      ENDDO

      CALL trace_end("compute_pvort")

   END SUBROUTINE compute_pvort

END MODULE caldyn_gcm_mod
