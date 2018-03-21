MODULE caldyn_gcm_mod
   use mod_misc
   PRIVATE
   TYPE(t_message) :: req_ps, req_mass, req_theta_rhodz, req_u, req_qu

   public compute_caldyn_horiz

CONTAINS

   SUBROUTINE compute_caldyn_horiz(u,rhodz,qu,theta,pk,geopot, hflux,convm, dtheta_rhodz, du)
      use prec
      use mod_misc
      IMPLICIT NONE
      REAL(rstd),INTENT(IN)  :: u(iim*3*jjm,llm)    ! prognostic "velocity"
      REAL(rstd),INTENT(IN)  :: rhodz(iim*jjm,llm)
      REAL(rstd),INTENT(IN)  :: qu(iim*3*jjm,llm)
      REAL(rstd),INTENT(IN)  :: theta(iim*jjm,llm)  ! potential temperature
      REAL(rstd),INTENT(INOUT) :: pk(iim*jjm,llm) ! Exner function
      REAL(rstd),INTENT(IN)  :: geopot(iim*jjm,llm+1)    ! geopotential

      REAL(rstd),INTENT(INOUT) :: hflux(iim*3*jjm,llm) ! hflux in kg/s
      REAL(rstd),INTENT(INOUT) :: convm(iim*jjm,llm)  ! mass flux convergence
      REAL(rstd),INTENT(INOUT) :: dtheta_rhodz(iim*jjm,llm)
      REAL(rstd),INTENT(INOUT) :: du(iim*3*jjm,llm)

      REAL(rstd) :: cor_NT(iim*jjm,llm)  ! NT coriolis force u.(du/dPhi)
      REAL(rstd) :: urel(3*iim*jjm,llm)  ! relative velocity
      REAL(rstd) :: Ftheta(3*iim*jjm,llm) ! theta flux
      REAL(rstd) :: berni(iim*jjm,llm)  ! Bernoulli function

      INTEGER :: i,j,ij,l
      REAL(rstd) :: ww,uu

      CALL trace_start("compute_caldyn_horiz")

      !    CALL wait_message(req_theta_rhodz) 

      DO l = ll_begin, ll_end
!!!  Compute mass and theta fluxes
         IF (caldyn_conserv==energy) CALL test_message(req_qu) 
         !DIR$ SIMD
         DO ij=ij_begin_ext,ij_end_ext
            hflux(ij+u_right,l)=0.5*(rhodz(ij,l)+rhodz(ij+t_right,l))*u(ij+u_right,l)*le(ij+u_right)
            hflux(ij+u_lup,l)=0.5*(rhodz(ij,l)+rhodz(ij+t_lup,l))*u(ij+u_lup,l)*le(ij+u_lup)
            hflux(ij+u_ldown,l)=0.5*(rhodz(ij,l)+rhodz(ij+t_ldown,l))*u(ij+u_ldown,l)*le(ij+u_ldown)

            Ftheta(ij+u_right,l)=0.5*(theta(ij,l)+theta(ij+t_right,l))*hflux(ij+u_right,l)
            Ftheta(ij+u_lup,l)=0.5*(theta(ij,l)+theta(ij+t_lup,l))*hflux(ij+u_lup,l)
            Ftheta(ij+u_ldown,l)=0.5*(theta(ij,l)+theta(ij+t_ldown,l))*hflux(ij+u_ldown,l)
         ENDDO

!!! compute horizontal divergence of fluxes
         !DIR$ SIMD
         DO ij=ij_begin,ij_end         
            ! convm = -div(mass flux), sign convention as in Ringler et al. 2012, eq. 21
            convm(ij,l)= -1./Ai(ij)*(ne_right*hflux(ij+u_right,l)   +  &
                 ne_rup*hflux(ij+u_rup,l)       +  &  
                 ne_lup*hflux(ij+u_lup,l)       +  &  
                 ne_left*hflux(ij+u_left,l)     +  &  
                 ne_ldown*hflux(ij+u_ldown,l)   +  &  
                 ne_rdown*hflux(ij+u_rdown,l))    

            ! signe ? attention d (rho theta dz)
            ! dtheta_rhodz =  -div(flux.theta)
            dtheta_rhodz(ij,l)=-1./Ai(ij)*(ne_right*Ftheta(ij+u_right,l)    +  &
                 ne_rup*Ftheta(ij+u_rup,l)        +  &  
                 ne_lup*Ftheta(ij+u_lup,l)        +  &  
                 ne_left*Ftheta(ij+u_left,l)      +  &  
                 ne_ldown*Ftheta(ij+u_ldown,l)    +  &  
                 ne_rdown*Ftheta(ij+u_rdown,l))
         ENDDO

      END DO

!!! Compute potential vorticity (Coriolis) contribution to du

      SELECT CASE(caldyn_conserv)
      CASE(energy) ! energy-conserving TRiSK

         CALL wait_message(req_qu)

         DO l=ll_begin,ll_end
            !DIR$ SIMD
            DO ij=ij_begin,ij_end         

               !             if ( de(ij+u_right) > 1.0 ) then
               uu = wee(ij+u_right,1,1)*hflux(ij+u_rup,l)*(qu(ij+u_right,l)+qu(ij+u_rup,l))+                             &
                    wee(ij+u_right,2,1)*hflux(ij+u_lup,l)*(qu(ij+u_right,l)+qu(ij+u_lup,l))+                             &
                    wee(ij+u_right,3,1)*hflux(ij+u_left,l)*(qu(ij+u_right,l)+qu(ij+u_left,l))+                           &
                    wee(ij+u_right,4,1)*hflux(ij+u_ldown,l)*(qu(ij+u_right,l)+qu(ij+u_ldown,l))+                         &
                    wee(ij+u_right,5,1)*hflux(ij+u_rdown,l)*(qu(ij+u_right,l)+qu(ij+u_rdown,l))+                         & 
                    wee(ij+u_right,1,2)*hflux(ij+t_right+u_ldown,l)*(qu(ij+u_right,l)+qu(ij+t_right+u_ldown,l))+         &
                    wee(ij+u_right,2,2)*hflux(ij+t_right+u_rdown,l)*(qu(ij+u_right,l)+qu(ij+t_right+u_rdown,l))+         &
                    wee(ij+u_right,3,2)*hflux(ij+t_right+u_right,l)*(qu(ij+u_right,l)+qu(ij+t_right+u_right,l))+         &
                    wee(ij+u_right,4,2)*hflux(ij+t_right+u_rup,l)*(qu(ij+u_right,l)+qu(ij+t_right+u_rup,l))+             &
                    wee(ij+u_right,5,2)*hflux(ij+t_right+u_lup,l)*(qu(ij+u_right,l)+qu(ij+t_right+u_lup,l))
               du(ij+u_right,l) = .5*uu/de(ij+u_right)
               !             endif

               if ( de(ij+u_lup) > 1.0 ) then
                  uu = wee(ij+u_lup,1,1)*hflux(ij+u_left,l)*(qu(ij+u_lup,l)+qu(ij+u_left,l)) +                        &
                       wee(ij+u_lup,2,1)*hflux(ij+u_ldown,l)*(qu(ij+u_lup,l)+qu(ij+u_ldown,l)) +                      &
                       wee(ij+u_lup,3,1)*hflux(ij+u_rdown,l)*(qu(ij+u_lup,l)+qu(ij+u_rdown,l)) +                      &
                       wee(ij+u_lup,4,1)*hflux(ij+u_right,l)*(qu(ij+u_lup,l)+qu(ij+u_right,l)) +                      &
                       wee(ij+u_lup,5,1)*hflux(ij+u_rup,l)*(qu(ij+u_lup,l)+qu(ij+u_rup,l)) +                          & 
                       wee(ij+u_lup,1,2)*hflux(ij+t_lup+u_right,l)*(qu(ij+u_lup,l)+qu(ij+t_lup+u_right,l)) +          &
                       wee(ij+u_lup,2,2)*hflux(ij+t_lup+u_rup,l)*(qu(ij+u_lup,l)+qu(ij+t_lup+u_rup,l)) +              &
                       wee(ij+u_lup,3,2)*hflux(ij+t_lup+u_lup,l)*(qu(ij+u_lup,l)+qu(ij+t_lup+u_lup,l)) +              &
                       wee(ij+u_lup,4,2)*hflux(ij+t_lup+u_left,l)*(qu(ij+u_lup,l)+qu(ij+t_lup+u_left,l)) +            &
                       wee(ij+u_lup,5,2)*hflux(ij+t_lup+u_ldown,l)*(qu(ij+u_lup,l)+qu(ij+t_lup+u_ldown,l))
                  du(ij+u_lup,l) = .5*uu/de(ij+u_lup)
               endif

               if ( de(ij+u_ldown) > 1.0 ) then
                  uu = wee(ij+u_ldown,1,1)*hflux(ij+u_rdown,l)*(qu(ij+u_ldown,l)+qu(ij+u_rdown,l)) +                      &
                       wee(ij+u_ldown,2,1)*hflux(ij+u_right,l)*(qu(ij+u_ldown,l)+qu(ij+u_right,l)) +                      &
                       wee(ij+u_ldown,3,1)*hflux(ij+u_rup,l)*(qu(ij+u_ldown,l)+qu(ij+u_rup,l)) +                          &
                       wee(ij+u_ldown,4,1)*hflux(ij+u_lup,l)*(qu(ij+u_ldown,l)+qu(ij+u_lup,l)) +                          &
                       wee(ij+u_ldown,5,1)*hflux(ij+u_left,l)*(qu(ij+u_ldown,l)+qu(ij+u_left,l)) +                        & 
                       wee(ij+u_ldown,1,2)*hflux(ij+t_ldown+u_lup,l)*(qu(ij+u_ldown,l)+qu(ij+t_ldown+u_lup,l)) +          &
                       wee(ij+u_ldown,2,2)*hflux(ij+t_ldown+u_left,l)*(qu(ij+u_ldown,l)+qu(ij+t_ldown+u_left,l)) +        &
                       wee(ij+u_ldown,3,2)*hflux(ij+t_ldown+u_ldown,l)*(qu(ij+u_ldown,l)+qu(ij+t_ldown+u_ldown,l)) +      &
                       wee(ij+u_ldown,4,2)*hflux(ij+t_ldown+u_rdown,l)*(qu(ij+u_ldown,l)+qu(ij+t_ldown+u_rdown,l)) +      &
                       wee(ij+u_ldown,5,2)*hflux(ij+t_ldown+u_right,l)*(qu(ij+u_ldown,l)+qu(ij+t_ldown+u_right,l))
                  du(ij+u_ldown,l) = .5*uu/de(ij+u_ldown)
               endif

            ENDDO
         ENDDO

      CASE(enstrophy) ! enstrophy-conserving TRiSK

         DO l=ll_begin,ll_end
            !DIR$ SIMD
            DO ij=ij_begin,ij_end         

               !             if ( de(ij+u_right) > 1.0 ) then
               uu = wee(ij+u_right,1,1)*hflux(ij+u_rup,l)+                           &
                    wee(ij+u_right,2,1)*hflux(ij+u_lup,l)+                           &
                    wee(ij+u_right,3,1)*hflux(ij+u_left,l)+                          &
                    wee(ij+u_right,4,1)*hflux(ij+u_ldown,l)+                         &
                    wee(ij+u_right,5,1)*hflux(ij+u_rdown,l)+                         & 
                    wee(ij+u_right,1,2)*hflux(ij+t_right+u_ldown,l)+                 &
                    wee(ij+u_right,2,2)*hflux(ij+t_right+u_rdown,l)+                 &
                    wee(ij+u_right,3,2)*hflux(ij+t_right+u_right,l)+                 &
                    wee(ij+u_right,4,2)*hflux(ij+t_right+u_rup,l)+                   &
                    wee(ij+u_right,5,2)*hflux(ij+t_right+u_lup,l)
               du(ij+u_right,l) = qu(ij+u_right,l)*uu/de(ij+u_right)
               !             endif

               if ( de(ij+u_lup) > 1.0 ) then
                  uu = wee(ij+u_lup,1,1)*hflux(ij+u_left,l)+                        &
                       wee(ij+u_lup,2,1)*hflux(ij+u_ldown,l)+                       &
                       wee(ij+u_lup,3,1)*hflux(ij+u_rdown,l)+                       &
                       wee(ij+u_lup,4,1)*hflux(ij+u_right,l)+                       &
                       wee(ij+u_lup,5,1)*hflux(ij+u_rup,l)+                         & 
                       wee(ij+u_lup,1,2)*hflux(ij+t_lup+u_right,l)+                 &
                       wee(ij+u_lup,2,2)*hflux(ij+t_lup+u_rup,l)+                   &
                       wee(ij+u_lup,3,2)*hflux(ij+t_lup+u_lup,l)+                   &
                       wee(ij+u_lup,4,2)*hflux(ij+t_lup+u_left,l)+                  &
                       wee(ij+u_lup,5,2)*hflux(ij+t_lup+u_ldown,l)
                  du(ij+u_lup,l) = qu(ij+u_lup,l)*uu/de(ij+u_lup)
               endif

               if ( de(ij+u_ldown) > 1.0 ) then
                  uu = wee(ij+u_ldown,1,1)*hflux(ij+u_rdown,l)+                      &
                       wee(ij+u_ldown,2,1)*hflux(ij+u_right,l)+                      &
                       wee(ij+u_ldown,3,1)*hflux(ij+u_rup,l)+                        &
                       wee(ij+u_ldown,4,1)*hflux(ij+u_lup,l)+                        &
                       wee(ij+u_ldown,5,1)*hflux(ij+u_left,l)+                       & 
                       wee(ij+u_ldown,1,2)*hflux(ij+t_ldown+u_lup,l)+                &
                       wee(ij+u_ldown,2,2)*hflux(ij+t_ldown+u_left,l)+               &
                       wee(ij+u_ldown,3,2)*hflux(ij+t_ldown+u_ldown,l)+              &
                       wee(ij+u_ldown,4,2)*hflux(ij+t_ldown+u_rdown,l)+              &
                       wee(ij+u_ldown,5,2)*hflux(ij+t_ldown+u_right,l)
                  du(ij+u_ldown,l) = qu(ij+u_ldown,l)*uu/de(ij+u_ldown)
               endif

            ENDDO
         ENDDO

      CASE DEFAULT
         STOP
      END SELECT

!!! Compute bernouilli term = Kinetic Energy + geopotential
      IF(boussinesq) THEN
         ! first use hydrostatic balance with theta*rhodz to find pk (Lagrange multiplier=pressure) 
         ! uppermost layer
         !DIR$ SIMD
         DO ij=ij_begin_ext,ij_end_ext         
            pk(ij,llm) = ptop + (.5*g)*theta(ij,llm)*rhodz(ij,llm)
         END DO
         ! other layers
         DO l = llm-1, 1, -1
            !          !$OMP DO SCHEDULE(STATIC) 
            !DIR$ SIMD
            DO ij=ij_begin_ext,ij_end_ext         
               pk(ij,l) = pk(ij,l+1) + (.5*g)*(theta(ij,l)*rhodz(ij,l)+theta(ij,l+1)*rhodz(ij,l+1))
            END DO
         END DO
         ! surface pressure (for diagnostics) FIXME
         ! DO ij=ij_begin_ext,ij_end_ext         
         !   ps(ij) = pk(ij,1) + (.5*g)*theta(ij,1)*rhodz(ij,1)
         ! END DO
         ! now pk contains the Lagrange multiplier (pressure)


         DO l=ll_begin,ll_end
            DO ij=ij_begin_ext,ij_end_ext
               berni(ij,l) = pk(ij,l)
            END DO
            !DIR$ SIMD
            DO ij=ij_begin,ij_end
               if ( Ai(ij) > 1.e-30 ) then
                  berni(ij,l) = berni(ij,l) + &
                       1/(4*Ai(ij))*(le(ij+u_right)*de(ij+u_right)*u(ij+u_right,l)**2 +    &
                       le(ij+u_rup)*de(ij+u_rup)*u(ij+u_rup,l)**2 +          &
                       le(ij+u_lup)*de(ij+u_lup)*u(ij+u_lup,l)**2 +          &
                       le(ij+u_left)*de(ij+u_left)*u(ij+u_left,l)**2 +       &
                       le(ij+u_ldown)*de(ij+u_ldown)*u(ij+u_ldown,l)**2 +    &
                       le(ij+u_rdown)*de(ij+u_rdown)*u(ij+u_rdown,l)**2 )  
               endif
            ENDDO

            ! from now on pk contains the vertically-averaged geopotential
            DO ij=ij_begin_ext,ij_end_ext
               pk(ij,l) = .5*(geopot(ij,l)+geopot(ij,l+1))
            END DO
         ENDDO

      ELSE ! compressible

         berni(:,:) = 0.0

         DO l=ll_begin,ll_end
            DO ij=ij_begin_ext,ij_end_ext
               berni(ij,l) = .5*(geopot(ij,l)+geopot(ij,l+1))
            END DO
            !DIR$ SIMD
            DO ij=ij_begin,ij_end
               if ( Ai(ij) > 1.e-30 ) then
                  berni(ij,l) = berni(ij,l) &
                       + 1/(4*Ai(ij))*(le(ij+u_right)*de(ij+u_right)*u(ij+u_right,l)**2 +    &
                       le(ij+u_rup)*de(ij+u_rup)*u(ij+u_rup,l)**2 +          &
                       le(ij+u_lup)*de(ij+u_lup)*u(ij+u_lup,l)**2 +          &
                       le(ij+u_left)*de(ij+u_left)*u(ij+u_left,l)**2 +       &
                       le(ij+u_ldown)*de(ij+u_ldown)*u(ij+u_ldown,l)**2 +    &
                       le(ij+u_rdown)*de(ij+u_rdown)*u(ij+u_rdown,l)**2 )  
               endif
            ENDDO
         ENDDO

      END IF ! Boussinesq/compressible

!!! Add gradients of Bernoulli and Exner functions to du
      DO l=ll_begin,ll_end
         !DIR$ SIMD
         DO ij=ij_begin,ij_end

            if ( de(ij+u_right) > 1.0 ) then
               du(ij+u_right,l) = du(ij+u_right,l) + 1/de(ij+u_right) * (             &
                    0.5*(theta(ij,l)+theta(ij+t_right,l))                &
                    *( ne_right*pk(ij,l)+ne_left*pk(ij+t_right,l))    &
                    + ne_right*berni(ij,l)+ne_left*berni(ij+t_right,l) )
            endif

            if ( de(ij+u_lup) > 1.0 ) then
               du(ij+u_lup,l) = du(ij+u_lup,l) +  1/de(ij+u_lup) * (                  &
                    0.5*(theta(ij,l)+theta(ij+t_lup,l))                  &
                    *( ne_lup*pk(ij,l)+ne_rdown*pk(ij+t_lup,l))       &
                    + ne_lup*berni(ij,l)+ne_rdown*berni(ij+t_lup,l) )
            endif

            if ( de(ij+u_ldown) > 1.0 ) then
               du(ij+u_ldown,l) = du(ij+u_ldown,l) + 1/de(ij+u_ldown) * (             &
                    0.5*(theta(ij,l)+theta(ij+t_ldown,l))                &
                    *( ne_ldown*pk(ij,l)+ne_rup*pk(ij+t_ldown,l))     &
                    + ne_ldown*berni(ij,l)+ne_rup*berni(ij+t_ldown,l) )
            endif

         ENDDO
      ENDDO

      CALL trace_end("compute_caldyn_horiz")

   END SUBROUTINE compute_caldyn_horiz
END MODULE caldyn_gcm_mod
