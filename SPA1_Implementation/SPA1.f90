! Calculates the Saddle Point Approximation 1st order term for 1D toy model (G. Granucci, M. Persico and G. Spighi, J. Chem. Phys., 2012, 137, 22A501)
! Analytical Computation only valid for specific model only
! SPA 1st order term:
!    Interstate:
!        dSO_{I,J}^{Msi,Msj}(X_c)/dx * < X_i^I | x-x_c | X_j^I >
!
!        (r_sigma-1 < x_c < r_sigma +1) :
!                  dSO_{I,J}^{Msi,Msj}(X_c) / dx = 12 * ( (x_c-r_sigma)^2) / dr_sigma^3 ) - 3/dr_sigma
!         otherwise:
!                  dSO_{I,J}^{Msi,Msj}(X_c) / dx = 0
!
!    Intrastate:
!        dV(X_c)/dx * < X_i^I | x-x_c | X_j^I >
!
!        Singlet state:
!                  dV(X_c)/dx =  -a_1 * alpha_1 * e^(-alpha_1 * x_c)
!        Triplet state:
!                  dV(X_c)/dx =  -a_2 * alpha_2 * e^(-alpha_2 * x_c)
!
!
!  < X_i^I | x-x_c | X_j^I > = -i * (p_1-p_2)/4 * alpha_I * < X_i^I | X_j^I > (M. Šulc, H. Hernández, T. J. Martínez, J. Vaníček, J. Chem. Phys., 2013, 139, 034112)

   function SPA1_SOC_model(T1, T2, S_ij) result(PotEn)
      type(T_Trajectory), intent(in) :: T1, T2
      complex(kind=DefComp), intent(in) :: S_ij
      real(kind=DefReal) :: x_cent, x_1, x_2, p_1, p_2, sigma_G, alpha_1, alpha_2
      complex(kind=DefComp) :: PotEn, dE, dSOC, roe

      if (glIMethod /= 4) then
         call FMS_DieError("SPA1 only available for GAIMS 1D Toy Model")
      end if

      x_1 = T1%Particle(1)%get_pos(1)
      x_2 = T2%Particle(1)%get_pos(1)
      p_1 = T1%Particle(1)%get_mom(1)
      p_2 = T2%Particle(1)%get_mom(1)
      alpha_1 = T1%Particle(1)%width
      alpha_2 = T2%Particle(1)%width

      x_cent = (alpha_1 * x_1 + alpha_2 * x_2) / (alpha_1 + alpha_2) !centroid position
      ! roe from appendix A of M. Šulc, H. Hernández, T. J. Martínez, J. Vaníček, J. Chem. Phys., 2013, 139, 034112
      roe = -c1i * (p_1 - p_2) / (4 * (alpha_1))

      if ((T1%StateID == T2%StateID) .and. (T1%Ms == T2%Ms)) then !for intrastate case A
         !Derivative of potential energy
         if (T1%StateID == 2) then !T
            dE = -0.25d0 * 0.5d0 * exp(-0.25d0 * x_cent) !Derivative of V from G. Granucci, M. Persico and G. Spighi, J. Chem. Phys., 2012, 137, 22A501
         else !S
            dE = -0.35d0 * 0.03452d0 * exp(-0.35d0 * x_cent)
         end if
         PotEn = dE * S_ij * roe
      elseif ((T1%StateID == T2%StateID) .and. (T1%Ms /= T2%Ms)) then !for interstate case F
         dSOC = (0.0, 0.0)
         PotEn = dSOC * S_ij * roe
      else !for interstate case D
         !Derivative of SO_{I,J}^{Msi,Msj}(X_c) in Eq(11) from G. Granucci, M. Persico and G. Spighi, J. Chem. Phys., 2012, 137, 22A501
         if (((glGrsigma - 1.d0) < x_cent) .and. (x_cent < (glGrsigma + 1.d0))) then
            sigma_G = 12.d0 * (((x_cent - glGrsigma)**2) / (2.d0**3)) - 3.d0 / 2.d0
         else
            sigma_G = 0.d0
         end if

         if (T1%StateID == 1) then
            if (T2%Ms == 1) then
               dSOC = conjg((0.0005, 0.0005) * sigma_G) ! z*  : S,T-1
            elseif (T2%Ms == 2) then
               dSOC = c1i * 0.001d0 * sigma_G ! ib  : S,T0
            else
               dSOC = (0.0005, 0.0005) * sigma_G ! z   : S,T1
            end if
         else
            if (T1%Ms == 1) then
               dSOC = (0.0005, 0.0005) * sigma_G ! z   : T-1,S
            elseif (T1%Ms == 2) then
               dSOC = -c1i * 0.001d0 * sigma_G ! -ib : T0,S
            else
               dSOC = conjg((0.0005, 0.0005) * sigma_G) ! z*  : T1,S
            end if
         end if
         PotEn = dSOC * S_ij * roe
      end if

   end function SPA1_SOC_model
