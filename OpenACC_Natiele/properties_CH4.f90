subroutine properties
    use comum
    implicit none
    REAL(8)        :: cp_air, rho_air, k_air, nu_air, mu_air, alpha_air
    REAL(8)        :: cp_met, rho_met, k_met, mu_met, nu_met, alpha_met
    REAL(8), parameter     :: X_met=0.25d0, X_air=0.75d0

    rho_air = 356.6123426031d0 * Too ** (- 1.0013371777d0)

    cp_air = 4.73251788718275d-11 * Too ** 4.d0 - 2.5966626136914d-7 * &
             Too ** 3.d0 + 0.0004309589d0 * Too ** 2.d0 - 0.0748376208d0 * Too + 993.1365094124d0

    k_air = 4.63985599705762d-12 * Too ** 3.d0 - 0.000000029d0 * &
            Too ** 2.d0 + 9.04946015844389d-5 * Too + 0.0009128206d0

    alpha_air =  - 1.50263396909776d-14 * Too ** 3.d0 + 1.0946104054489d-10 * Too**2.d0 &
                 + 8.25165965131297d-8 * Too - 1.20726065874961d-5

    mu_air = 4.40778022994337d-15 * Too ** 3.d0 - 2.31633149512214d-11 * Too ** 2.d0 &
             + 5.78517134041232d-8 * Too + 2.98966575978449d-6

    nu_air = mu_air / rho_air

!!!!NIST webbook for Methane (C H 4 ) at 1atm

    cp_met = (3.42925052210868d-6*Too**2.d0 + 0.0003463592d0*Too + 1.8482748861d0) * 10**3.d0![J/kgK]

    rho_met =  2.60191696337177d-11*Too**4.d0 - 5.26952306479639d-8*Too**3.d0 + &
               4.12957366234466d-5*Too**2.d0 - 0.0155954244d0*Too + 2.8262044923d0 ![kg/m3]

    k_met =  1.30126895841081d-7*Too**2.d0 + 0.000064197d0*Too + 0.0037188796d0 ![W/mK]

    mu_met = 2.84900946271051d-8 * Too + 2.63131387329591d-6 ![Ns/m2 ]
        
    
    nu_met = mu_met/rho_met 

    alpha_met = (k_met)/(cp_met*rho_met)


    rho_tot = X_met * rho_met + X_air * rho_air

    cp_tot  = X_met * cp_met +  X_air * cp_air

    k_tot   = X_met * k_met +   X_air * k_air
! 
    nu_tot  = X_met * nu_met +  X_air * nu_air

    alpha_tot  = X_met * alpha_met +  X_air * alpha_air

end subroutine properties

!################# OLD ##################
!subroutine properties
!    use comum
!    implicit none
!    REAL(8)        :: cp_met, rho_met, k_met, mu_met, nu_met, alpha_met
!    REAL(8)        :: cp_N2, rho_N2, k_N2, nu_N2, alpha_N2
!    REAL(8)        :: cp_O2, rho_O2, k_O2, nu_O2, alpha_O2
!    REAL(8), parameter     :: X_met=0.25d0, X_N2=0.50d0, X_O2 = 0.25d0

!!    !Stephen R Turns et al. An introduction to combustion 1996. Table B.3

!!    cp_met = (-7.40308d2 + 1.0893537d1*Too - 1.265124d-2*Too**2.d0 + 9.843763d-6*Too**3.d0 &
!!            - 4.3228296d-9*Too**4.d0 + 7.863665d-13*Too**5.d0 )![J/kgK]

!!    rho_met = 5.10490319d-6 *Too**3.d0 - 0.0045335654*Too**2.d0 + 1.3537197365d0*Too - 135.5452328229d0 ![kg/m3]


!!    k_met = - 4.606147d-2 + 5.95652224d-4*Too - 2.98893153d-6*Too**2.d0 + 8.44612876d-9*Too**3.d0 &
!!            - 1.22927d-11*Too**4.d0 + 9.0127d-15*Too**5.d0 - 2.62961d-18*Too**6.d0![W/mK]


!!    mu_met = (1.540087d0 + 1.095157d-2*Too + 1.800664d-5*Too**2.d0 - 1.36379d-8*Too**3.d0 )*10**6.d0 ![Ns/m2 ]
!!    

!!    nu_met = mu_met/rho_met 

!    !NIST webbook for n-mettane (C 7 H 16 ) at 1atm

!!    cp_met = (0.0037638842d0*Too + 0.6385320938d0) * 10**3.d0![J/kgK]

!!    rho_met = 0.000015878d0*Too**2.d0 - 0.0214227788d0*Too + 9.2184793364d0 ![kg/m3]

!!    k_met = 0.0001099183d0*Too - 0.02319995272d0 ![W/mK]

!!    mu_met = 1.88414068667895d-8*Too + 6.16915491028679d-7 ![Ns/m2 ]
!!    
!!    nu_met = mu_met/rho_met 

!!    alpha_met = (k_met)/(cp_met*rho_met)


!    cp_air = 1167.d0 ![J/kgK]

!    rho_N2 = 0.3368d0 ![kg/m3]

!    k_N2 = 64.7d-3 ![W/mK]

!    nu_N2 = 118.7d-6 ![m2/s]

!    alpha_N2 = (k_N2)/(cp_N2*rho_N2)


!    cp_O2 = 1090.d0 ![J/kgK]

!    rho_O2 = 0.3848d0 ![kg/m3]

!    k_O2 = 71.d-3 ![W/mK]

!    nu_O2 = 124.d-6 ![m2/s]

!    alpha_O2 = (k_O2)/(cp_O2*rho_O2)



!    rho_tot = X_met * rho_met + X_N2 * rho_N2 + X_O2 * rho_O2

!    cp_tot  = X_met * cp_met +  X_N2 * cp_N2 +  X_O2 * cp_O2

!    k_tot   = X_met * k_met +   X_N2 * k_N2 +   X_O2 * k_O2
! 
!    nu_tot  = X_met * nu_met +  X_N2 * nu_N2 +  X_O2 * nu_O2

!    alpha_tot  = X_met * alpha_met +  X_N2 * alpha_N2 +  X_O2 * alpha_O2



!end subroutine properties

