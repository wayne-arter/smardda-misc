module clcoef_h

  use const_kind_m


!> parameters describing how to construct object
  type, public :: ccnumerics_t
     character(len=80) :: formula !< clcoef formula
     real(kr8) :: f !< power split (math variable name allowed)
     integer(ki4) :: nrpams !< number of real parameters
     integer(ki4) :: nipams !< number of integer parameters
     real(kr8), dimension(:), allocatable   :: rpar !< general real parameters
     integer(ki4), dimension(:), allocatable   :: npar !< general integer parameters
  end type ccnumerics_t

  type, public :: clcoef_t
     real(kr8) :: pow !< power
     real(kr8) :: lambda !< Coulomb logarithm used
     real(kr8) :: tau_e !< electron collision time
     real(kr8) :: tau_i !< ion collision time
     real(kr8) :: omega_ce !< electron gyrofrequency
     real(kr8) :: omega_ci !< ion gyrofrequency
     real(kr8) :: m_i !< ion mass (not used)
     real(kr8) :: kappab_epara !< parallel electron thermal transport (Braginskii)
     real(kr8) :: kappab_eperp !< perpendicular electron thermal transport (Braginskii)
     real(kr8) :: kappab_ipara !< parallel ion thermal transport (Braginskii)
     real(kr8) :: kappab_iperp !< perpendicular ion thermal transport (Braginskii)
     real(kr8) :: x_e !< electron gyrofrequency * electron collision time
     real(kr8) :: x_i !< ion gyrofrequency * ion collision time
     real(kr8) :: nu_epara !< parallel electron viscosity
     real(kr8) :: nu_eperp !< perpendicular electron viscosity
     real(kr8) :: nu_ipara !< parallel ion viscosity
     real(kr8) :: nu_iperp !< perpendicular ion viscosity
     real(kr8) :: kappa_epara !< parallel electron thermal diffusivity
     real(kr8) :: kappa_eperp !< perpendicular electron thermal diffusivity
     real(kr8) :: kappa_ipara !< parallel ion thermal diffusivity
     real(kr8) :: kappa_iperp !< perpendicular ion thermal diffusivity
     real(kr8) :: eta !< parallel resistive diffusion
     real(kr8) :: ra !<  Rayleigh number
     real(kr8) :: chandraq !<  Chandrasekhar number
     real(kr8) :: beta !<  plasma beta
     real(kr8) :: lunds !<  Lundquist number

     real(kr8) :: c_tau_e !< numerical constant for electron collision time
     real(kr8) :: c_tau_i !< numerical constant for ion collision time
     real(kr8) :: c_omega_ce !< numerical constant for electron gyrofrequency
     real(kr8) :: c_omega_ci !< numerical constant for ion gyrofrequency
     real(kr8) :: c_kappab_epara !< numerical constant for parallel electron thermal transport (Braginskii)
     real(kr8) :: c_kappab_eperp !< numerical constant for perpendicular electron thermal transport (Braginskii)
     real(kr8) :: c_kappab_ipara !< numerical constant for parallel ion thermal transport (Braginskii)
     real(kr8) :: c_kappab_iperp !< numerical constant for perpendicular ion thermal transport (Braginskii)
     real(kr8) :: c_x_e !< numerical constant for electron gyrofrequency * electron collision time
     real(kr8) :: c_x_i !< numerical constant for ion gyrofrequency * ion collision time
     real(kr8) :: c_nu_epara !< numerical constant for parallel electron viscosity
     real(kr8) :: c_nu_eperp !< numerical constant for perpendicular electron viscosity
     real(kr8) :: c_nu_ipara !< numerical constant for parallel ion viscosity
     real(kr8) :: c_nu_iperp !< numerical constant for perpendicular ion viscosity
     real(kr8) :: c_kappa_epara !< numerical constant for parallel electron thermal diffusivity
     real(kr8) :: c_kappa_eperp !< numerical constant for perpendicular electron thermal diffusivity
     real(kr8) :: c_kappa_ipara !< numerical constant for parallel ion thermal diffusivity
     real(kr8) :: c_kappa_iperp !< numerical constant for perpendicular ion thermal diffusivity
     real(kr8) :: c_eta !< numerical constant for parallel resistive diffusion
     real(kr8) :: c_ra !< numerical constant for Rayleigh number
     real(kr8) :: c_chandraq !< numerical constant for Chandrasekhar number
     real(kr8) :: c_beta !< numerical constant for plasma beta
     real(kr8) :: c_lunds !< numerical constant for Lundquist number
     type(ccnumerics_t) :: n !< control  parameters
  end type clcoef_t

end module clcoef_h
