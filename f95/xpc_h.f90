module xpc_h

  use const_kind_m
  use clcoef_h


!> parameters describing how to construct object
  type, public :: xpnumerics_t
     real(kr8) :: a !< atomic_number
     real(kr8) :: z !< effective_charge
     real(kr8) :: b !< imposed_Bfield
     real(kr8) :: t_e !< electron_temperature
     real(kr8) :: t_i !< ion_temperature
     real(kr8) :: c_n !< number_density_factor
     real(kr8) :: n !< electron_number_density
     real(kr8) :: c_lambda !< Coulomb_logarithm_factor
     real(kr8) :: lambda !< Coulomb_logarithm
     character(len=80) :: formula !< xpc formula
     real(kr8) :: f !< power split (math variable name allowed)
     integer(ki4) :: nrpams !< number of real parameters
     integer(ki4) :: nipams !< number of integer parameters
     real(kr8), dimension(:), allocatable   :: rpar !< general real parameters
     integer(ki4), dimension(:), allocatable   :: npar !< general integer parameters
  end type xpnumerics_t

  type, public :: xpc_t
     real(kr8) :: pow !< power
     type(clcoef_t) :: clcoef !< classical coefficients
     type(xpnumerics_t) :: n !< control  parameters
  end type xpc_t

end module xpc_h
