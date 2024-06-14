&progfiles
/
&miscparameters
/
&plotselections
/
&xpcparameters
 ! These are used to define Braginskii coefficients
 atomic_number = 1.0,
 effective_charge = 1.0,
 imposed_Bfield = 1.0,
 electron_temperature = 10.,
 ion_temperature = 10.,
 electron_number_density = 1.0E+18,
 !! Set this negative and electron-ion formula will be used
 Coulomb_logarithm = 14.,
 ! Usage to define numerical factors in formulae as functions
 ! of n, Te and Lambda (Te factor always unity implying eV)
 number_density_factor = 1.0E+18,
 Coulomb_logarithm_factor = 1.0,
 ! Usage in convection type formulae
 layer_depth = 0.01,
 pressure_lengthscale = 0.01,
 minor_radius = 0.1
 effective_Bfield = 1.0,
 !! Determine effective gravity
 poloidal_angle = 0.0,
 major_radius = 1.0,
 ! For future expansion
 effective_field_formula="null",
 power_split = 0.50,
 xpc_formula="unset",
 general_real_parameters = 10*0.0,
 number_of_real_parameters = 0,
 general_integer_parameters = 10*0,
 number_of_integer_parameters = 0,
 /
&clcoefparameters
/
