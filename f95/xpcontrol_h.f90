module xpcontrol_h

  use const_kind_m
  use clcoef_h

!! public types

!> run control parameters
  type, public :: xpparams_t
     character(len=80) :: control !< option control parameter
     real(kr8) :: realpar !< real control parameter
     integer(ki4) :: intpar !< integer control parameter
     logical :: logicpar !< logical control parameter
  end type xpparams_t

!> file names
  type, public :: xpfiles_t
     character(len=80)  :: out       !< output data
     character(len=80)  :: log           !< log file
     character(len=80)  :: xpcdata         !< xpc input data file
     character(len=80)  :: xpcout         !< xpc output data file
     character(len=80)  :: vtk   !< vtk file
     character(len=80)  :: gnu !< gnuplot file
  end type xpfiles_t

!> plot output selectors
  type, public :: xpplots_t
     logical  :: xpcout !< xpc output data selector
     logical  :: vtk   !< vtk plot selector
     logical  :: gnu !< gnuplot plot selector
  end type xpplots_t

end module xpcontrol_h
