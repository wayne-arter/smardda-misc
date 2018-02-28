module xpc_m

  use xpc_h
  use log_m
  use const_numphys_h
  use const_kind_m
  use clcoef_h
  use clcoef_m

  implicit none
  private

! public subroutines
  public :: &
  xpc_initfile,  & !< open file
  xpc_readcon,  & !< read data from file
  xpc_solve,  & !< generic subroutine
  xpc_userdefined,  & !< user-defined function
  xpc_fn, &  !< general external function call
  xpc_dia, &  !< object diagnostics to log file
  xpc_initwrite, & !< open new file, making up name
  xpc_write, &  !< write out object
  xpc_writeg, &  !< write out object as gnuplot
  xpc_writev, &  !< write out object as vtk
  xpc_delete, & !< delete object
  xpc_close, & !< close file
  xpc_closewrite !< close write file

! private variables
  character(*), parameter :: m_name='xpc_m' !< module name
  integer(ki4)  :: status   !< error status
  integer(ki4), save  :: ninxp=5     !< control file unit number
  integer(ki4), save  :: noutxp=6      !< output file unit number
  character(len=80), save :: controlfile !< control file name
  character(len=80), save :: outputfile !< output file name
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: ij !< loop counter
  integer(ki4)  :: ilog      !< for namelist dump after error

  contains
!---------------------------------------------------------------------
!> open file
subroutine xpc_initfile(file,channel)

  !! arguments
  character(*), intent(in) :: file !< file name
  integer(ki4), intent(out),optional :: channel   !< input channel for object data structure
  !! local
  character(*), parameter :: s_name='xpc_initfile' !< subroutine name
  logical :: unitused !< flag to test unit is available

  if (trim(file)=='null') then
     call log_error(m_name,s_name,1,log_info,'null filename ignored')
     return
  end if

  !! get file unit
  do i=99,1,-1
     inquire(i,opened=unitused)
     if(.not.unitused)then
        ninxp=i
        if (present(channel)) channel=i
        exit
     end if
  end do

  !! open file
  controlfile=trim(file)
  call log_value("Control data file",trim(controlfile))
  open(unit=ninxp,file=controlfile,status='OLD',iostat=status)
  if(status/=0)then
     !! error opening file
     print '("Fatal error: Unable to open control file, ",a)',controlfile
     call log_error(m_name,s_name,2,error_fatal,'Cannot open control data file')
     stop
  end if

end subroutine xpc_initfile
!---------------------------------------------------------------------
!> read data from file
subroutine xpc_readcon(selfn,channel)

  !! arguments
  type(xpnumerics_t), intent(out) :: selfn !< type which data will be assigned to
  integer(ki4), intent(in),optional :: channel   !< input channel for object data structure

  !! local
  character(*), parameter :: s_name='xpc_readcon' !< subroutine name
  character(len=80) :: xpc_formula !< formula to be used
  integer(ki4), parameter :: MAX_NUMBER_OF_PARAMETERS=10 !< maximum number of parameters allowed
  real(kr8) :: power_split !< variable with meaningful name
  real(kr8) :: atomic_number !< self-explanatory local
  real(kr8) :: effective_charge !< self-explanatory local
  real(kr8) :: imposed_Bfield !< self-explanatory local
  real(kr8) :: electron_temperature !< self-explanatory local
  real(kr8) :: ion_temperature !< self-explanatory local
  real(kr8) :: number_density_factor !< self-explanatory local
  real(kr8) :: electron_number_density !< self-explanatory local
  real(kr8) :: Coulomb_logarithm_factor !< self-explanatory local
  real(kr8) :: Coulomb_logarithm !< self-explanatory local
  real(kr8) :: layer_depth !< self-explanatory local
  real(kr8) :: pressure_lengthscale !< self-explanatory local
  real(kr8) :: poloidal_angle !< self-explanatory local
  real(kr8) :: major_radius !< self-explanatory local
  real(kr8) :: minor_radius !< self-explanatory local
  real(kr8) :: effective_Bfield !< self-explanatory local
  character(80) :: effective_field_formula !< self-explanatory local

  real(kr8), dimension(MAX_NUMBER_OF_PARAMETERS) :: general_real_parameters  !< local variable
  integer(ki4), dimension(MAX_NUMBER_OF_PARAMETERS) :: general_integer_parameters  !< local variable
  integer(ki4) :: number_of_real_parameters  !< local variable
  integer(ki4) :: number_of_integer_parameters  !< local variable

  !! xpc parameters
  namelist /xpcparameters/ &
 & atomic_number , &
 & effective_charge , &
 & imposed_Bfield , &
 & electron_temperature , &
 & ion_temperature , &
 & number_density_factor , &
 & electron_number_density , &
 & Coulomb_logarithm_factor , &
 & Coulomb_logarithm , &
 & layer_depth , &
 & pressure_lengthscale , &
 & poloidal_angle , &
 & major_radius , &
 & minor_radius , &
 & effective_Bfield , &
 & effective_field_formula , &
 &power_split, xpc_formula, &
 &general_real_parameters, number_of_real_parameters, &
 &general_integer_parameters, number_of_integer_parameters

  !! set default xpc parameters
  power_split=0.5_kr8

  xpc_formula='unset'
  general_real_parameters=0
  general_integer_parameters=0
  number_of_real_parameters=0
  number_of_integer_parameters=0
  atomic_number = 1
  effective_charge = 1
  imposed_Bfield = 1
  electron_temperature = 10
  ion_temperature = 10
  number_density_factor = 1.e+18
  electron_number_density = 1.e+18
  Coulomb_logarithm_factor = 1
  Coulomb_logarithm = 14
  layer_depth = 0.01
  pressure_lengthscale = 0.01
  poloidal_angle = 0.0
  major_radius = 1.0
  minor_radius = 0.1
  effective_Bfield = 1.0
  effective_field_formula = 'null'

  if(present(channel).AND.channel/=0) then
     !! assume unit already open and reading infile
     ninxp=channel
  end if

  !!read xpc parameters
  read(ninxp,nml=xpcparameters,iostat=status)
  if(status/=0) then
     !!dump namelist contents to logfile to assist error location
     print '("Fatal error reading xpc parameters")'
     call log_getunit(ilog)
     write(ilog,nml=xpcparameters)
     call log_error(m_name,s_name,1,error_fatal,'Error reading xpc parameters')
  end if

  call lowor(xpc_formula,1,len_trim(xpc_formula))
  !! check for valid data

  formula_chosen: select case (xpc_formula)
  case('unset','exp')
     if(power_split<0.OR.power_split>1) &
 &   call log_error(m_name,s_name,11,error_fatal,'power_split must be >=0 and <=1')

  case('expdouble')
     if(power_split<0.OR.power_split>1) &
 &   call log_error(m_name,s_name,21,error_fatal,'power_split must be >=0 and <=1')

  case('userdefined')
     if(number_of_real_parameters<0) &
 &   call log_error(m_name,s_name,44,error_fatal,'number of real parameters must be >=0')
     if(number_of_real_parameters>MAX_NUMBER_OF_PARAMETERS) then
        call log_value("max number of real parameters",MAX_NUMBER_OF_PARAMETERS)
        call log_error(m_name,s_name,45,error_fatal,'too many parameters: increase MAX_NUMBER_OF_PARAMETERS')
     end if
     if(number_of_integer_parameters<0) &
 &   call log_error(m_name,s_name,46,error_fatal,'number of integer parameters must be >=0')
     if(number_of_integer_parameters>MAX_NUMBER_OF_PARAMETERS) then
        call log_value("max number of integer parameters",MAX_NUMBER_OF_PARAMETERS)
        call log_error(m_name,s_name,47,error_fatal,'too many parameters: increase MAX_NUMBER_OF_PARAMETERS')
     end if
     if(number_of_integer_parameters==0.AND.number_of_real_parameters==0) &
 &   call log_error(m_name,s_name,48,error_fatal,'no parameters set')

  end select formula_chosen

  !! store values
  selfn%formula=xpc_formula

  selfn%f=power_split

  !! allocate arrays and assign

  selfn%nrpams=number_of_real_parameters
  selfn%nipams=number_of_integer_parameters

  formula_allocate: select case (xpc_formula)

  case('userdefined')
     if (number_of_real_parameters>0) then
        allocate(selfn%rpar(number_of_real_parameters), stat=status)
        call log_alloc_check(m_name,s_name,65,status)
        selfn%rpar=general_real_parameters(:number_of_real_parameters)
     end if
     if (number_of_integer_parameters>0) then
        allocate(selfn%npar(number_of_integer_parameters), stat=status)
        call log_alloc_check(m_name,s_name,66,status)
        selfn%npar=general_integer_parameters(:number_of_integer_parameters)
     end if
  case default
  end select formula_allocate
  selfn%a = atomic_number
  selfn%z = effective_charge
  selfn%b = imposed_Bfield
  selfn%t_e = electron_temperature
  selfn%t_i = ion_temperature
  selfn%c_n = number_density_factor
  selfn%n = electron_number_density
  selfn%c_lambda = Coulomb_logarithm_factor
  selfn%lambda = Coulomb_logarithm
  selfn%depth = layer_depth
  selfn%lpscale = pressure_lengthscale
  selfn%polang = poloidal_angle
  selfn%rmajor = major_radius
  selfn%rminor = minor_radius
  selfn%b1 = effective_Bfield
  selfn%b_formula = effective_field_formula

end  subroutine xpc_readcon
!---------------------------------------------------------------------
!> generic subroutine
subroutine xpc_solve(self)

  !! arguments
  type(xpc_t), intent(inout) :: self !< module object
  !! local
  character(*), parameter :: s_name='xpc_solve' !< subroutine name

  call clcoef_solve(self%clcoef,self%n)

end subroutine xpc_solve
!---------------------------------------------------------------------
!> output to log file
subroutine xpc_dia(self)

  !! arguments
  type(xpc_t), intent(inout) :: self !< module object
  !! local
  character(*), parameter :: s_name='xpc_dia' !< subroutine name

  call log_value("power ",self%pow)

end subroutine xpc_dia
!---------------------------------------------------------------------
!> userdefined function
function xpc_userdefined(self,psi)

  !! arguments
  type(xpc_t), intent(in) :: self !< module object
  real(kr8) :: xpc_userdefined !< local variable
  real(kr8), intent(in) :: psi !< position in \f$ \psi \f$

  !! local variables
  character(*), parameter :: s_name='xpc_userdefined' !< subroutine name
  real(kr8) :: pow !< local variable
  real(kr8) :: zpos !< position
  integer(ki4) :: ilocal !< local integer variable

  zpos=psi
  pow=0._kr8
  !> user defines \Tt{pow} here
  !! .....
  !! return xpc
  xpc_userdefined=pow

end function xpc_userdefined
!---------------------------------------------------------------------
!> general external function call
function xpc_fn(self,psi)

  !! arguments
  type(xpc_t), intent(in) :: self !< module object
  real(kr8) :: xpc_fn !< local variable
  real(kr8), intent(in) :: psi !< position in \f$ \psi \f$

  !! local variables
  character(*), parameter :: s_name='xpc_fn' !< subroutine name
  real(kr8) :: pow !< local variable

  pow=0._kr8
  !! select xpc
  formula_chosen: select case (self%n%formula)
  case('userdefined')
     pow=xpc_userdefined(self,psi)
  end select formula_chosen

  !! return xpc
  xpc_fn=pow

end function xpc_fn
!---------------------------------------------------------------------
!> open new file, making up name
subroutine xpc_initwrite(fileroot,channel)

  !! arguments
  character(*), intent(in) :: fileroot !< file root
  integer(ki4), intent(out),optional :: channel   !< output channel for object data structure
  !! local
  character(*), parameter :: s_name='xpc_initwrite' !< subroutine name
  logical :: unitused !< flag to test unit is available
  character(len=80) :: outputfile !< output file name

  !! get file unit
  do i=99,1,-1
     inquire(i,opened=unitused)
     if(.not.unitused)then
        if (present(channel)) channel=i
        exit
     end if
  end do
  noutxp=i

  !! open file
  outputfile=trim(fileroot)//"_xpc.out"
  call log_value("Control data file",trim(outputfile))
  open(unit=noutxp,file=outputfile,status='NEW',iostat=status)
  if(status/=0)then
     open(unit=noutxp,file=outputfile,status='REPLACE',iostat=status)
  end if
  if(status/=0)then
     !! error opening file
     print '("Fatal error: Unable to open output file, ",a)',outputfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot open output data file')
     stop
  end if

end subroutine xpc_initwrite
!---------------------------------------------------------------------
!> write xpc data
subroutine xpc_write(self,channel)

  !! arguments
  type(xpc_t), intent(in) :: self   !< xpc data structure
  integer(ki4), intent(in), optional :: channel   !< output channel for xpc data structure

  !! local
  character(*), parameter :: s_name='xpc_write' !< subroutine name
  integer(ki4) :: iout   !< output channel for xpc data structure

  !! sort out unit
  if(present(channel)) then
     iout=channel
  else
     iout=noutxp
  end if

  call clcoef_write(self%clcoef,iout)
 
  call clcoef_numwrite(self%clcoef,iout)

end subroutine xpc_write
!---------------------------------------------------------------------
!> write object data as gnuplot
subroutine xpc_writeg(self,select,channel)

  !! arguments
  type(xpc_t), intent(in) :: self   !< object data structure
  character(*), intent(in) :: select  !< case
  integer(ki4), intent(in), optional :: channel   !< output channel for xpc data structure

  !! local
  character(*), parameter :: s_name='xpc_writeg' !< subroutine name
  integer(ki4) :: iout   !< output channel for xpc data structure

  call log_error(m_name,s_name,1,log_info,'gnuplot file produced')

  plot_type: select case(select)
  case('cartesian')

  case default

  end select plot_type

end subroutine xpc_writeg
!---------------------------------------------------------------------
!> write object data as vtk
subroutine xpc_writev(self,select,channel)

  !! arguments
  type(xpc_t), intent(in) :: self   !< object data structure
  character(*), intent(in) :: select  !< case
  integer(ki4), intent(in), optional :: channel   !< output channel for xpc data structure

  !! local
  character(*), parameter :: s_name='xpc_writev' !< subroutine name
  integer(ki4) :: iout   !< output channel for xpc data structure

  call log_error(m_name,s_name,1,log_info,'vtk file produced')

  plot_type: select case(select)
  case('cartesian')

  case default

  end select plot_type

end subroutine xpc_writev
!---------------------------------------------------------------------
!> close write file
subroutine xpc_closewrite

  !! local
  character(*), parameter :: s_name='xpc_closewrite' !< subroutine name

  !! close file
  close(unit=noutxp,iostat=status)
  if(status/=0)then
     !! error closing file
     print '("Fatal error: Unable to close output file, ",a)',outputfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot close output data file')
     stop
  end if

end subroutine xpc_closewrite
!---------------------------------------------------------------------
!> delete object
subroutine xpc_delete(self)

  !! arguments
  type(xpc_t), intent(inout) :: self !< module object
  !! local
  character(*), parameter :: s_name='xpc_delete' !< subroutine name

  formula_deallocate: select case (self%n%formula)
  case('userdefined')
     if (self%n%nrpams>0) deallocate(self%n%rpar)
     if (self%n%nipams>0) deallocate(self%n%npar)
  case default
  end select formula_deallocate

end subroutine xpc_delete
!---------------------------------------------------------------------
!> close file
subroutine xpc_close

  !! local
  character(*), parameter :: s_name='xpc_close' !< subroutine name

  !! close file
  close(unit=ninxp,iostat=status)
  if(status/=0)then
     !! error closing file
     print '("Fatal error: Unable to close control file, ",a)',controlfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot close control data file')
     stop
  end if

end subroutine xpc_close

end module xpc_m
