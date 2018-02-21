module clcoef_m

  use clcoef_h
  use xpc_h
  use log_m
  use const_numphys_h
  use const_kind_m

  implicit none
  private

! public subroutines
  public :: &
  clcoef_initfile,  & !< open file
  clcoef_readcon,  & !< read data from file
  clcoef_solve,  & !< control subroutine
  clcoef_calc, & !< calculate classical coefficients
  clcoef_write, & !< write out classical coefficients
  clcoef_numcalc, & !< calculate classical numerical coefficients
  clcoef_numwrite, & !< write out classical numerical coefficients
  clcoef_userdefined,  & !< user-defined function
  clcoef_fn, &  !< general external function call
  clcoef_dia, &  !< object diagnostics to log file
  clcoef_initwrite, & !< open new file, making up name
  clcoef_writeg, &  !< write out object as gnuplot
  clcoef_writev, &  !< write out object as vtk
  clcoef_delete, & !< delete object
  clcoef_close, & !< close file
  clcoef_closewrite !< close write file

! private variables
  character(*), parameter :: m_name='clcoef_m' !< module name
  integer(ki4)  :: status   !< error status
  integer(ki4), save  :: nincc=5     !< control file unit number
  integer(ki4), save  :: noutcc=6      !< output file unit number
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
subroutine clcoef_initfile(file,channel)

  !! arguments
  character(*), intent(in) :: file !< file name
  integer(ki4), intent(out),optional :: channel   !< input channel for object data structure
  !! local
  character(*), parameter :: s_name='clcoef_initfile' !< subroutine name
  logical :: unitused !< flag to test unit is available

  !! get file unit
  do i=99,1,-1
     inquire(i,opened=unitused)
     if(.not.unitused)then
        nincc=i
        if (present(channel)) channel=i
        exit
     end if
  end do

  !! open file
  controlfile=trim(file)
  call log_value("Control data file",trim(controlfile))
  open(unit=nincc,file=controlfile,status='OLD',iostat=status)
  if(status/=0)then
     !! error opening file
     print '("Fatal error: Unable to open control file, ",a)',controlfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot open control data file')
     stop
  end if

end subroutine clcoef_initfile
!---------------------------------------------------------------------
!> read data from file
subroutine clcoef_readcon(selfn,channel)

  !! arguments
  type(ccnumerics_t), intent(out) :: selfn !< type which data will be assigned to
  integer(ki4), intent(in),optional :: channel   !< input channel for object data structure

  !! local
  character(*), parameter :: s_name='clcoef_readcon' !< subroutine name
  character(len=80) :: clcoef_formula !< formula to be used
  integer(ki4), parameter :: MAX_NUMBER_OF_PARAMETERS=10 !< maximum number of parameters allowed
  real(kr8) :: power_split !< variable with meaningful name

  real(kr8), dimension(MAX_NUMBER_OF_PARAMETERS) :: general_real_parameters  !< local variable
  integer(ki4), dimension(MAX_NUMBER_OF_PARAMETERS) :: general_integer_parameters  !< local variable
  integer(ki4) :: number_of_real_parameters  !< local variable
  integer(ki4) :: number_of_integer_parameters  !< local variable

  !! clcoef parameters
  namelist /clcoefparameters/ &
 &power_split, clcoef_formula, &
 &general_real_parameters, number_of_real_parameters, &
 &general_integer_parameters, number_of_integer_parameters

  !! set default clcoef parameters
  power_split=0.5_kr8

  clcoef_formula='unset'
  general_real_parameters=0
  general_integer_parameters=0
  number_of_real_parameters=0
  number_of_integer_parameters=0

  if(present(channel).AND.channel/=0) then
     !! assume unit already open and reading infile
     nincc=channel
  end if

  !!read clcoef parameters
  read(nincc,nml=clcoefparameters,iostat=status)
  if(status/=0) then
     !!dump namelist contents to logfile to assist error location
     print '("Fatal error reading clcoef parameters")'
     call log_getunit(ilog)
     write(ilog,nml=clcoefparameters)
     call log_error(m_name,s_name,1,error_fatal,'Error reading clcoef parameters')
  end if

  call lowor(clcoef_formula,1,len_trim(clcoef_formula))
  !! check for valid data

  formula_chosen: select case (clcoef_formula)
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
  selfn%formula=clcoef_formula

  selfn%f=power_split

  !! allocate arrays and assign

  selfn%nrpams=number_of_real_parameters
  selfn%nipams=number_of_integer_parameters

  formula_allocate: select case (clcoef_formula)

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

end  subroutine clcoef_readcon
!---------------------------------------------------------------------
!> generic subroutine
subroutine clcoef_solve(self,xpnumerics)

  !! arguments
  type(clcoef_t), intent(inout) :: self !< module object
  type(xpnumerics_t), intent(in) :: xpnumerics !< control  parameters

  !! local
  character(*), parameter :: s_name='clcoef_solve' !< subroutine name

  call clcoef_numcalc(self,xpnumerics)

  call clcoef_calc(self,xpnumerics)

end subroutine clcoef_solve
!---------------------------------------------------------------------
!> calculate classical coefficients
subroutine clcoef_calc(self,xpnumerics)

  !! arguments
  type(clcoef_t), intent(inout) :: self   !< object data structure
  type(xpnumerics_t), intent(in) :: xpnumerics !< control  parameters

  !! local
  character(*), parameter :: s_name='clcoef_calc' !< subroutine name
  real(kr8) :: m_i !< ion mass
  real(kr8) :: z_te !< convert temperature to energy
  real(kr8) :: a !< local variable
  real(kr8) :: z !< local variable
  real(kr8) :: b !< local variable
  real(kr8) :: n !< local variable
  real(kr8) :: t_e !< local variable
  real(kr8) :: t_i !< local variable
  real(kr8) :: c_n !< local variable
  real(kr8) :: lambda !< local variable

  a = xpnumerics%a
  z = xpnumerics%z
  b = xpnumerics%b
  n = xpnumerics%n
  t_e = xpnumerics%t_e
  t_i = xpnumerics%t_i
  c_n = xpnumerics%c_n
  lambda = xpnumerics%lambda

  !> temperatures in eV
  z_te = const_charge

  m_i=a*const_massp
  !5: omega_ce=e/const_masse*B;
  self%omega_ce=(b*const_charge)/const_masse
  !8: omega_ci=Z*const_charge*B/(A*m_p);
  self%omega_ci=(z*b*const_charge)/m_i
  !13: tau_e=6*sqrt(2*const_pid^3)*const_epsilon0^2*sqrt(const_masse)/e^4*(k_B*T_e)^(3/2)/N/Lambda;
  self%tau_e=(6*sqrt(const_pid)*sqrt(const_masse)*sqrt(z_te*t_e)*sqrt(2._kr8)*const_epsilon0**2 *&
 &z_te*const_pid*t_e)/(z**2*const_charge**4*lambda*n)

  !14: tau_i=12*sqrt(const_pid^3)*const_epsilon0^2*sqrt(const_massp)/e^4*(k_B*T_i)^(3/2)*sqrt(A)/N/Lambda;
  self%tau_i=(12*sqrt(const_pid)*sqrt(const_massp)*sqrt(a)*sqrt(z_te*t_i)*const_epsilon0**2 * &
 &z_te*const_pid*t_i)/(z**4*const_charge**4*lambda*n)

  !2: kappae_Bpara=3.2*N*k_B*T_e/const_masse*tau_e;;
  self%kappab_epara=(3.2*z_te*n*t_e*self%tau_e)/const_masse
  !4: kappae_Bperp=4.7*N*k_B*Te/const_masse*tau_e/(omega_ce*tau_e)^2;
  self%kappab_eperp=(4.7*z_te*n*t_e)/(const_masse*self%omega_ce**2*self%tau_e)
  !6: kappai_Bpara=3.9*N*k_B*Ti/m_i*tau_i;
  self%kappab_ipara=(3.9*z_te*n*self%tau_i*t_i)/m_i
  !7: kappai_Bperp=2*N*k_B*Ti/m_i*tau_i/(omega_ci*tau_i)^2;
  self%kappab_iperp=(2*z_te*n*t_i)/(m_i*self%omega_ci**2*self%tau_i)
  !9: R_Be=1/(omega_ce*tau_e);
  self%x_e=self%omega_ce*self%tau_e
  !10: R_Bi=1/(omega_ci*tau_i);
  self%x_i=self%omega_ci*self%tau_i
  !11: nu_epara=0.73*k_B*T_e/const_masse*tau_e;
  self%nu_epara=(0.73*z_te*t_e*self%tau_e)/const_masse
  !12: nu_ipara=0.96*k_B*T_i/m_i*tau_i;
  self%nu_ipara=(0.96*z_te*self%tau_i)/m_i
  !15: kappa_epara=13*sqrt(2*const_pid^3)/sqrt(const_masse)*const_epsilon0^2/e^4*(k_B*T_e)^(5/2)/N/Lambda;
  self%kappa_epara=(13*sqrt(const_pid)*sqrt(z_te*t_e)*sqrt(2._kr8)*const_epsilon0**2*  &
 &z_te**2*const_pid*t_e**2)/(z**2*sqrt(const_masse)*const_charge**4*lambda*n)

  !16: kappa_ipara=16*sqrt(const_pid^3)/sqrt(const_massp)*const_epsilon0^2/e^4*sqrt(A)*(k_B*T_i)^(5/2)/N/Lambda;
  self%kappa_ipara=(16*sqrt(const_pid)*sqrt(a)*sqrt(z_te*t_i)*const_epsilon0**2*  &
 &z_te**2*const_pid*t_i**2)/(z**4*sqrt(const_massp)*const_charge**4*lambda*n)

  !17: eta1=0.51*sqrt(const_masse)*e^2/(6*sqrt(2*const_pid^3)*const_mu0*const_epsilon0^2)*Z*Lambda;
  !18: eta=eta1/sqrt((k_B*T_e)**3);
  self%eta=(0.51*sqrt(const_masse)*const_charge**2*lambda*z)/(6*sqrt(const_pid)*sqrt(2._kr8)*  &
 &const_epsilon0**2*const_mu0*const_pid)/(sqrt(z_te*t_e)*abs(z_te)*abs(t_e))

  self%nu_eperp=self%nu_epara*(0.51/0.73)/self%x_e**2
  self%nu_iperp=self%nu_ipara*(0.3/0.96)/self%x_i**2

  self%kappa_iperp=(sqrt(const_massp)*const_charge**2*z**2*lambda*n*sqrt(a))/ &
 &(9*sqrt(const_pid)*const_pi*sqrt(z_te*t_i)*const_epsilon0**2*b**2)
  self%kappa_eperp=self%kappa_epara*(4.7/3.2)/self%x_e**2

end subroutine clcoef_calc
!---------------------------------------------------------------------
!> write out classical coefficients
subroutine clcoef_write(self,channel)

  !! arguments
  type(clcoef_t), intent(in) :: self   !< object data structure
  integer(ki4), intent(in) :: channel !< name of output file

  !! local
  character(*), parameter :: s_name='clcoef_write' !< subroutine name
  logical :: unitused !< flag to test unit is available
  character(30), parameter :: zcfmt='(A,1P,G14.5)' !< format sta
  
  write(channel,zcfmt) 'omega_ce = ', self%omega_ce
  write(channel,zcfmt) 'omega_ci = ', self%omega_ci
  write(channel,zcfmt) 'tau_e = ', self%tau_e
  write(channel,zcfmt) 'tau_i = ', self%tau_i
  write(channel,zcfmt) 'kappab_epara = ', self%kappab_epara
  write(channel,zcfmt) 'kappab_eperp = ', self%kappab_eperp
  write(channel,zcfmt) 'kappab_ipara = ', self%kappab_ipara
  write(channel,zcfmt) 'kappab_iperp = ', self%kappab_iperp
  write(channel,zcfmt) 'x_e = ', self%x_e
  write(channel,zcfmt) 'x_i = ', self%x_i
  write(channel,zcfmt) 'nu_epara = ', self%nu_epara
  write(channel,zcfmt) 'nu_eperp = ', self%nu_eperp
  write(channel,zcfmt) 'nu_ipara = ', self%nu_ipara
  write(channel,zcfmt) 'nu_iperp = ', self%nu_iperp
  write(channel,zcfmt) 'kappa_ipara = ', self%kappa_ipara
  write(channel,zcfmt) 'kappa_iperp = ', self%kappa_iperp
  write(channel,zcfmt) 'kappa_epara = ', self%kappa_epara
  write(channel,zcfmt) 'kappa_eperp = ', self%kappa_eperp
  write(channel,zcfmt) 'eta = ', self%eta

end subroutine clcoef_write
!---------------------------------------------------------------------
!> calculate classical numerical coefficients
subroutine clcoef_numcalc(self,xpnumerics)

  !! arguments
  type(clcoef_t), intent(inout) :: self   !< object data structure
  type(xpnumerics_t), intent(in) :: xpnumerics !< control  parameters

  !! local
  character(*), parameter :: s_name='clcoef_numcalc' !< subroutine name
  real(kr8) :: c_n !< local variable
  real(kr8) :: lambda !< local variable
  real(kr8) :: c_lambda !< local variable

  c_n = xpnumerics%c_n
  c_lambda = xpnumerics%c_lambda

  self%c_omega_ce=const_charge/const_masse

  self%c_omega_ci=const_charge/const_massp

  self%c_tau_e=(6*sqrt(const_pid)*sqrt(const_masse)*sqrt(const_charge)*sqrt(2._kr8)*const_epsilon0**2 *&
 &const_charge*const_pid)/(const_charge**4*c_lambda*c_n)
  self%c_tau_i=(12*sqrt(const_pid)*sqrt(const_massp)*sqrt(const_charge)*const_epsilon0**2 * &
 &const_charge*const_pid)/(const_charge**4*c_lambda*c_n)
  self%c_kappab_epara=(3.2*const_charge*c_n*self%c_tau_e)/const_masse

  self%c_kappab_eperp=(4.7*const_charge*c_n)/(const_masse*self%c_omega_ce**2*self%c_tau_e)

  self%c_kappab_ipara=(3.9*const_charge*c_n*self%c_tau_i)/const_massp

  self%c_kappab_iperp=(2*const_charge*c_n)/(const_massp*self%c_omega_ci**2*self%c_tau_i)

  self%c_x_e=self%c_omega_ce*self%c_tau_e

  self%c_x_i=self%c_omega_ci*self%c_tau_i

  self%c_nu_epara=(0.73*const_charge*self%c_tau_e)/const_masse

  self%c_nu_ipara=(0.96*const_charge*self%c_tau_i)/const_massp

  self%c_eta=(0.51*sqrt(const_masse)*const_charge**2*c_lambda)/(6*sqrt(const_pid)*sqrt(2._kr8)*  &
 &const_epsilon0**2*const_mu0*const_pid)/(sqrt(const_charge)*abs(const_charge))
  self%c_kappa_ipara=(16*sqrt(const_pid)*sqrt(const_charge)*const_epsilon0**2*  &
 &const_charge**2*const_pid)/(sqrt(const_massp)*const_charge**4*c_lambda*c_n)
  self%c_kappa_epara=(13*sqrt(const_pid)*sqrt(const_charge)*sqrt(2._kr8)*const_epsilon0**2*  &
 &const_charge**2*const_pid)/(sqrt(const_masse)*const_charge**4*c_lambda*c_n)

  self%c_nu_eperp=self%c_nu_epara*(0.51/0.73)/self%c_x_e**2
  self%c_nu_iperp=self%c_nu_ipara*(0.3/0.96)/self%c_x_i**2

  self%c_kappa_iperp=(sqrt(const_massp*const_charge)*const_charge*c_lambda*c_n)/ &
 &(9*sqrt(const_pid)*const_pi*const_epsilon0**2)
  self%c_kappa_eperp=self%c_kappa_epara*(4.7/3.2)/self%c_x_e**2

end subroutine clcoef_numcalc
!---------------------------------------------------------------------
!> write out classical numerical coefficients
subroutine clcoef_numwrite(self,channel)

  !! arguments
  type(clcoef_t), intent(in) :: self   !< object data structure
  integer(ki4), intent(in) :: channel !< name of output file

  !! local
  character(*), parameter :: s_name='clcoef_numwrite' !< subroutine name
  character(30), parameter :: zcfmt='(A,1P,G14.5)' !< format sta

  write(channel,zcfmt) 'c_omega_ce = ', self%c_omega_ce
  write(channel,zcfmt) 'c_omega_ci = ', self%c_omega_ci
  write(channel,zcfmt) 'c_tau_e = ', self%c_tau_e
  write(channel,zcfmt) 'c_tau_i = ', self%c_tau_i
  write(channel,zcfmt) 'c_kappab_epara = ', self%c_kappab_epara
  write(channel,zcfmt) 'c_kappab_eperp = ', self%c_kappab_eperp
  write(channel,zcfmt) 'c_kappab_ipara = ', self%c_kappab_ipara
  write(channel,zcfmt) 'c_kappab_iperp = ', self%c_kappab_iperp
  write(channel,zcfmt) 'c_x_e = ', self%c_x_e
  write(channel,zcfmt) 'c_x_i = ', self%c_x_i
  write(channel,zcfmt) 'c_nu_epara = ', self%c_nu_epara
  write(channel,zcfmt) 'c_nu_eperp = ', self%c_nu_eperp
  write(channel,zcfmt) 'c_nu_ipara = ', self%c_nu_ipara
  write(channel,zcfmt) 'c_nu_iperp = ', self%c_nu_iperp
  write(channel,zcfmt) 'c_kappa_ipara = ', self%c_kappa_ipara
  write(channel,zcfmt) 'c_kappa_iperp = ', self%c_kappa_iperp
  write(channel,zcfmt) 'c_kappa_epara = ', self%c_kappa_epara
  write(channel,zcfmt) 'c_kappa_eperp = ', self%c_kappa_eperp
  write(channel,zcfmt) 'c_eta = ', self%c_eta

end subroutine clcoef_numwrite
!---------------------------------------------------------------------
!> output to log file
subroutine clcoef_dia(self)

  !! arguments
  type(clcoef_t), intent(inout) :: self !< module object
  !! local
  character(*), parameter :: s_name='clcoef_dia' !< subroutine name

  call log_value("power ",self%pow)

end subroutine clcoef_dia
!---------------------------------------------------------------------
!> userdefined function
function clcoef_userdefined(self,psi)

  !! arguments
  type(clcoef_t), intent(in) :: self !< module object
  real(kr8) :: clcoef_userdefined !< local variable
  real(kr8), intent(in) :: psi !< position in \f$ \psi \f$

  !! local variables
  character(*), parameter :: s_name='clcoef_userdefined' !< subroutine name
  real(kr8) :: pow !< local variable
  real(kr8) :: zpos !< position
  integer(ki4) :: ilocal !< local integer variable

  zpos=psi
  pow=0._kr8
  !> user defines \Tt{pow} here
  !! .....
  !! return clcoef
  clcoef_userdefined=pow

end function clcoef_userdefined
!---------------------------------------------------------------------
!> general external function call
function clcoef_fn(self,psi)

  !! arguments
  type(clcoef_t), intent(in) :: self !< module object
  real(kr8) :: clcoef_fn !< local variable
  real(kr8), intent(in) :: psi !< position in \f$ \psi \f$

  !! local variables
  character(*), parameter :: s_name='clcoef_fn' !< subroutine name
  real(kr8) :: pow !< local variable

  !! select clcoef
  formula_chosen: select case (self%n%formula)
  case('userdefined')
     pow=clcoef_userdefined(self,psi)
  end select formula_chosen

  !! return clcoef
  clcoef_fn=pow

end function clcoef_fn
!---------------------------------------------------------------------
!> open new file, making up name
subroutine clcoef_initwrite(fileroot,channel)

  !! arguments
  character(*), intent(in) :: fileroot !< file root
  integer(ki4), intent(out),optional :: channel   !< output channel for object data structure
  !! local
  character(*), parameter :: s_name='clcoef_initwrite' !< subroutine name
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
  noutcc=i

  !! open file
  outputfile=trim(fileroot)//"_clcoef.out"
  call log_value("Control data file",trim(outputfile))
  open(unit=noutcc,file=outputfile,status='NEW',iostat=status)
  if(status/=0)then
     open(unit=noutcc,file=outputfile,status='REPLACE',iostat=status)
  end if
  if(status/=0)then
     !! error opening file
     print '("Fatal error: Unable to open output file, ",a)',outputfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot open output data file')
     stop
  end if

end subroutine clcoef_initwrite
!---------------------------------------------------------------------
!> write object data as gnuplot
subroutine clcoef_writeg(self,select,channel)

  !! arguments
  type(clcoef_t), intent(in) :: self   !< object data structure
  character(*), intent(in) :: select  !< case
  integer(ki4), intent(in), optional :: channel   !< output channel for clcoef data structure

  !! local
  character(*), parameter :: s_name='clcoef_writeg' !< subroutine name
  integer(ki4) :: iout   !< output channel for clcoef data structure

  call log_error(m_name,s_name,1,log_info,'gnuplot file produced')

  plot_type: select case(select)
  case('cartesian')

  case default

  end select plot_type

end subroutine clcoef_writeg
!---------------------------------------------------------------------
!> write object data as vtk
subroutine clcoef_writev(self,select,channel)

  !! arguments
  type(clcoef_t), intent(in) :: self   !< object data structure
  character(*), intent(in) :: select  !< case
  integer(ki4), intent(in), optional :: channel   !< output channel for clcoef data structure

  !! local
  character(*), parameter :: s_name='clcoef_writev' !< subroutine name
  integer(ki4) :: iout   !< output channel for clcoef data structure

  call log_error(m_name,s_name,1,log_info,'vtk file produced')

  plot_type: select case(select)
  case('cartesian')

  case default

  end select plot_type

end subroutine clcoef_writev
!---------------------------------------------------------------------
!> close write file
subroutine clcoef_closewrite

  !! local
  character(*), parameter :: s_name='clcoef_closewrite' !< subroutine name

  !! close file
  close(unit=noutcc,iostat=status)
  if(status/=0)then
     !! error closing file
     print '("Fatal error: Unable to close output file, ",a)',outputfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot close output data file')
     stop
  end if

end subroutine clcoef_closewrite
!---------------------------------------------------------------------
!> delete object
subroutine clcoef_delete(self)

  !! arguments
  type(clcoef_t), intent(inout) :: self !< module object
  !! local
  character(*), parameter :: s_name='clcoef_delete' !< subroutine name

  formula_deallocate: select case (self%n%formula)
  case('userdefined')
     if (self%n%nrpams>0) deallocate(self%n%rpar)
     if (self%n%nipams>0) deallocate(self%n%npar)
  case default
  end select formula_deallocate

end subroutine clcoef_delete
!---------------------------------------------------------------------
!> close file
subroutine clcoef_close

  !! local
  character(*), parameter :: s_name='clcoef_close' !< subroutine name

  !! close file
  close(unit=nincc,iostat=status)
  if(status/=0)then
     !! error closing file
     print '("Fatal error: Unable to close control file, ",a)',controlfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot close control data file')
     stop
  end if

end subroutine clcoef_close

end module clcoef_m
