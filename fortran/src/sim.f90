module sim_m
use babu_m
use jsonx_m, only: jsonx_load, jsonx_get
use json_m, only: json_value, json_get, json_failed, json_clear_exceptions
use linalg_mod, only: solve_linear_system
use connection_m
use udp_m, only: udp_initialize, udp_finalize
implicit none

real :: mass
real :: I(3,3)
real :: Iinv(3,3)
real :: FM(6)
real :: h(3)
real :: y_init(13)
real :: controls(4)
real :: rho0
real, allocatable :: aero_ref_location(:)

type(json_value), pointer :: j_main

type(connection) :: graphics 
type(connection) :: controls_in

type stall_settings_t
  real :: alpha_0, alpha_s, lambda_b, minval
end type stall_settings_t

! aerodynamic coefficients
real :: sref , long_ref, lat_ref
real :: CL0 , CL_alpha , CL_qbar , CL_alphacap , CL_deltae
real :: CS_beta , CS_alphapbar , CS_pbar , CS_rbar , CS_deltaa , CS_deltar
real :: CD_L0 , CD_L , CD_L2 , CD_S2 , CD_qbar , CD_alpha_qbar , CD_deltae , &
        CD_alpha_deltae, CD_deltae2
real :: Cl_beta, Cl_pbar, Cl_rbar , Cl_alpha_rbar, Cl_deltaa , Cl_deltar
real :: Cm0 , Cm_alpha, Cm_qbar , Cm_alphacap, Cm_deltae
real :: Cn_beta , Cn_pbar , Cn_alphapbar , Cn_rbar , Cn_deltaa , Cn_alpha_deltaa, &
        Cn_deltar
real :: T0, Ta

! stall model constants
logical :: include_stall
type(stall_settings_t) :: CLstall, CDstall, Cmstall

! actuator dynamics 
logical :: use_actuator_dynamics
real :: actuator_states(4)        ! δ(t) actual actuator positions [aileron, elevator, rudder, throttle]
real :: sigma_actuator(4)         ! σ  damping rates (1/time_constant) for each actuator
real :: actuator_rate_limits(4)   ! Maximum rates [δ_dot_max]
real :: actuator_pos_limits(4,2)  ! Position limits [min, max] for each actuator

logical :: rk4_verbose
integer :: rk4_file = 20
real :: dt, tf

! Trim variables
character(10) :: trim_type
real :: trim_V, trim_alt, trim_phi, trim_theta, trim_gamma
logical :: trim_verbose, print_state_aero, use_climb_angle

contains

function json_key_exists(j_ptr, keypath) result(exists)
  use json_m, only: json_value, json_get, json_failed, json_clear_exceptions
  implicit none
  type(json_value), pointer, intent(in) :: j_ptr
  character(*), intent(in) :: keypath
  logical :: exists
  type(json_value), pointer :: temp_ptr
  logical :: found

  call json_get(j_ptr, keypath, temp_ptr, found)
  exists = found .and. .not. json_failed()
  call json_clear_exceptions()

end function json_key_exists


function compute_theta_from_gamma(gamma_rad, u, v, w, phi_rad) result(theta_rad)
  implicit none
  real, intent(in) :: gamma_rad  ! Climb angle in radians
  real, intent(in) :: u, v, w     ! Body-axis velocity components
  real, intent(in) :: phi_rad     ! Bank angle in radians
  real :: theta_rad               ! Elevation angle in radians 
  real :: Vtotal, Sphi, Cphi, Stheta_plus, Stheta_minus
  real :: numerator, denom, discriminant, theta_plus, theta_minus
  real :: error_plus, error_minus, Sgamma


  Vtotal = sqrt(u**2 + v**2 + w**2)


  Sphi = sin(phi_rad)
  Cphi = cos(phi_rad)
  Sgamma = sin(gamma_rad)

  ! quadratic equation (7.2.7) for S_theta using equation (7.2.8)
  ! The equation is:
  ! S_theta = [u*V*Sgamma ± (v*Sphi + w*Cphi)*sqrt(u^2 + (v*Sphi + w*Cphi)^2 - V^2*Sgamma^2)] / [u^2 + (v*Sphi + w*Cphi)^2]

  numerator = u * Vtotal * Sgamma
  denom = u**2 + (v*Sphi + w*Cphi)**2

  discriminant = u**2 + (v*Sphi + w*Cphi)**2 - Vtotal**2 * Sgamma**2

  if (discriminant < 0.0) then
    write(*,'(A)') 'Warning: Invalid discriminant in theta calculation'
    theta_rad = 0.0
    return
  end if

  ! Two solutions from quadratic formula
  Stheta_plus  = (numerator + (v*Sphi + w*Cphi)*sqrt(discriminant)) / denom
  Stheta_minus = (numerator - (v*Sphi + w*Cphi)*sqrt(discriminant)) / denom


  if (abs(Stheta_plus) <= 1.0) then
    theta_plus = asin(Stheta_plus)
  else
    theta_plus = 0.0
  end if

  if (abs(Stheta_minus) <= 1.0) then
    theta_minus = asin(Stheta_minus)
  else
    theta_minus = 0.0
  end if

  error_plus  = abs(Vtotal*Sgamma - (u*sin(theta_plus)  - (v*Sphi + w*Cphi)*cos(theta_plus)))
  error_minus = abs(Vtotal*Sgamma - (u*sin(theta_minus) - (v*Sphi + w*Cphi)*cos(theta_minus)))

  if (error_plus < error_minus) then
    theta_rad = theta_plus
  else
    theta_rad = theta_minus
  end if
end function compute_theta_from_gamma

function rk4(t0, y0, dt) result(ans)
  implicit none
  real, intent(in) :: t0, dt
  real, intent(in) :: y0(13)
  real :: ans(13)
  real :: k1(13), k2(13), k3(13), k4(13)
  k1 = diff_eq(t0, y0,1)
  k2 = diff_eq(t0 + 0.5*dt, y0 + 0.5*dt*k1,2)
  k3 = diff_eq(t0 + 0.5*dt, y0 + 0.5*dt*k2,3)
  k4 = diff_eq(t0 + dt, y0 + dt*k3,4)
  ans = y0 + (dt/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4)
end function rk4

function diff_eq(t, y,call_num) result(ans)
  implicit none
  real, intent(in) :: t
  real, intent(in) :: y(13)
  integer, intent(in) :: call_num
  real :: ans(13)
  real :: dydt(13), Mt(3)
  real :: u, v, w, p, q, r, e0, ex, ey, ez
  real :: Ixx, Iyy, Izz, Ixy, Ixz, Iyz
  real :: g, hx, hy,hz

  call pseudo_aero(y)

  if(rk4_verbose .or. (call_num == 1 .and. abs(y(1)-330.0) < 0.1 .and. abs(y(4)) < 0.1)) then
    if (.not. rk4_verbose) then
      write(*,'(A,6(ES24.15,2X))') ' FM (body) = ', FM(:)
      write(*,'(A,ES24.15)') ' mass = ', mass
    else
      write(rk4_file,*)
      write(rk4_file,*) 'diff eq function called'
      write(rk4_file,*) 'RK call number = ', call_num
      write(rk4_file,'(A,ES20.12)') ' TIME [s]=',t
      write(rk4_file,'(A,13ES20.12)') ' State in =',y(:)
      write(rk4_file,'(A,6ES20.12)') ' FM (body)=',FM(:)
    end if
  end if

  hx = h(1)
  hy = h(2)
  hz = h(3)

  g = gravity_f(-y(9))

  Ixx = I(1,1)
  Iyy = I(2,2)
  Izz = I(3,3)
  Ixy = -I(1,2)
  Ixz = -I(1,3)
  Iyz = -I(2,3)

  u = y(1); v = y(2); w = y(3)
  p = y(4); q = y(5); r = y(6)
  e0 = y(10); ex = y(11); ey = y(12); ez = y(13)
 
  ans(1) = 2.0*g*(ex*ez - ey*e0) + FM(1)/mass + r*v - q*w
  ans(2) = 2.0*g*(ey*ez + ex*e0) + FM(2)/mass + p*w - r*u
  ans(3) = g*( e0*e0 + ez*ez - ex*ex - ey*ey ) + FM(3)/mass + q*u - p*v
  
  Mt(1) = -hz*q + hy*r + FM(4) + (Iyy - Izz)*q*r + Iyz*(q*q-r*r)+Ixz*p*q-Ixy*p*r
  Mt(2) = hz*p - hx*r + FM(5) + (Izz - Ixx)*p*r + Ixz*(r*r-p*p)+Ixy*q*r-Iyz*p*q
  Mt(3) = -hy*p + hx*q + FM(6) + (Ixx - Iyy)*p*q + Ixy*(p*p-q*q)+Iyz*p*r-Ixz*q*r
  
  ans(4) = Iinv(1,1)*Mt(1) + Iinv(1,2)*Mt(2) + Iinv(1,3)*Mt(3)
  ans(5) = Iinv(2,1)*Mt(1) + Iinv(2,2)*Mt(2) + Iinv(2,3)*Mt(3)
  ans(6) = Iinv(3,1)*Mt(1) + Iinv(3,2)*Mt(2) + Iinv(3,3)*Mt(3)
  ans(7:9) = quat_dependent_to_base(y(1:3), y(10:13))
  ans(10) = 0.5*( -p*ex - q*ey - r*ez )
  ans(11) = 0.5*( p*e0 + r*ey - q*ez )
  ans(12) = 0.5*( q*e0 - r*ex + p*ez )
  ans(13) = 0.5*( r*e0 + q*ex - p*ey )
  if(rk4_verbose) then
    write(rk4_file,'(A,13ES20.12)') ' differential equation results =',ans(:)
  end if
  return
end function diff_eq

subroutine update_actuators(commanded, dt_step)
  implicit none
  real, intent(in) :: commanded(4)  ! u(t) Commanded control positions
  real, intent(in) :: dt_step       ! Time step for integration
  real :: error(4), delta_dot(4), delta_new(4)
  integer :: i
  
  if (.not. use_actuator_dynamics) then
    ! No dynamics so instant response (δ = u)
    actuator_states = commanded
    return
  end if
  
  ! For each actuator: aileron(1), elevator(2), rudder(3), throttle(4)
  do i = 1, 4
    ! Computing error: u(t) - δ(t)
    error(i) = commanded(i) - actuator_states(i)
    
    ! First-order dynamics using Eq. 10.2.1 from Chapter 10
    delta_dot(i) = sigma_actuator(i) * error(i)
    
    ! Apply rate limiting
    if (abs(delta_dot(i)) > actuator_rate_limits(i)) then
      delta_dot(i) = sign(actuator_rate_limits(i), delta_dot(i))
    end if
    
    ! Euler integration
    delta_new(i) = actuator_states(i) + delta_dot(i) * dt_step
    
    ! Apply position saturation limits: δ_min ≤ δ ≤ δ_max
    if (delta_new(i) < actuator_pos_limits(i,1)) then
      delta_new(i) = actuator_pos_limits(i,1)
    else if (delta_new(i) > actuator_pos_limits(i,2)) then
      delta_new(i) = actuator_pos_limits(i,2)
    end if
    
    ! To prevent overshoot 
    if (abs(commanded(i) - actuator_states(i)) < abs(delta_dot(i) * dt_step)) then
      actuator_states(i) = commanded(i)
    else
      actuator_states(i) = delta_new(i)
    end if
  end do
  
end subroutine update_actuators

subroutine pseudo_aero(y)
  implicit none
  real, intent(in) :: y(13)
  real :: da, de, dr, tau
  real :: V, Vlat, alpha, beta, pbar, qbar, rbar, ahat
  real :: CL1, CLift, CSide, CDrag, Croll, Cpitch, Cyaw
  real :: CLnewt, CDnewt,Cmnewt, pos, neg, sigma, signa
  real :: sa, ca, sb, cb
  real :: Z,h_ft,T,P,rho,a,mu
  
  ahat = 0.0
  da = controls(1)
  de =controls(2)
  dr = controls(3)
  tau = controls(4)

  call standard_atmospheric_English(-y(9),Z,T,P,rho,a)
  
  V = sqrt(y(1)**2 + y(2)**2 + y(3)**2)
  Vlat = sqrt(max(y(1)**2 + y(3)**2, 1.0e-16))
  alpha = atan2(y(3), y(1))
  beta = atan2(y(2), Vlat)

    pbar = 0.5*y(4)*lat_ref / V
    qbar = 0.5*y(5)*long_ref / V
    rbar = 0.5*y(6)*lat_ref / V

  CL1 = CL0 + CL_alpha*alpha
  CLift = CL1 + CL_qbar*qbar + CL_alphacap*ahat + CL_deltae*de
  CSide = CS_beta*beta + (CS_pbar + CS_alphapbar*alpha)*pbar + CS_rbar*rbar + &
          CS_deltaa*da + CS_deltar*dr
  CDrag = CD_L0 + CD_L*CL1 + CD_L2*CL1**2 + CD_S2*CSide**2 + (CD_qbar + &
          CD_alpha_qbar*alpha)*qbar + (CD_deltae + CD_alpha_deltae*alpha)*de + &
          CD_deltae2*de**2
  Croll = Cl_beta*beta + Cl_pbar*pbar + (Cl_rbar + Cl_alpha_rbar*alpha)*rbar + &
          Cl_deltaa*da + Cl_deltar*dr
  Cpitch = Cm0 + Cm_alpha*alpha + Cm_qbar*qbar + Cm_alphacap*ahat + Cm_deltae*de
  Cyaw = Cn_beta*beta + (Cn_pbar + Cn_alphapbar*alpha)*pbar + Cn_rbar*rbar + &
          (Cn_deltaa + Cn_alpha_deltaa*alpha)*da + Cn_deltar*dr

  sa = sin(alpha); ca = cos(alpha); sb = sin(beta); cb = cos(beta)
  signa = sign(1.0,alpha)

  if (include_stall) then
    CLnewt = 2.0*signa*sa*sa*ca
    pos = exp(CLstall%lambda_b*(alpha - CLstall%alpha_0 + CLstall%alpha_s))
    neg = exp(-CLstall%lambda_b*(alpha - CLstall%alpha_0 - CLstall%alpha_s))
    sigma = (1.0 + neg + pos)/((1.0 + neg)*(1.0 + pos))
    CLift = (1.0 - sigma)*CLift + sigma*CLnewt

    CDnewt = 2.0*sin(abs(alpha))**3
    pos = exp(CDstall%lambda_b*(alpha - CDstall%alpha_0 + CDstall%alpha_s))
    neg = exp(-CDstall%lambda_b*(alpha - CDstall%alpha_0 - CDstall%alpha_s))
    sigma = (1.0 + neg + pos)/((1.0 + neg)*(1.0 + pos))
    CDrag = (1.0 - sigma)*CDrag + sigma*CDnewt

    Cmnewt = 2.0*Cmstall%minval*signa*sa*sa*ca
    pos = exp(Cmstall%lambda_b*(alpha - Cmstall%alpha_0 + Cmstall%alpha_s))
    neg = exp(-Cmstall%lambda_b*(alpha - Cmstall%alpha_0 - Cmstall%alpha_s))
    sigma = (1.0 + neg + pos)/((1.0 + neg)*(1.0 + pos))
    Cpitch = (1.0 - sigma)*Cpitch + sigma*Cmnewt
  end if

  FM(1) = -(CDrag*ca*cb + CSide*ca*sb - CLift*sa)
  FM(2) = (CSide*cb - CDrag*sb)
  FM(3) = -(CDrag*sa*cb + CSide*sa*sb + CLift*ca)
 
  FM(4) = lat_ref * Croll
  FM(5) = long_ref * Cpitch
  FM(6) = lat_ref * Cyaw

  ! Debug output for initial state
  ! if (abs(y(1)-330.0) < 0.1 .and. abs(y(4)) < 0.1 .and. abs(controls(1)) < 0.1) then
  !   write(*,'(A)') ' --- AERO DEBUG ---'
  !   write(*,'(A,ES24.15)') ' Altitude = ', -y(9)
  !   write(*,'(A,ES24.15)') ' rho = ', rho
  !   write(*,'(A,ES24.15)') ' V = ', V
  !   write(*,'(A,ES24.15)') ' alpha = ', alpha
  !   write(*,'(A,ES24.15)') ' CL1 = ', CL1
  !   write(*,'(A,ES24.15)') ' CLift = ', CLift
  !   write(*,'(A,ES24.15)') ' CDrag = ', CDrag
  !   write(*,'(A,ES24.15)') ' CSide = ', CSide
  !   write(*,'(A,ES24.15)') ' Cpitch = ', Cpitch
  !   write(*,'(A,ES24.15)') ' FM(5) before scaling = ', long_ref * Cpitch
  !   write(*,'(A)') ' ------------------'
  ! end if

  FM = 0.5*rho*V**2*sref*FM
  
  call standard_atmospheric_English(0.0,z,T,P,rho0,a)
  FM(1) = FM(1) + tau*T0*(rho/rho0)**Ta
  
  FM(4:6) = FM(4:6) + cross_product(aero_ref_location, FM(1:3))

  return
end subroutine pseudo_aero

pure function cross_product(r, f) result(m)
  implicit none
  real, intent(in) :: r(3), f(3)
  real :: m(3)
  m(1) = r(2)*f(3) - r(3)*f(2)
  m(2) = r(3)*f(1) - r(1)*f(3)
  m(3) = r(1)*f(2) - r(2)*f(1)
end function cross_product

subroutine mass_inertia()
  implicit none
  real :: denom
  call jsonx_get(j_main, 'vehicle.mass.weight[lbf]', mass)
  mass = mass / gravity_f(0.0)

  call jsonx_get(j_main, 'vehicle.mass.Ixx[slug-ft^2]', I(1,1))
  call jsonx_get(j_main, 'vehicle.mass.Iyy[slug-ft^2]', I(2,2))
  call jsonx_get(j_main, 'vehicle.mass.Izz[slug-ft^2]', I(3,3))
  call jsonx_get(j_main, 'vehicle.mass.Ixy[slug-ft^2]', I(1,2))
  call jsonx_get(j_main, 'vehicle.mass.Ixz[slug-ft^2]', I(1,3))
  call jsonx_get(j_main, 'vehicle.mass.Iyz[slug-ft^2]', I(2,3))

  I(1,2) = -I(1,2); I(2,1) = I(1,2)
  I(1,3) = -I(1,3); I(3,1) = I(1,3)
  I(2,3) = -I(2,3); I(3,2) = I(2,3)

  denom = I(1,1)*I(2,2)*I(3,3) + I(1,2)*I(2,3)*I(3,1) + I(1,3)*I(2,1)*I(3,2) &
        - I(1,3)*I(2,2)*I(3,1) - I(1,2)*I(2,1)*I(3,3) - I(1,1)*I(2,3)*I(3,2)

  Iinv(1,1) = I(2,2)*I(3,3) - I(2,3)*I(3,2)
  Iinv(1,2) = I(1,3)*I(3,2) - I(1,2)*I(3,3)
  Iinv(1,3) = I(1,2)*I(2,3) - I(1,3)*I(2,2)
  Iinv(2,1) = I(2,3)*I(3,1) - I(2,1)*I(3,3)
  Iinv(2,2) = I(1,1)*I(3,3) - I(1,3)*I(3,1)
  Iinv(2,3) = I(1,3)*I(2,1) - I(1,1)*I(2,3)
  Iinv(3,1) = I(2,1)*I(3,2) - I(2,2)*I(3,1)
  Iinv(3,2) = I(1,2)*I(3,1) - I(1,1)*I(3,2)
  Iinv(3,3) = I(1,1)*I(2,2) - I(1,2)*I(2,1)   
  Iinv = Iinv / denom   

  call jsonx_get(j_main, 'vehicle.mass.hx[slug-ft^2/s]', h(1))
  call jsonx_get(j_main, 'vehicle.mass.hy[slug-ft^2/s]', h(2))
  call jsonx_get(j_main, 'vehicle.mass.hz[slug-ft^2/s]', h(3))

  write(*,*)"I=",I
  write(*,*) "Iinv = ",Iinv

  return

end subroutine mass_inertia

!helper function
subroutine get_vec3_from_json(j, basekey, v)
  use json_m, only: json_value, json_get, json_value_get, json_failed, json_clear_exceptions, json_value_count
  implicit none
  type(json_value), pointer, intent(in) :: j
  character(*), intent(in) :: basekey
  real, intent(inout) :: v(3)
  type(json_value), pointer :: array_ptr, element_ptr
  logical :: found
  integer :: i, n_children

  ! Initialize output to zero
  v = 0.0

  ! Get pointer to the array object
  call json_get(j, basekey, array_ptr, found)
  if (.not. found .or. json_failed()) then
    call json_clear_exceptions()
    return
  end if

  ! Get number of elements in array
  n_children = json_value_count(array_ptr)

  ! Read each element (up to 3)
  do i = 1, min(n_children, 3)
    call json_value_get(array_ptr, i, element_ptr)
    if (associated(element_ptr)) then
      call json_get(element_ptr, value=v(i), found=found)
      if (.not. found .or. json_failed()) then
        call json_clear_exceptions()
        v(i) = 0.0
      end if
    end if
  end do

end subroutine get_vec3_from_json

subroutine solve_trim_sct()
  implicit none
  real :: x(6)  ! All 6 unknowns: alpha, beta, aileron, elevator, rudder, throttle
  real :: Res(6), Jac(6,6), dx(6), R_plus(6), R_minus(6)
  real :: eps, tol, relax, max_error
  real :: alpha, beta, u, v, w, p, q, r  ! State variables
  real :: phi, theta, cphi, sphi, ctheta, stheta, ttheta  ! Euler angles and trig
  real :: Omega, g_val, denom, numerator  ! For rotation rate calculations
  integer :: iter, max_iter, i, j
  logical :: converged

  ! Default solver settings
  eps    = 0.01
  relax  = 0.9
  tol    = 1.0e-13  
  max_iter = 1000
  trim_verbose = .true.

  ! Read solver parameters from JSON 
  call jsonx_get(j_main, 'vehicle.initial.trim.solver.finite_difference_step_size', eps, 0.01)
  call jsonx_get(j_main, 'vehicle.initial.trim.solver.relaxation_factor',           relax, 0.9)
  call jsonx_get(j_main, 'vehicle.initial.trim.solver.tolerance',                   tol, 1.0e-13)
  call jsonx_get(j_main, 'vehicle.initial.trim.solver.max_iterations',              max_iter, 1000)
  call jsonx_get(j_main, 'vehicle.initial.trim.solver.verbose',                     trim_verbose, .false.)

  ! Initial guess from JSON (all 6 variables)
  x = 0.0
  ! call jsonx_get(j_main, 'vehicle.initial.trim.guess.alpha[deg]',    x(1)); x(1) = x(1) * PI / 180.0
  ! call jsonx_get(j_main, 'vehicle.initial.trim.guess.beta[deg]',     x(2)); x(2) = x(2) * PI / 180.0
  ! call jsonx_get(j_main, 'vehicle.initial.trim.guess.aileron[deg]',  x(3)); x(3) = x(3) * PI / 180.0
  ! call jsonx_get(j_main, 'vehicle.initial.trim.guess.elevator[deg]', x(4)); x(4) = x(4) * PI / 180.0
  ! call jsonx_get(j_main, 'vehicle.initial.trim.guess.rudder[deg]',   x(5)); x(5) = x(5) * PI / 180.0
  ! call jsonx_get(j_main, 'vehicle.initial.trim.guess.throttle',      x(6))
  write(*,'(A)') 'Trimming Aircraft for sct'

  ! Check if climb angle is specified 
  use_climb_angle = .false.
  if (json_key_exists(j_main, 'vehicle.initial.trim.climb_angle[deg]')) then
    use_climb_angle = .true.
    call jsonx_get(j_main, 'vehicle.initial.trim.climb_angle[deg]', trim_gamma)
    trim_gamma = trim_gamma * PI / 180.0
    write(*,'(A,F15.6)') '  --> Climb angle set to gamma [deg] = ', trim_gamma*180.0/PI
  else
    write(*,'(A,F15.6)') '  --> Elevation angle set to theta [deg] = ', trim_theta*180.0/PI
  end if

  write(*,'(A,F15.6)') '  --> Azimuth angle set to psi [deg] = ', 0.0

  if (.not. use_climb_angle) then
    write(*,'(A,F15.6)') '  --> Elevation angle set to theta [deg] = ', trim_theta*180.0/PI
  end if

  write(*,'(A,F15.6)') '  --> Bank angle set to phi [deg] = ', trim_phi*180.0/PI
  write(*,'(A)') ''
  write(*,'(A,ES24.15)') 'Initial theta [deg] = ', trim_theta*180.0/PI
  if (use_climb_angle) then
    write(*,'(A,ES24.15)') 'Initial gamma [deg] = ', trim_gamma*180.0/PI
  else
    write(*,'(A,ES24.15)') 'Initial gamma [deg] = ', 0.0
  end if
  write(*,'(A,ES24.15)') 'Initial phi [deg]   = ', trim_phi*180.0/PI
  write(*,'(A,ES24.15)') 'Initial beta [deg]  = ', x(2)*180.0/PI
  write(*,'(A)') ''
  
  write(*,'(A)') 'Newton Solver Settings:'
  write(*,'(A,ES24.15)') 'Finite Difference Step Size = ', eps
  write(*,'(A,ES24.15)') '          Relaxation Factor = ', relax
  write(*,'(A,ES24.15)') '                  Tolerance = ', tol
  ! Newton's method iteration loop
  do iter = 1, max_iter
    ! If using climb angle, compute elevation angle from current state
    if (use_climb_angle) then
      ! Compute u, v, w from current alpha, beta
      alpha = x(1)
      beta = x(2)
      u = trim_V * cos(alpha) * cos(beta)
      v = trim_V * sin(beta)
      w = trim_V * sin(alpha) * cos(beta)
      
      ! Compute theta from gamma using Eqs. 7.2.8 and 7.2.5
      trim_theta = compute_theta_from_gamma(trim_gamma, u, v, w, trim_phi)
      
      if (trim_verbose) then
        block
          real :: theta_plus, theta_minus, Stheta_plus, Stheta_minus
          real :: Sphi, Cphi, Sgamma, numerator, denom, discriminant, Vstotal
          
          write(*,'(A)') ''
          write(*,'(A)') 'Solving for elevation angle given a climb angle:'
          
          Vstotal = sqrt(u**2 + v**2 + w**2)
          Sphi = sin(trim_phi)
          Cphi = cos(trim_phi)
          Sgamma = sin(trim_gamma)
          numerator = u * Vstotal * Sgamma
          denom = u**2 + (v*Sphi + w*Cphi)**2
          discriminant = u**2 + (v*Sphi + w*Cphi)**2 - Vstotal**2 * Sgamma**2
          
          if (discriminant >= 0.0 .and. abs(denom) > 1.0e-10) then
            Stheta_plus  = (numerator + (v*Sphi + w*Cphi)*sqrt(discriminant)) / denom
            Stheta_minus = (numerator - (v*Sphi + w*Cphi)*sqrt(discriminant)) / denom
            
            if (abs(Stheta_plus) <= 1.0) then
              theta_plus = asin(Stheta_plus)
              write(*,'(A,ES24.15)') '        theta 1 [deg] = ', theta_plus * 180.0 / PI
            end if
            if (abs(Stheta_minus) <= 1.0) then
              theta_minus = asin(Stheta_minus)
              write(*,'(A,ES24.15)') '        theta 2 [deg] = ', theta_minus * 180.0 / PI
            end if
          end if
          
          write(*,'(A,ES24.15)') '  Correct theta [deg] = ', trim_theta * 180.0 / PI
          write(*,'(A,ES24.15)') '  Correct theta [rad] = ', trim_theta
        end block
      end if
    end if
    
    ! Update rotation rates based on current state
    
    if (trim_verbose) then
      write(*,'(A)') ''
      write(*,'(A)') 'Updating rotation rates for sct'
      
      ! Computing rotation rates using kinematic constraints for SCT
      
      alpha = x(1)
      beta = x(2)
      
      phi = trim_phi
      write(*,*) "phi =", phi
      theta = trim_theta
      write(*,* ) "theta =", theta
      cphi = cos(phi)
      sphi = sin(phi)
      ctheta = cos(theta)
      stheta = sin(theta)
      
      ! Compute u, v, w from alpha, beta
      u = trim_V * cos(alpha) * cos(beta)
      v = trim_V * sin(beta)
      w = trim_V * sin(alpha) * cos(beta)
      
      if (abs(cphi) < 1.0e-10) then
        p = 0.0
        q = 0.0
        r = 0.0
      else
        ! Equation (6.2.3): rotation rates for steady-coordinated turn

        g_val = gravity_f(trim_alt)
        denom = u * ctheta * cphi + w * stheta
        if (abs(denom) < 1.0e-10) then
          p = 0.0
          q = 0.0
          r = 0.0
        else
          numerator = g_val * sphi * ctheta / denom
          p = numerator * (-stheta)
          q = numerator * (sphi * ctheta)
          r = numerator * (cphi * ctheta)
        end if
      end if
      write(*,'(A,ES24.15)') 'altitude = ', trim_alt
      write(*,'(A,ES24.15)') 'gravity = ', g_val
      write(*,'(A,ES24.15)') 'u  = ', u 
      write(*,'(A,ES24.15)') 'v  = ', v
      write(*,'(A,ES24.15)') 'w  = ', w
      write(*,'(A,ES24.15)') 'p [deg/s] = ', p * 180.0 / PI
      write(*,'(A,ES24.15)') 'q [deg/s] = ', q * 180.0 / PI
      write(*,'(A,ES24.15)') 'r [deg/s] = ', r * 180.0 / PI
    end if
    
    ! Compute residual R(x) - all 6 equations
    Res = calc_residual_sct(x,p,q,r)
    
    max_error = maxval(abs(Res))
    
    if (trim_verbose) then
      write(*,'(A)') ''
      write(*,'(A)') 'G defined as G = [alpha, beta, aileron, elevator, rudder, throttle]'
      write(*,'(A,6(ES24.15,A))') '      G = ', x(1), ', ', x(2), ', ', x(3), ', ', x(4), ', ', x(5), ', ', x(6), ','
      write(*,'(A,6(ES24.15,A))') '      R = ', Res(1), ', ', Res(2), ', ', Res(3), ', ', Res(4), ', ', Res(5), ', ', Res(6), ','
    end if
    
    ! Check convergence
    if (max_error < tol) then
      converged = .true.
      exit
    end if
    
    if (trim_verbose) then
      write(*,'(A)') ''
      write(*,'(A)') 'Building Jacobian Matrix:'
    end if
    
    ! Build Jacobian using central finite differences (6x6)
    do j = 1, 6
      if (trim_verbose) then
        write(*,'(A)') ''
        write(*,'(A,I2,A)') 'Computing gradient relative to G[', j-1, ' ]'
      end if
      
      ! Central difference (for all variables including throttle)
      ! Forward perturbation
      x(j) = x(j) + eps
      R_plus = calc_residual_sct(x,p,q,r)
      
      if (trim_verbose) then
        write(*,'(A)') '   Positive Finite-Difference Step '
        write(*,'(A,6(ES24.15,A))') '      G = ', x(1), ', ', x(2), ', ', x(3), ', ', x(4), ', ', x(5), ', ', x(6), ','
        write(*,'(A,6(ES24.15,A))') '      R = ', R_plus(1), ', ', R_plus(2), ', ', R_plus(3), ', ', R_plus(4), ', ', R_plus(5), ', ', R_plus(6), ','
      end if
      
      ! Backward perturbation
      x(j) = x(j) - 2.0*eps
      R_minus = calc_residual_sct(x,p,q,r)
      
      if (trim_verbose) then
        write(*,'(A)') '   Negative Finite-Difference Step '
        write(*,'(A,6(ES24.15,A))') '      G = ', x(1), ', ', x(2), ', ', x(3), ', ', x(4), ', ', x(5), ', ', x(6), ','
        write(*,'(A,6(ES24.15,A))') '      R = ', R_minus(1), ', ', R_minus(2), ', ', R_minus(3), ', ', R_minus(4), ', ', R_minus(5), ', ', R_minus(6), ','
      end if
      
      ! Restore
      x(j) = x(j) + eps
      
      ! Central difference for all 6 rows
      do i = 1, 6
        Jac(i,j) = (R_plus(i) - R_minus(i)) / (2.0*eps)
      end do
    end do
    
    if (trim_verbose) then
      write(*,'(A)') ''
      write(*,'(A)') 'Jacobian Matrix ='
      do i = 1, 6
        write(*,'(2X,6(ES24.15,2X))') (Jac(i,j), j=1,6)
      end do
    end if

    ! Solve J * dx = -R using LU decomposition
    call solve_linear_system(Jac, -Res, dx)
    
    if (trim_verbose) then
      write(*,'(A)') ''
      write(*,'(A,6(ES24.15,A))') 'Delta G = ', dx(1), ', ', dx(2), ', ', dx(3), ', ', dx(4), ', ', dx(5), ', ', dx(6), ','
    end if
    
    x = x + relax * dx
    
    if (x(6) < 0.0) then
      x(6) = 0.0
    end if
    
    if (trim_verbose) then
      write(*,'(A)') 'New G : '
      ! Recompute residual after update to show new state
      Res = calc_residual_sct(x,p,q,r)
      write(*,'(A,6(ES24.15,A))') '      G = ', x(1), ', ', x(2), ', ', x(3), ', ', x(4), ', ', x(5), ', ', x(6), ','
      write(*,'(A,6(ES24.15,A))') '      R = ', Res(1), ', ', Res(2), ', ', Res(3), ', ', Res(4), ', ', Res(5), ', ', Res(6), ','
      write(*,'(A)') ''
      write(*,'(A,I8,3X,12(A,ES24.15,3X))') &
        'Iteration', iter, 'Residual', maxval(abs(Res)), &
        'alpha[deg]', x(1)*180.0/PI, 'beta[deg]', x(2)*180.0/PI, &
        'p[deg/s]', p*180.0/PI, 'q[deg/s]', q*180.0/PI, 'r[deg/s]', r*180.0/PI, &
        'phi[deg]', trim_phi*180.0/PI, 'theta[deg]', trim_theta*180.0/PI, &
        'aileron[deg]', x(3)*180.0/PI, 'elevator[deg]', x(4)*180.0/PI, &
        'rudder[deg]', x(5)*180.0/PI, 'throttle[]', x(6)
    end if
  end do

  if (.not. converged) then
    write(*,'(A,I5,A)') 'Warning: Did not converge after ', max_iter, ' iterations'
    write(*,'(A,ES12.5)') 'Final max_error = ', max_error
  end if

  ! Set initial state from converged (or final) trim solution
  call set_initial_state_from_trim(x)

  ! Print trim results
  call print_trim_results(x)

end subroutine solve_trim_sct

function calc_residual_sct(x, p_in, q_in, r_in) result(Res)
  implicit none

  real, intent(in) :: x(6)      ! [alpha, beta, da, de, dr, tau]
  real, intent(in) :: p_in, q_in, r_in   
  real :: Res(6)
  real :: y(13), dydt(13)
  real :: alpha, beta, u, v, w
  real :: phi, theta, psi
  real :: eul(3)
  real :: p, q, r  
  real :: cphi, sphi, ctheta, stheta
  real :: g_val, denom, numerator
  real :: matmul_h_pqr(3), I_stuff(3)

  ! Extract G
  alpha = x(1)
  beta  = x(2)
  controls(1) = x(3)  ! aileron
  controls(2) = x(4)  ! elevator
  controls(3) = x(5)  ! rudder
  controls(4) = x(6)  ! throttle

  ! Kinematics: body-axis velocities from V, alpha, beta
  u = trim_V * cos(alpha) * cos(beta)
  v = trim_V * sin(beta)
  w = trim_V * sin(alpha) * cos(beta)

  ! Euler angles (phi, theta specified; psi arbitrary)
  phi   = trim_phi
  theta = trim_theta
  psi   = 0.0
  
  ! Computing rotation rates based on SCT kinematic constraints
  cphi = cos(phi)
  sphi = sin(phi)
  ctheta = cos(theta)
  stheta = sin(theta)
  
  if (abs(cphi) < 1.0e-10) then
    ! Special case: phi = 90 degrees (wings vertical)
    p = 0.0
    q = 0.0
    r = 0.0
  else
    g_val = gravity_f(trim_alt)
    denom = u * ctheta * cphi + w * stheta
    if (abs(denom) < 1.0e-10) then
      p = 0.0
      q = 0.0
      r = 0.0
    else
      numerator = g_val * sphi * ctheta / denom
      p = numerator * (-stheta)
      q = numerator * (sphi * ctheta)
      r = numerator * (cphi * ctheta)
    end if
  end if

  ! Build state vector y
  y = 0.0
  y(1) = u
  y(2) = v
  y(3) = w
  y(4) = p
  y(5) = q
  y(6) = r
  y(7) = 0.0
  y(8) = 0.0
  y(9) = -trim_alt

  eul = (/ phi, theta, psi /)
  y(10:13) = euler_to_quat(eul)

  ! Debug output for first call only (when all x values are zero)
  if (trim_verbose .and. abs(alpha) < 1.0e-10 .and. abs(beta) < 1.0e-10 .and. &
      abs(controls(1)) < 1.0e-10 .and. abs(controls(2)) < 1.0e-10) then
    write(*,'(A)') ''
    write(*,'(A,13(ES24.15,2X))') 'y going into diff_eq in calc_R  ', y
    
    ! Compute matmul [h][pqr] for debug
    matmul_h_pqr(1) = h(2)*r - h(3)*q
    matmul_h_pqr(2) = h(3)*p - h(1)*r  
    matmul_h_pqr(3) = h(1)*q - h(2)*p
    write(*,'(A,3(ES24.15,2X))') '  matmul [h][pqr] = ', matmul_h_pqr
    
    ! Compute I stuff for debug
    I_stuff(1) = (I(2,2) - I(3,3))*q*r + I(2,3)*(q**2 - r**2) + I(1,3)*p*q - I(1,2)*p*r
    I_stuff(2) = (I(3,3) - I(1,1))*p*r + I(1,3)*(r**2 - p**2) + I(1,2)*q*r - I(2,3)*p*q
    I_stuff(3) = (I(1,1) - I(2,2))*p*q + I(1,2)*(p**2 - q**2) + I(2,3)*p*r - I(1,3)*q*r
    write(*,'(A,ES24.15)') ' I stuff(1) = ', I_stuff(1)
    write(*,'(A,ES24.15)') ' I stuff(2) = ', I_stuff(2)
    write(*,'(A,ES24.15)') ' I stuff(3) = ', I_stuff(3)
    write(*,'(A)') ''
  end if

  ! Dynamics and residual
  dydt = diff_eq(0.0, y, 1)
  Res  = dydt(1:6)
  
  ! Debug output for dydt
  if (trim_verbose .and. abs(alpha) < 1.0e-10 .and. abs(beta) < 1.0e-10 .and. &
      abs(controls(1)) < 1.0e-10 .and. abs(controls(2)) < 1.0e-10) then
    write(*,'(A,13(ES24.15,2X))') 'dydt = ', dydt
    write(*,'(A,6(ES24.15,3X))') '   R_ = ', Res
  end if
end function calc_residual_sct

subroutine set_initial_state_from_trim(x)
  implicit none
  real, intent(in) :: x(6)
  real :: alpha, beta, u, v, w, p, q, r, phi, theta, psi
  real :: cphi, sphi, ctheta, stheta
  real :: Omega, g_val, denom, numerator

  alpha = x(1)
  beta = x(2)
  controls(1) = x(3)  ! aileron
  controls(2) = x(4)  ! elevator
  controls(3) = x(5)  ! rudder
  controls(4) = x(6)  ! throttle

  ! Compute velocities
  u = trim_V * cos(alpha) * cos(beta)
  v = trim_V * sin(beta)
  w = trim_V * sin(alpha) * cos(beta)

  ! Euler angles
  phi = trim_phi
  theta = trim_theta
  psi = 0.0

  cphi = cos(phi)
  sphi = sin(phi)
  ctheta = cos(theta)
  stheta = sin(theta)

    g_val = gravity_f(trim_alt)
    
    if (abs(cphi) < 1.0e-10) then
      p = 0.0
      q = 0.0
      r = 0.0
    else
      ! Equation (6.2.3): rotation rates for steady-coordinated turn
      denom = u * ctheta * cphi + w * stheta
      if (abs(denom) < 1.0e-10) then
        p = 0.0
        q = 0.0
        r = 0.0
      else
        numerator = g_val * sphi * ctheta / denom
        p = numerator * (-stheta)
        q = numerator * (sphi * ctheta)
        r = numerator * (cphi * ctheta)
      end if
    end if

  ! Set initial state
  y_init(1) = u
  y_init(2) = v
  y_init(3) = w
  y_init(4) = p
  y_init(5) = q
  y_init(6) = r
  y_init(7) = 0.0
  y_init(8) = 0.0
  y_init(9) = -trim_alt
  y_init(10:13) = euler_to_quat([phi, theta, psi])

end subroutine set_initial_state_from_trim

subroutine print_trim_results(x)
  implicit none
  real, intent(in) :: x(6)
  real :: alpha_deg, beta_deg, da_deg, de_deg, dr_deg
  real :: p_deg, q_deg, r_deg, phi_deg, theta_deg
  real :: y(13), V, dydt(13)

  alpha_deg = x(1) * 180.0 / PI
  beta_deg = x(2) * 180.0 / PI
  da_deg = x(3) * 180.0 / PI
  de_deg = x(4) * 180.0 / PI
  dr_deg = x(5) * 180.0 / PI

  phi_deg = trim_phi * 180.0 / PI
  theta_deg = trim_theta * 180.0 / PI

  p_deg = y_init(4) * 180.0 / PI
  q_deg = y_init(5) * 180.0 / PI
  r_deg = y_init(6) * 180.0 / PI

  V = sqrt(y_init(1)**2 + y_init(2)**2 + y_init(3)**2)

  write(*,'(A)') ''
  write(*,'(A)') '============================================'
  write(*,'(A)') 'Trim Solution for Steady Coordinated Turn'
  write(*,'(A)') '============================================'
  write(*,'(A,F10.2,A)') 'Velocity: ', V, ' ft/s'
  write(*,'(A,F10.2,A)') 'Altitude: ', trim_alt, ' ft'
  write(*,'(A,F16.13,A)') 'Bank angle (phi): ', phi_deg, ' deg'
  write(*,'(A,F16.13,A)') 'Elevation angle (theta): ', theta_deg, ' deg'
  write(*,'(A)') '--------------------------------------------'
  write(*,'(A,F16.13,A)') 'alpha = ', alpha_deg, ' deg'
  write(*,'(A,F16.13,A)') 'beta  = ', beta_deg, ' deg'
  write(*,'(A,F16.13,A)') 'p     = ', p_deg, ' deg/s'
  write(*,'(A,F16.13,A)') 'q     = ', q_deg, ' deg/s'
  write(*,'(A,F17.13,A)') 'r     = ', r_deg, ' deg/s'
  write(*,'(A,F17.13,A)') 'aileron  = ', da_deg, ' deg'
  write(*,'(A,F17.13,A)') 'elevator = ', de_deg, ' deg'
  write(*,'(A,F17.13,A)') 'rudder   = ', dr_deg, ' deg'
  write(*,'(A,F17.14)') 'throttle = ', x(6)
  write(*,'(A)') '============================================'

  ! Verify trim by computing residual
  y = y_init
  dydt = diff_eq(0.0, y, 1)
  write(*,'(A)') 'Residual verification (should be near zero):'
  write(*,'(A,6ES14.6)') 'dydt(1:6) = ', dydt(1:6)
  write(*,'(A)') ''

  if (print_state_aero) then
    call print_state_and_aero()
  end if

end subroutine print_trim_results

subroutine solve_trim_shss()
  use json_m, only: json_get, json_failed, json_clear_exceptions
  implicit none
  real :: x(6) 
  real :: Res(6), Jac(6,6), dx(6), R_plus(6), R_minus(6)
  real :: eps, tol, relax, max_error
  integer :: iter, max_iter, i, j
  logical :: converged
  logical :: sideslip_specified
  real :: specified_beta
  logical :: found

  ! Default solver settings
  eps    = 0.01
  relax  = 0.9
  tol    = 1.0e-6  
  max_iter = 1000
  trim_verbose = .false.

  ! Read solver parameters from JSON
  call jsonx_get(j_main, 'vehicle.initial.trim.solver.finite_difference_step_size', eps, 0.01)
  call jsonx_get(j_main, 'vehicle.initial.trim.solver.relaxation_factor',           relax, 0.9)
  call jsonx_get(j_main, 'vehicle.initial.trim.solver.tolerance',                   tol, 1.0e-13)
  call jsonx_get(j_main, 'vehicle.initial.trim.solver.max_iterations',              max_iter, 1000)
  call jsonx_get(j_main, 'vehicle.initial.trim.solver.verbose',                     trim_verbose, .false.)

  ! Check if sideslip angle is specified
  sideslip_specified = .false.
  specified_beta = 0.0
    
  call json_get(j_main, 'vehicle.initial.trim.sideslip[deg]', specified_beta, found)
  if (found .and. .not. json_failed()) then
    sideslip_specified = .true.
    specified_beta = specified_beta * PI / 180.0
  end if
  call json_clear_exceptions()

  ! Initial guess from JSON (all 6 variables)
  x = 0.0
  call jsonx_get(j_main, 'vehicle.initial.trim.guess.alpha[deg]',    x(1)); x(1) = x(1) * PI / 180.0

  if (sideslip_specified) then
    ! x(2) is phi (bank angle) - use specified phi as initial guess
    call jsonx_get(j_main, 'vehicle.initial.trim.guess.phi[deg]',  x(2)); x(2) = x(2) * PI / 180.0
  else
    ! x(2) is beta (sideslip angle) - use specified beta as initial guess
    call jsonx_get(j_main, 'vehicle.initial.trim.guess.beta[deg]', x(2)); x(2) = x(2) * PI / 180.0
  end if

  call jsonx_get(j_main, 'vehicle.initial.trim.guess.aileron[deg]',  x(3)); x(3) = x(3) * PI / 180.0
  call jsonx_get(j_main, 'vehicle.initial.trim.guess.elevator[deg]', x(4)); x(4) = x(4) * PI / 180.0
  call jsonx_get(j_main, 'vehicle.initial.trim.guess.rudder[deg]',   x(5)); x(5) = x(5) * PI / 180.0
  call jsonx_get(j_main, 'vehicle.initial.trim.guess.throttle',      x(6))

  write(*,'(A)') 'Trimming Aircraft for shss'

  ! Check if climb angle is specified instead of using the elevation angle directly
  use_climb_angle = .false.
  if (json_key_exists(j_main, 'vehicle.initial.trim.climb_angle[deg]')) then
    use_climb_angle = .true.
    call jsonx_get(j_main, 'vehicle.initial.trim.climb_angle[deg]', trim_gamma)
    trim_gamma = trim_gamma * PI / 180.0
    write(*,'(A,F15.6)') '  --> Climb angle set to gamma [deg] = ', trim_gamma*180.0/PI
  else
    write(*,'(A,F15.6)') '  --> Elevation angle set to theta [deg] = ', trim_theta*180.0/PI
  end if

  write(*,'(A,F15.6)') '  --> Azimuth angle set to psi [deg] = ', 0.0

  if (.not. use_climb_angle) then
    write(*,'(A,F15.6)') '  --> Elevation angle set to theta [deg] = ', trim_theta*180.0/PI
  end if

  if (sideslip_specified) then
    write(*,'(A,F15.6)') '  --> Sideslip angle set to beta [deg] = ', specified_beta*180.0/PI
  else
    write(*,'(A,F15.6)') '  --> Bank angle set to phi [deg] = ', trim_phi*180.0/PI
  end if

  write(*,'(A)') ''
  write(*,'(A,ES24.15)') 'Initial theta [deg] = ', trim_theta*180.0/PI
  if (use_climb_angle) then
    write(*,'(A,ES24.15)') 'Initial gamma [deg] = ', trim_gamma*180.0/PI
  else
    write(*,'(A,ES24.15)') 'Initial gamma [deg] = ', 0.0
  end if

  if (sideslip_specified) then
    write(*,'(A,ES24.15)') 'Initial phi [deg]   = ', x(2)*180.0/PI
    write(*,'(A,ES24.15)') 'Initial beta [deg]  = ', specified_beta*180.0/PI
  else
    write(*,'(A,ES24.15)') 'Initial phi [deg]   = ', trim_phi*180.0/PI
    write(*,'(A,ES24.15)') 'Initial beta [deg]  = ', x(2)*180.0/PI
  end if

  write(*,'(A)') ''
  write(*,'(A)') 'Newton Solver Settings:'
  write(*,'(A,ES24.15)') 'Finite Difference Step Size = ', eps
  write(*,'(A,ES24.15)') '          Relaxation Factor = ', relax
  write(*,'(A,ES24.15)') '                  Tolerance = ', tol
  write(*,'(A)') ''

  converged = .false.

  ! Newton's method iteration loop
  do iter = 1, max_iter
    ! If using climb angle, compute elevation angle from current state
    if (use_climb_angle) then
      block
        real :: phi_for_calc, alpha_temp, beta_temp, u_temp, v_temp, w_temp
        
        ! Compute u, v, w from current alpha, beta (or phi if beta is specified)
        alpha_temp = x(1)
        if (sideslip_specified) then
          beta_temp = specified_beta
          phi_for_calc = x(2)  ! x(2) is phi when beta is specified
        else
          beta_temp = x(2)
          phi_for_calc = trim_phi  ! otherwise use the specified phi
        end if
        
        u_temp = trim_V * cos(alpha_temp) * cos(beta_temp)
        v_temp = trim_V * sin(beta_temp)
        w_temp = trim_V * sin(alpha_temp) * cos(beta_temp)
        
        ! Compute theta from gamma using Eqs. 7.2.8 and 7.2.5
        trim_theta = compute_theta_from_gamma(trim_gamma, u_temp, v_temp, w_temp, phi_for_calc)
        
        if (trim_verbose) then
          block
            real :: theta_plus, theta_minus, Stheta_plus, Stheta_minus
            real :: Sphi, Cphi, Sgamma, numerator, denom, discriminant, Vstotal
            
            write(*,'(A)') ''
            write(*,'(A)') 'Solving for elevation angle given a climb angle:'
            
            Vstotal = sqrt(u_temp**2 + v_temp**2 + w_temp**2)
            Sphi = sin(phi_for_calc)
            Cphi = cos(phi_for_calc)
            Sgamma = sin(trim_gamma)
            numerator = u_temp * Vstotal * Sgamma
            denom = u_temp**2 + (v_temp*Sphi + w_temp*Cphi)**2
            discriminant = u_temp**2 + (v_temp*Sphi + w_temp*Cphi)**2 - Vstotal**2 * Sgamma**2
            
            if (discriminant >= 0.0 .and. abs(denom) > 1.0e-10) then
              Stheta_plus  = (numerator + (v_temp*Sphi + w_temp*Cphi)*sqrt(discriminant)) / denom
              Stheta_minus = (numerator - (v_temp*Sphi + w_temp*Cphi)*sqrt(discriminant)) / denom
              
              if (abs(Stheta_plus) <= 1.0) then
                theta_plus = asin(Stheta_plus)
                write(*,'(A,ES24.15)') '        theta 1 [deg] = ', theta_plus * 180.0 / PI
              end if
              if (abs(Stheta_minus) <= 1.0) then
                theta_minus = asin(Stheta_minus)
                write(*,'(A,ES24.15)') '        theta 2 [deg] = ', theta_minus * 180.0 / PI
              end if
            end if
            
            write(*,'(A,ES24.15)') '  Correct theta [deg] = ', trim_theta * 180.0 / PI
            write(*,'(A,ES24.15)') '  Correct theta [rad] = ', trim_theta
          end block
        end if
      end block
    end if
    
    ! Compute residual R(x) - all 6 equations
    if (sideslip_specified) then
      Res = calc_residual_shss_beta_spec(x, specified_beta)
    else
      Res = calc_residual_shss(x)
    end if
    
    max_error = maxval(abs(Res))
    
    if (trim_verbose) then
      write(*,'(A)') ''
      if (sideslip_specified) then
        write(*,'(A)') 'G defined as G = [alpha, phi, aileron, elevator, rudder, throttle]'
      else
        write(*,'(A)') 'G defined as G = [alpha, beta, aileron, elevator, rudder, throttle]'
      end if
      write(*,'(A,6(ES24.15,A))') '      G = ', x(1), ', ', x(2), ', ', x(3), ', ', x(4), ', ', x(5), ', ', x(6), ','
      write(*,'(A,6(ES24.15,A))') '      R = ', Res(1), ', ', Res(2), ', ', Res(3), ', ', Res(4), ', ', Res(5), ', ', Res(6), ','
    end if
    
    ! Check convergence
    if (max_error < tol) then
      converged = .true.
      exit
    end if
    
    if (trim_verbose) then
      write(*,'(A)') ''
      write(*,'(A)') 'Building Jacobian Matrix:'
    end if
    
    ! Build Jacobian using central finite differences (6x6)
    do j = 1, 6
      if (trim_verbose) then
        write(*,'(A)') ''
        write(*,'(A,I2,A)') 'Computing gradient relative to G[', j-1, ' ]'
      end if
      
      ! Central difference (for all variables including throttle)
      ! Forward perturbation
      x(j) = x(j) + eps
      if (sideslip_specified) then
        R_plus = calc_residual_shss_beta_spec(x, specified_beta)
      else
        R_plus = calc_residual_shss(x)
      end if
      
      if (trim_verbose) then
        write(*,'(A)') '   Positive Finite-Difference Step '
        write(*,'(A,6(ES24.15,A))') '      G = ', x(1), ', ', x(2), ', ', x(3), ', ', x(4), ', ', x(5), ', ', x(6), ','
        write(*,'(A,6(ES24.15,A))') '      R = ', R_plus(1), ', ', R_plus(2), ', ', R_plus(3), ', ', R_plus(4), ', ', R_plus(5), ', ', R_plus(6), ','
      end if
      
      ! Backward perturbation
      x(j) = x(j) - 2.0*eps
      if (sideslip_specified) then
        R_minus = calc_residual_shss_beta_spec(x, specified_beta)
      else
        R_minus = calc_residual_shss(x)
      end if
      
      if (trim_verbose) then
        write(*,'(A)') '   Negative Finite-Difference Step '
        write(*,'(A,6(ES24.15,A))') '      G = ', x(1), ', ', x(2), ', ', x(3), ', ', x(4), ', ', x(5), ', ', x(6), ','
        write(*,'(A,6(ES24.15,A))') '      R = ', R_minus(1), ', ', R_minus(2), ', ', R_minus(3), ', ', R_minus(4), ', ', R_minus(5), ', ', R_minus(6), ','
      end if
      
      ! Restore
      x(j) = x(j) + eps
      
      ! Central difference for all 6 rows
      do i = 1, 6
        Jac(i,j) = (R_plus(i) - R_minus(i)) / (2.0*eps)
      end do
    end do
    
    if (trim_verbose) then
      write(*,'(A)') ''
      write(*,'(A)') 'Jacobian Matrix ='
      do i = 1, 6
        write(*,'(2X,6(ES24.15,2X))') (Jac(i,j), j=1,6)
      end do
    end if
    
    ! Solve J * dx = -R using LU decomposition
    call solve_linear_system(Jac, -Res, dx)
    
    if (trim_verbose) then
      write(*,'(A)') ''
      write(*,'(A,6(ES24.15,A))') 'Delta G = ', dx(1), ', ', dx(2), ', ', dx(3), ', ', dx(4), ', ', dx(5), ', ', dx(6), ','
    end if
    
    ! Update with relaxation: x = x + relax * dx
    x = x + relax * dx
    
    ! Clamp throttle to [0, 1] (physical constraint)
    if (x(6) < 0.0) then
      x(6) = 0.0
    end if
    if (x(6) > 1.0) then
      x(6) = 1.0
    end if
    
    if (trim_verbose) then
      write(*,'(A)') 'New G : '
      ! Recompute residual after update to show new state
      if (sideslip_specified) then
        Res = calc_residual_shss_beta_spec(x, specified_beta)
      else
        Res = calc_residual_shss(x)
      end if
      write(*,'(A,6(ES24.15,A))') '      G = ', x(1), ', ', x(2), ', ', x(3), ', ', x(4), ', ', x(5), ', ', x(6), ','
      write(*,'(A,6(ES24.15,A))') '      R = ', Res(1), ', ', Res(2), ', ', Res(3), ', ', Res(4), ', ', Res(5), ', ', Res(6), ','
      write(*,'(A)') ''
      if (sideslip_specified) then
        write(*,'(A,I8,3X,12(A,ES24.15,3X))') &
          'Iteration', iter, 'Residual', maxval(abs(Res)), &
          'alpha[deg]', x(1)*180.0/PI, 'beta[deg]', specified_beta*180.0/PI, &
          'p[deg/s]', 0.0, 'q[deg/s]', 0.0, 'r[deg/s]', 0.0, &
          'phi[deg]', x(2)*180.0/PI, 'theta[deg]', trim_theta*180.0/PI, &
          'aileron[deg]', x(3)*180.0/PI, 'elevator[deg]', x(4)*180.0/PI, &
          'rudder[deg]', x(5)*180.0/PI, 'throttle[]', x(6)
      else
        write(*,'(A,I8,3X,12(A,ES24.15,3X))') &
          'Iteration', iter, 'Residual', maxval(abs(Res)), &
          'alpha[deg]', x(1)*180.0/PI, 'beta[deg]', x(2)*180.0/PI, &
          'p[deg/s]', 0.0, 'q[deg/s]', 0.0, 'r[deg/s]', 0.0, &
          'phi[deg]', trim_phi*180.0/PI, 'theta[deg]', trim_theta*180.0/PI, &
          'aileron[deg]', x(3)*180.0/PI, 'elevator[deg]', x(4)*180.0/PI, &
          'rudder[deg]', x(5)*180.0/PI, 'throttle[]', x(6)
      end if
    end if
  end do

  if (.not. converged) then
    write(*,'(A,I5,A)') 'Warning: Did not converge after ', max_iter, ' iterations'
    write(*,'(A,ES12.5)') 'Final max_error = ', max_error
  end if

  ! Set initial state from converged (or final) trim solution
  if (sideslip_specified) then
    call set_initial_state_from_trim_shss_beta_spec(x, specified_beta)
  else
    call set_initial_state_from_trim_shss(x)
  end if

  ! Print trim results
  if (sideslip_specified) then
    call print_trim_results_shss_beta_spec(x, specified_beta)
  else
    call print_trim_results_shss(x)
  end if

end subroutine solve_trim_shss

function calc_residual_shss(x) result(Res)
  implicit none
  real, intent(in) :: x(6)
  real :: Res(6)
  real :: y(13), dydt(13)
  real :: alpha, beta, u, v, w
  real :: phi, theta, psi

  alpha = x(1)
  beta = x(2)
  controls(1) = x(3)  ! aileron
  controls(2) = x(4)  ! elevator
  controls(3) = x(5)  ! rudder
  controls(4) = max(0.0, x(6))  ! throttle 

  ! Compute body-axis velocities from V, alpha, beta
  u = trim_V * cos(alpha) * cos(beta)
  v = trim_V * sin(beta)
  w = trim_V * sin(alpha) * cos(beta)

  ! Euler angles (phi and theta specified; psi arbitrary)
  phi = trim_phi
  theta = trim_theta
  psi = 0.0  ! arbitrary for dynamics

  ! For SHSS
  y(1) = u
  y(2) = v
  y(3) = w
  y(4) = 0.0  ! p = 0 for SHSS
  y(5) = 0.0  ! q = 0 for SHSS
  y(6) = 0.0  ! r = 0 for SHSS
  y(7) = 0.0  
  y(8) = 0.0  
  y(9) = -trim_alt  ! z position (negative altitude)

  ! Quaternion from Euler angles
  y(10:13) = euler_to_quat([phi, theta, psi])

  ! Compute equations of motion dydt = f(y)
  dydt = diff_eq(0.0, y, 1)

  ! Residual is first 6 components (should be zero at trim)
  Res = dydt(1:6)

  end function calc_residual_shss

  
  function calc_residual_shss_beta_spec(x, beta_specified) result(Res)
  implicit none
  real, intent(in) :: x(6)
  real, intent(in) :: beta_specified
  real :: Res(6)
  real :: y(13), dydt(13)
  real :: alpha, beta, u, v, w
  real :: phi, theta, psi

  ! Extract variables: x = [alpha, phi, aileron, elevator, rudder, throttle]
  alpha = x(1)
  phi = x(2)  ! phi is now an unknown
  beta = beta_specified  ! beta is specified
  controls(1) = x(3)  ! aileron
  controls(2) = x(4)  ! elevator
  controls(3) = x(5)  ! rudder
  controls(4) = max(0.0, x(6))  ! throttle (clamped to non-negative)

  ! Compute body-axis velocities from V, alpha, beta
  u = trim_V * cos(alpha) * cos(beta)
  v = trim_V * sin(beta)
  w = trim_V * sin(alpha) * cos(beta)

  ! Euler angles (theta specified; psi arbitrary)
  theta = trim_theta
  psi = 0.0  ! arbitrary for dynamics

  ! For SHSS: [p, q, r] = [0, 0, 0]

  ! Build state vector
  y(1) = u
  y(2) = v
  y(3) = w
  y(4) = 0.0  ! p = 0 for SHSS
  y(5) = 0.0  ! q = 0 for SHSS
  y(6) = 0.0  ! r = 0 for SHSS
  y(7) = 0.0  ! x position (irrelevant for dynamics)
  y(8) = 0.0  ! y position (irrelevant for dynamics)
  y(9) = -trim_alt  ! z position (negative altitude)

  ! Quaternion from Euler angles
  y(10:13) = euler_to_quat([phi, theta, psi])

  ! Compute equations of motion dydt = f(y)
  dydt = diff_eq(0.0, y, 1)

  ! Residual is first 6 components (should be zero at trim)
  Res = dydt(1:6)

end function calc_residual_shss_beta_spec

subroutine set_initial_state_from_trim_shss(x)
  implicit none
  real, intent(in) :: x(6)
  real :: alpha, beta, u, v, w, phi, theta, psi

  alpha = x(1)
  beta = x(2)
  controls(1) = x(3)  ! aileron
  controls(2) = x(4)  ! elevator
  controls(3) = x(5)  ! rudder
  controls(4) = x(6)  ! throttle

  ! Compute velocities
  u = trim_V * cos(alpha) * cos(beta)
  v = trim_V * sin(beta)
  w = trim_V * sin(alpha) * cos(beta)

  ! Euler angles
  phi = trim_phi
  theta = trim_theta
  psi = 0.0

  ! For SHSS: [p, q, r] = [0, 0, 0]

  ! Set initial state
  y_init(1) = u
  y_init(2) = v
  y_init(3) = w
  y_init(4) = 0.0  ! p = 0 for SHSS
  y_init(5) = 0.0  ! q = 0 for SHSS
  y_init(6) = 0.0  ! r = 0 for SHSS
  y_init(7) = 0.0
  y_init(8) = 0.0
  y_init(9) = -trim_alt
  y_init(10:13) = euler_to_quat([phi, theta, psi])

  end subroutine set_initial_state_from_trim_shss

  subroutine set_initial_state_from_trim_shss_beta_spec(x, beta_specified)
  implicit none
  real, intent(in) :: x(6)
  real, intent(in) :: beta_specified
  real :: alpha, beta, u, v, w, phi, theta, psi

  alpha = x(1)
  phi = x(2)  ! phi is now an unknown
  beta = beta_specified  ! beta is specified
  controls(1) = x(3)  ! aileron
  controls(2) = x(4)  ! elevator
  controls(3) = x(5)  ! rudder
  controls(4) = x(6)  ! throttle

  ! Compute velocities
  u = trim_V * cos(alpha) * cos(beta)
  v = trim_V * sin(beta)
  w = trim_V * sin(alpha) * cos(beta)

  ! Euler angles
  theta = trim_theta
  psi = 0.0

  ! For SHSS: [p, q, r] = [0, 0, 0]

  ! Set initial state
  y_init(1) = u
  y_init(2) = v
  y_init(3) = w
  y_init(4) = 0.0  ! p = 0 for SHSS
  y_init(5) = 0.0  ! q = 0 for SHSS
  y_init(6) = 0.0  ! r = 0 for SHSS
  y_init(7) = 0.0
  y_init(8) = 0.0
  y_init(9) = -trim_alt
  y_init(10:13) = euler_to_quat([phi, theta, psi])

end subroutine set_initial_state_from_trim_shss_beta_spec

subroutine print_trim_results_shss(x)
  implicit none
  real, intent(in) :: x(6)
  real :: alpha_deg, beta_deg, da_deg, de_deg, dr_deg
  real :: phi_deg, theta_deg
  real :: y(13), V, dydt(13)

  alpha_deg = x(1) * 180.0 / PI
  beta_deg = x(2) * 180.0 / PI
  da_deg = x(3) * 180.0 / PI
  de_deg = x(4) * 180.0 / PI
  dr_deg = x(5) * 180.0 / PI

  phi_deg = trim_phi * 180.0 / PI
  theta_deg = trim_theta * 180.0 / PI

  V = sqrt(y_init(1)**2 + y_init(2)**2 + y_init(3)**2)

  write(*,'(A)') ''
  write(*,'(A)') '============================================'
  write(*,'(A)') 'Trim Solution for Steady-Heading Sideslip'
  write(*,'(A)') '============================================'
  write(*,'(A,F10.2,A)') 'Velocity: ', V, ' ft/s'
  write(*,'(A,F10.2,A)') 'Altitude: ', trim_alt, ' ft'
  write(*,'(A,F16.13,A)') 'Bank angle (phi): ', phi_deg, ' deg'
  write(*,'(A,F16.13,A)') 'Elevation angle (theta): ', theta_deg, ' deg'
  write(*,'(A)') '--------------------------------------------'
  write(*,'(A,F16.13,A)') 'alpha = ', alpha_deg, ' deg'
  write(*,'(A,F16.13,A)') 'beta  = ', beta_deg, ' deg'
  write(*,'(A,F16.13,A)') 'p     = ', 0.0, ' deg/s'
  write(*,'(A,F16.13,A)') 'q     = ', 0.0, ' deg/s'
  write(*,'(A,F16.13,A)') 'r     = ', 0.0, ' deg/s'
  write(*,'(A,F16.13,A)') 'aileron  = ', da_deg, ' deg'
  write(*,'(A,F18.13,A)') 'elevator = ', de_deg, ' deg'
  write(*,'(A,F18.13,A)') 'rudder   = ', dr_deg, ' deg'
  write(*,'(A,F17.14)') 'throttle = ', x(6)
  write(*,'(A)') '============================================'

  y = y_init
  dydt = diff_eq(0.0, y, 1)
  write(*,'(A)') 'Residual verification (should be near zero):'
  write(*,'(A,6ES14.6)') 'dydt(1:6) = ', dydt(1:6)
  write(*,'(A)') ''

  if (print_state_aero) then
    call print_state_and_aero()
  end if

end subroutine print_trim_results_shss

subroutine print_trim_results_shss_beta_spec(x, beta_specified)
  implicit none
  real, intent(in) :: x(6)
  real, intent(in) :: beta_specified
  real :: alpha_deg, beta_deg, da_deg, de_deg, dr_deg
  real :: phi_deg, theta_deg
  real :: y(13), V, dydt(13)

  alpha_deg = x(1) * 180.0 / PI
  beta_deg = beta_specified * 180.0 / PI
  phi_deg = x(2) * 180.0 / PI
  da_deg = x(3) * 180.0 / PI
  de_deg = x(4) * 180.0 / PI
  dr_deg = x(5) * 180.0 / PI

  theta_deg = trim_theta * 180.0 / PI

  V = sqrt(y_init(1)**2 + y_init(2)**2 + y_init(3)**2)

  write(*,'(A)') ''
  write(*,'(A)') '============================================'
  write(*,'(A)') 'Trim Solution for Steady-Heading Sideslip'
  write(*,'(A)') '============================================'
  write(*,'(A,F10.2,A)') 'Velocity: ', V, ' ft/s'
  write(*,'(A,F10.2,A)') 'Altitude: ', trim_alt, ' ft'
  write(*,'(A,F16.13,A)') 'Bank angle (phi): ', phi_deg, ' deg'
  write(*,'(A,F16.13,A)') 'Elevation angle (theta): ', theta_deg, ' deg'
  write(*,'(A)') '--------------------------------------------'
  write(*,'(A,F16.13,A)') 'alpha = ', alpha_deg, ' deg'
  write(*,'(A,F16.13,A)') 'beta  = ', beta_deg, ' deg'
  write(*,'(A,F16.13,A)') 'p     = ', 0.0, ' deg/s'
  write(*,'(A,F16.13,A)') 'q     = ', 0.0, ' deg/s'
  write(*,'(A,F16.13,A)') 'r     = ', 0.0, ' deg/s'
  write(*,'(A,F16.13,A)') 'aileron  = ', da_deg, ' deg'
  write(*,'(A,F18.13,A)') 'elevator = ', de_deg, ' deg'
  write(*,'(A,F16.13,A)') 'rudder   = ', dr_deg, ' deg'
  write(*,'(A,F17.14)') 'throttle = ', x(6)
  write(*,'(A)') '============================================'

  y = y_init
  dydt = diff_eq(0.0, y, 1)
  write(*,'(A)') 'Residual verification (should be near zero):'
  write(*,'(A,6ES14.6)') 'dydt(1:6) = ', dydt(1:6)
  write(*,'(A)') ''

  if (print_state_aero) then
    call print_state_and_aero()
  end if

end subroutine print_trim_results_shss_beta_spec

subroutine print_state_and_aero()
  implicit none
  real :: y(13), alpha, beta, V, Vlat
  real :: z_ft, h_ft, T_R, p_psf, rho_slug, a_fts, qdyn
  real :: pbar, qbar, rbar, CL1, CLift, CSide, CDrag, Croll, Cpitch, Cyaw
  real :: eul(3)

  y = y_init
  call pseudo_aero(y)

  V = sqrt(y(1)**2 + y(2)**2 + y(3)**2)
  Vlat = sqrt(max(y(1)**2 + y(3)**2, 1.0e-16))
  alpha = atan2(y(3), y(1))
  beta = atan2(y(2), Vlat)

  z_ft = max(-y(9), 0.0)
  call standard_atmospheric_English(z_ft, h_ft, T_R, p_psf, rho_slug, a_fts)

  if (V > 1.0e-12) then
    pbar = 0.5*y(4)*lat_ref / V
    qbar = 0.5*y(5)*long_ref / V
    rbar = 0.5*y(6)*lat_ref / V
  else
    pbar = 0.0; qbar = 0.0; rbar = 0.0
  end if

  CL1 = CL0 + CL_alpha*alpha
  CLift = CL1 + CL_qbar*qbar + CL_deltae*controls(2)
  CSide = CS_beta*beta + (CS_pbar + CS_alphapbar*alpha)*pbar + CS_rbar*rbar + &
          CS_deltaa*controls(1) + CS_deltar*controls(3)
  CDrag = CD_L0 + CD_L*CL1 + CD_L2*CL1**2 + CD_S2*CSide**2 + (CD_qbar + &
          CD_alpha_qbar*alpha)*qbar + (CD_deltae + CD_alpha_deltae*alpha)*controls(2) + &
          CD_deltae2*controls(2)**2
  Croll = Cl_beta*beta + Cl_pbar*pbar + (Cl_rbar + Cl_alpha_rbar*alpha)*rbar + &
          Cl_deltaa*controls(1) + Cl_deltar*controls(3)
  Cpitch = Cm0 + Cm_alpha*alpha + Cm_qbar*qbar + Cm_deltae*controls(2)
  Cyaw = Cn_beta*beta + (Cn_pbar + Cn_alphapbar*alpha)*pbar + Cn_rbar*rbar + &
          (Cn_deltaa + Cn_alpha_deltaa*alpha)*controls(1) + Cn_deltar*controls(3)

  qdyn = 0.5*rho_slug*V*V

  eul = quat_to_euler(y(10:13))

  write(*,'(A)') 'State and Aerodynamic Information:'
  write(*,'(A)') '--------------------------------------------'
  write(*,'(A,F12.4,A)') 'V       = ', V, ' ft/s'
  write(*,'(A,F12.6,A)') 'alpha   = ', alpha*180.0/PI, ' deg'
  write(*,'(A,F12.6,A)') 'beta    = ', beta*180.0/PI, ' deg'
  write(*,'(A,F12.6,A)') 'pbar    = ', pbar, ''
  write(*,'(A,F12.6,A)') 'qbar    = ', qbar, ''
  write(*,'(A,F12.6,A)') 'rbar    = ', rbar, ''
  write(*,'(A,F12.4,A)') 'rho     = ', rho_slug, ' slug/ft^3'
  write(*,'(A,F12.4,A)') 'qdyn    = ', qdyn, ' lbf/ft^2'
  write(*,'(A,F12.6)')   'CL      = ', CLift
  write(*,'(A,F12.6)')   'CD      = ', CDrag
  write(*,'(A,F12.6)')   'CS      = ', CSide
  write(*,'(A,F12.6)')   'Cl      = ', Croll
  write(*,'(A,F12.6)')   'Cm      = ', Cpitch
  write(*,'(A,F12.6)')   'Cn      = ', Cyaw
  write(*,'(A,6F12.4)')  'FM      = ', FM

end subroutine print_state_and_aero

subroutine vehicle_print_aero_table
  implicit none
  integer :: i , iunit
  real :: alpha, beta, states(13)
  real :: ca , cb , sa, sb
  real :: CL , CD , Cm
  real :: h_ft,T,P, rho,a_sound
  real :: alt, Vel, qbar

  alt = trim_alt
  call standard_atmospheric_English(alt,h_ft,T,P,rho,a_sound)

  Vel = trim_V
  qbar = 0.5*rho*Vel**2
  
  controls = 0.0
  states = 0.0

  open(newunit=iunit, file='aero_table.csv', status='REPLACE')
  write(iunit,*)'alpha[deg],CL,CD,Cm'
  do i =-180,180,1
    alpha = real(i)*PI/180.0
    beta = 0.0

    states(1) = Vel*cos(alpha)*cos(beta)
    states(2) = Vel*sin(beta)
    states(3) = Vel*sin(alpha)*cos(beta)
    states(9) = -alt 

    call pseudo_aero(states)
    ca = cos(alpha)
    sa = sin(alpha)
    cb = cos(beta)
    sb = sin(beta)
    
    CL = (-FM(3)*ca + FM(1)*sa) / (qbar*sref)
    CD = (-FM(1)*ca*cb - FM(2)*sb - FM(3)*sa*cb) / (qbar*sref)
    Cm = FM(5) / (qbar*sref*long_ref)

    write(iunit,*) alpha*180.0/PI,',',CL,',',CD,',',Cm
  end do

  close(iunit)
end subroutine vehicle_print_aero_table

subroutine init(filename)
  implicit none
  real :: alpha,beta,V,eul(3)
  character(100), intent(in) :: filename
  real :: h0, temp0, p0, a0
  character(:), allocatable :: init_type
  real :: eul_angles(3)
  real :: CG_shift(3)
  real :: tau_temp(4)  ! Time constants [sec]
  character(:), allocatable :: temp_type
  type(json_value), pointer :: j_connections, j_graphics
  type(json_value), pointer :: j_controls_in


  open(rk4_file,FILE='RK4.TXT', STATUS="REPLACE", ACTION="WRITE")
  call jsonx_load(filename, j_main)
  call jsonx_get(j_main, 'simulation.time_step[sec]', dt)
  call jsonx_get(j_main, 'simulation.total_time[sec]', tf)
  call jsonx_get(j_main, 'simulation.rk4_verbrose', rk4_verbose, .false.)

  ! Aerodynamic references & location
  call jsonx_get(j_main, 'vehicle.aerodynamics.reference.area[ft^2]', sref)
  call jsonx_get(j_main, 'vehicle.aerodynamics.reference.longitudinal_length[ft]', long_ref)
  call jsonx_get(j_main, 'vehicle.aerodynamics.reference.lateral_length[ft]', lat_ref)
  if (.not. allocated(aero_ref_location)) then
    allocate(aero_ref_location(3))
  end if
 
  block
    real :: CG_shift(3), temp_location(3)
    CG_shift = 0.0
    temp_location = 0.0
    call get_vec3_from_json(j_main, 'vehicle.CG_shift[ft]', CG_shift)
    
    ! Check if relative_location is specified (and non-zero)
    if (json_key_exists(j_main, 'vehicle.aerodynamics.reference.relative_location[ft]')) then
      call get_vec3_from_json(j_main, 'vehicle.aerodynamics.reference.relative_location[ft]', temp_location)
      if (abs(temp_location(1)) > 1e-10 .or. abs(temp_location(2)) > 1e-10 .or. abs(temp_location(3)) > 1e-10) then
        aero_ref_location = temp_location
      else
        ! Use negative of CG shift
        aero_ref_location = -CG_shift
      end if
    else
      ! Use negative of CG shift
      aero_ref_location = -CG_shift
    end if
  end block

  ! Coefficients
  call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CL.0', CL0)
  call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CL.alpha', CL_alpha)
  call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CL.alphahat', CL_alphacap)
  call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CL.qbar', CL_qbar)
  call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CL.elevator', CL_deltae)
  call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CS.beta', CS_beta)
  call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CS.pbar', CS_pbar)
  call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CS.alpha_pbar',CS_alphapbar)
  call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CS.rbar', CS_rbar)
  call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CS.aileron', CS_deltaa)
  call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CS.rudder', CS_deltar)
  call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CD.L0', CD_L0)
  call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CD.CL1', CD_L)
  call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CD.CL1_CL1', CD_L2)
  call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CD.CS_CS', CD_S2)
  call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CD.qbar', CD_qbar)
  call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CD.alpha_qbar', CD_alpha_qbar)
  call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CD.elevator', CD_deltae)
  call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CD.alpha_elevator', CD_alpha_deltae)
  call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.CD.elevator_elevator',CD_deltae2)
  call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.Cl.beta', Cl_beta)
  call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.Cl.pbar', Cl_pbar)
  call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.Cl.rbar', Cl_rbar)
  call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.Cl.alpha_rbar',Cl_alpha_rbar)
  call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.Cl.aileron', Cl_deltaa)
  call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.Cl.rudder', Cl_deltar)
  call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.Cm.0', Cm0)
  call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.Cm.alpha', Cm_alpha)
  call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.Cm.qbar', Cm_qbar)
  call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.Cm.alphahat', Cm_alphacap)
  call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.Cm.elevator', Cm_deltae)
  call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.Cn.beta', Cn_beta)
  call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.Cn.pbar', Cn_pbar)
  call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.Cn.alpha_pbar', Cn_alphapbar)
  call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.Cn.rbar', Cn_rbar)
  call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.Cn.aileron', Cn_deltaa)
  call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.Cn.alpha_aileron', Cn_alpha_deltaa)
  call jsonx_get(j_main, 'vehicle.aerodynamics.coefficients.Cn.rudder', Cn_deltar)
  call jsonx_get(j_main, 'vehicle.thrust.T0[lbf]', T0)
  call jsonx_get(j_main, 'vehicle.thrust.Ta', Ta)
  
  call jsonx_get(j_main, 'vehicle.aerodynamics.stall.include_stall', include_stall)
  call jsonx_get(j_main, 'vehicle.aerodynamics.stall.CL.alpha_0[deg]', CLstall%alpha_0)
  CLstall%alpha_0 = CLstall%alpha_0*PI/180.0
  call jsonx_get(j_main, 'vehicle.aerodynamics.stall.CL.alpha_s[deg]', CLstall%alpha_s)
  CLstall%alpha_s = CLstall%alpha_s*PI/180.0
  call jsonx_get(j_main, 'vehicle.aerodynamics.stall.CL.lambda_b', CLstall%lambda_b)
  CLstall%minval = 0.0  
  
  call jsonx_get(j_main, 'vehicle.aerodynamics.stall.CD.alpha_0[deg]', CDstall%alpha_0)
  CDstall%alpha_0 = CDstall%alpha_0*PI/180.0
  call jsonx_get(j_main, 'vehicle.aerodynamics.stall.CD.alpha_s[deg]', CDstall%alpha_s)
  CDstall%alpha_s = CDstall%alpha_s*PI/180.0
  call jsonx_get(j_main, 'vehicle.aerodynamics.stall.CD.lambda_b', CDstall%lambda_b)
  CDstall%minval = 0.0  
  
  call jsonx_get(j_main, 'vehicle.aerodynamics.stall.Cm.min', Cmstall%minval)
  call jsonx_get(j_main, 'vehicle.aerodynamics.stall.Cm.alpha_0[deg]', Cmstall%alpha_0)
  Cmstall%alpha_0 = Cmstall%alpha_0*PI/180.0
  call jsonx_get(j_main, 'vehicle.aerodynamics.stall.Cm.alpha_s[deg]', Cmstall%alpha_s)
  Cmstall%alpha_s = Cmstall%alpha_s*PI/180.0
  call jsonx_get(j_main, 'vehicle.aerodynamics.stall.Cm.lambda_b', Cmstall%lambda_b)

  call jsonx_get(j_main, 'vehicle.aerodynamics.stall.Cm.lambda_b', Cmstall%lambda_b)
  call jsonx_get(j_main, 'vehicle.actuators.enable', use_actuator_dynamics, .false.)
  
  if (use_actuator_dynamics) then
    
    write(*,'(A)') ''
    write(*,'(A)') '================================================'
    write(*,'(A)') 'Actuator Dynamics ENABLED (First-Order Model)'
    write(*,'(A)') '================================================'
    
    ! Load time constants τ for each actuator using Eq. 10.2.1
    call jsonx_get(j_main, 'vehicle.actuators.aileron_tau[s]', tau_temp(1), 0.05)
    call jsonx_get(j_main, 'vehicle.actuators.elevator_tau[s]', tau_temp(2), 0.05)
    call jsonx_get(j_main, 'vehicle.actuators.rudder_tau[s]', tau_temp(3), 0.05)
    call jsonx_get(j_main, 'vehicle.actuators.throttle_tau[s]', tau_temp(4), 1.0)
    
    ! Convert time constants to damping rates: σ = 1/τ
    sigma_actuator(1) = 1.0 / tau_temp(1)
    sigma_actuator(2) = 1.0 / tau_temp(2)
    sigma_actuator(3) = 1.0 / tau_temp(3)
    sigma_actuator(4) = 1.0 / tau_temp(4)
    
    ! Load rate limits [deg/s for control surfaces, 1/s for throttle]
    call jsonx_get(j_main, 'vehicle.actuators.aileron_rate_limit[deg/s]', &
                   actuator_rate_limits(1), 80.0)
    call jsonx_get(j_main, 'vehicle.actuators.elevator_rate_limit[deg/s]', &
                   actuator_rate_limits(2), 60.0)
    call jsonx_get(j_main, 'vehicle.actuators.rudder_rate_limit[deg/s]', &
                   actuator_rate_limits(3), 120.0)
    call jsonx_get(j_main, 'vehicle.actuators.throttle_rate_limit[1/s]', &
                   actuator_rate_limits(4), 1.0)
    
    actuator_rate_limits(1:3) = actuator_rate_limits(1:3) * PI / 180.0
    
    ! Load position limits
    call jsonx_get(j_main, 'vehicle.actuators.aileron_limits[deg]', actuator_pos_limits(1,1), -21.5)
    call jsonx_get(j_main, 'vehicle.actuators.aileron_limits[deg]', actuator_pos_limits(1,2), 21.5)
    call jsonx_get(j_main, 'vehicle.actuators.elevator_limits[deg]', actuator_pos_limits(2,1), -25.0)
    call jsonx_get(j_main, 'vehicle.actuators.elevator_limits[deg]', actuator_pos_limits(2,2), 25.0)
    call jsonx_get(j_main, 'vehicle.actuators.rudder_limits[deg]', actuator_pos_limits(3,1), -30.0)
    call jsonx_get(j_main, 'vehicle.actuators.rudder_limits[deg]', actuator_pos_limits(3,2), 30.0)
    call jsonx_get(j_main, 'vehicle.actuators.throttle_limits', actuator_pos_limits(4,1), 0.0)
    call jsonx_get(j_main, 'vehicle.actuators.throttle_limits', actuator_pos_limits(4,2), 1.0)
    
    actuator_pos_limits(1:3,:) = actuator_pos_limits(1:3,:) * PI / 180.0
    
    write(*,'(A)') 'Actuator Parameters:'
    write(*,'(A,4F8.4,A)') '  Time constants τ: ', tau_temp, ' [sec]'
    write(*,'(A,4F8.2,A)') '  Damping rates σ:  ', sigma_actuator, ' [1/sec]'
    write(*,'(A,3F8.2,F8.4,A)') '  Rate limits:      ', &
           actuator_rate_limits(1:3)*180.0/PI, actuator_rate_limits(4), ' [deg/s, 1/s]'
    write(*,'(A)') '================================================'
    write(*,'(A)') ''
  else
    write(*,'(A)') 'Actuator Dynamics DISABLED (instant response)'
  end if


  call jsonx_get(j_main, 'vehicle.initial.type', init_type)

  call mass_inertia()

  call standard_atmospheric_English(0.0, h0, temp0, p0, rho0, a0)

  call jsonx_get(j_main, 'analysis.print_state_aero', print_state_aero, .false.)


  if (trim(adjustl(init_type)) == 'trim') then
    ! Read trim parameters before solving
    eul_angles = 0.0
    
    ! Check if climb angle is specified instead of elevation angle
    if (json_key_exists(j_main, 'vehicle.initial.trim.climb_angle[deg]')) then
      ! Use climb angle - only read phi from Euler_angles
      call get_vec3_from_json(j_main, 'vehicle.initial.Euler_angles[deg]', eul_angles)
      trim_phi = eul_angles(1) * PI / 180.0
      trim_theta = 0.0  
      write(*,'(A)') 'Climb angle specified - elevation angle will be computed during trim'
    else
      ! Use elevation angle 
      call get_vec3_from_json(j_main, 'vehicle.initial.Euler_angles[deg]', eul_angles)
      trim_phi   = eul_angles(1) * PI / 180.0
      trim_theta = eul_angles(2) * PI / 180.0
      write(*,'(A)') 'Elevation angle specified directly'
    end if
    
    call jsonx_get(j_main, 'vehicle.initial.altitude[ft]', trim_alt)
    call jsonx_get(j_main, 'vehicle.initial.airspeed[ft/s]', trim_V)
    
    call jsonx_get(j_main, 'vehicle.initial.trim.type', temp_type, 'sct')
    trim_type = temp_type
    
    if (trim(adjustl(trim_type)) == 'sct') then
      call solve_trim_sct()
    else if (trim(adjustl(trim_type)) == 'shss') then
      call solve_trim_shss()
    else
      write(*,'(A,A)') 'Error: Unknown trim type: ', trim(trim_type)
      stop
    end if
  else
  
    ! INIT_TYPE == STATE 
    y_init = 0.0

    ! Basic initial conditions
    call jsonx_get(j_main, 'vehicle.initial.airspeed[ft/s]', V)
    call jsonx_get(j_main, 'vehicle.initial.altitude[ft]', trim_alt)
    trim_V = V  ! Store for aero table generation

    ! Euler angles [deg] → [rad]  (phi, theta, psi)
    call get_vec3_from_json(j_main, 'vehicle.initial.Euler_angles[deg]', eul_angles)
    eul = eul_angles * PI / 180.0

    ! Alpha / beta from the nested "state" block 
    call jsonx_get(j_main, 'vehicle.initial.state.angle_of_attack[deg]', alpha, 0.0)
    call jsonx_get(j_main, 'vehicle.initial.state.sideslip_angle[deg]',  beta,  0.0)
    alpha = alpha * PI / 180.0
    beta  = beta  * PI / 180.0

    ! Body velocity components from V, alpha, beta
    y_init(1) = V * cos(alpha) * cos(beta)
    y_init(2) = V * sin(beta)
    y_init(3) = V * sin(alpha) * cos(beta)

    call jsonx_get(j_main, 'vehicle.initial.state.p[deg/s]', y_init(4), 0.0)
    call jsonx_get(j_main, 'vehicle.initial.state.q[deg/s]', y_init(5), 0.0)
    call jsonx_get(j_main, 'vehicle.initial.state.r[deg/s]', y_init(6), 0.0)
    y_init(4:6) = y_init(4:6) * PI / 180.0

    call jsonx_get(j_main, 'vehicle.initial.altitude[ft]', y_init(9))
    y_init(9) = -y_init(9)

    y_init(10:13) = euler_to_quat(eul)

    call jsonx_get(j_main, 'vehicle.initial.state.aileron[deg]',  controls(1), 0.0)
    call jsonx_get(j_main, 'vehicle.initial.state.elevator[deg]', controls(2), 0.0)
    call jsonx_get(j_main, 'vehicle.initial.state.rudder[deg]',   controls(3), 0.0)
    call jsonx_get(j_main, 'vehicle.initial.state.throttle',      controls(4), 0.0)

    controls(1:3) = controls(1:3) * PI / 180.0
  end if

  actuator_states = controls

  ! Connections
  call jsonx_get(j_main, 'connections',j_connections)
  call jsonx_get(j_connections, 'graphics',j_graphics)
  call graphics%init(j_graphics)

  write(*,*)'Initialising vehicle..'
  
  write(*,*)'Generating aerodynamic table..'
  call vehicle_print_aero_table()

  call jsonx_get(j_connections, 'controls_in', j_controls_in)
  call controls_in%init(j_controls_in)


end subroutine init

subroutine run()
  implicit none
  real :: t, y(13), y1(13), eul(3), alpha, beta, s(27)
  real :: cpu_start_time, cpu_end_time, time1, time2,actual_time, integrated_time
  logical :: real_time
  real :: received(4)

  real_time = .false.
  
  if (abs(dt) < TOLERANCE) then
    real_time = .true.
    call cpu_time(time1)
    y1 = rk4(0.0, y_init, 0.0)
    call quat_norm(y1(10:13))
    call cpu_time(time2)
    dt = time2 - time1
  end if

  t = 0.0
  y = y_init

  write(*,*) 't[s], u, v, w, p, q, r, x, y, z, e0, ex, ey, ez, alpha[rad], beta[rad]'
  write(*,'(14ES20.12)') t, y(:)

  call cpu_time(cpu_start_time)
  time1 = cpu_start_time
  integrated_time = 0.0

  do while (t < tf)
  
    received = controls_in%recv()
    controls(1) = received(1) * PI / 180.0 ! aileron [rad]
    controls(2) = received(2) * PI / 180.0 ! elevator [rad]
    controls(3) = received(3) * PI / 180.0! rudder [rad]
    controls(4) = received(4) ! throttle
    
    call update_actuators(controls, dt)
    
    y1 = rk4(t, y, dt)
    call quat_norm(y1(10:13))
    y = y1
    t = t + dt
    integrated_time = integrated_time + dt

    if (rk4_verbose) then
      
      write(*,*)
      write(*,*) 't[s], u, v, w, p, q, r, x, y, z, e0, ex, ey, ez, alpha[rad], beta[rad]'
    end if
    write(*,'(15ES20.12)') t, dt, y(:)

    ! send graphics over connection
    ! Compute alpha and beta from current state
    call compute_alpha_beta(y, alpha, beta)
      
    ! Compute Euler angles from quaternion 
    eul = quat_to_euler(y(10:13))
    
    s(1)      = t           ! time
    s(2:14)   = y(1:13)     ! u,v,w,p,q,r,x,y,z,e0,ex,ey,ez  
    s(15)     = alpha       ! angle of attack
    s(16)     = beta        ! sideslip angle
    s(17:19)  = eul(1:3)    ! phi, theta, psi (Euler angles)
    s(20)     = controls(4) ! throttle commanded
    s(21)     = dt          ! time step
    ! Actuator commanded positions (controls = u(t))
    s(22)     = controls(1) ! aileron commanded [rad]
    s(23)     = controls(2) ! elevator commanded [rad]
    s(24)     = controls(3) ! rudder commanded [rad]
    ! Actuator actual positions (actuator_states = δ(t))
    s(25)     = actuator_states(1) ! aileron actual [rad]
    s(26)     = actuator_states(2) ! elevator actual [rad]
    s(27)     = actuator_states(3) ! rudder actual [rad]
   
    call graphics%send(s)
    
    if (real_time) then
      call cpu_time(time2)
      dt = time2 - time1
      time1 = time2
    end if
  end do

  call cpu_time(cpu_end_time)

  actual_time = cpu_end_time - cpu_start_time
  write(*,*)'Total integrated time[s] = ',integrated_time
  write(*,*)'Total actual elapse time[s] = ',actual_time
  write(*,*)'Total error in time[s] = ',integrated_time - actual_time
  
end subroutine run


! Compute aerodynamic angles from state
subroutine compute_alpha_beta(y, alpha, beta)
  implicit none
  real, intent(in) :: y(13)
  real, intent(out):: alpha, beta
  real :: V, Vlat
  V = sqrt(max(y(1)**2 + y(2)**2 + y(3)**2, 1.0e-16))
  Vlat = sqrt(max(y(1)**2 + y(3)**2, 1.0e-16))
  alpha = atan2(y(3), y(1))
  beta = atan2(y(2), Vlat)
end subroutine compute_alpha_beta

end module sim_m