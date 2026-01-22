module babu_m
  implicit none
  real, parameter :: PI = 3.14159265358979323846
  real, parameter :: TOLERANCE = 1.0e-12
  real, parameter :: alt_list(21) = (/ &
       0.0,  5000.0, 10000.0, 15000.0, 20000.0, 25000.0, 30000.0, 35000.0, &
   40000.0, 45000.0, 50000.0, 55000.0, 60000.0, 65000.0, 70000.0, 75000.0, &
   80000.0, 85000.0, 90000.0, 95000.0,100000.0 /)

contains
!---------------------------------------------------------------------
function quat_mult(q1,q2) result(ans)
  implicit none
  real, intent(in) :: q1(4), q2(4)
  real :: ans(4)
  ans(1) = q1(1)*q2(1) - q1(2)*q2(2) - q1(3)*q2(3) - q1(4)*q2(4)
  ans(2) = q1(1)*q2(2) - q1(2)*q2(1) - q1(3)*q2(4) - q1(4)*q2(3)
  ans(3) = q1(1)*q2(3) - q1(2)*q2(4) - q1(3)*q2(1) - q1(4)*q2(2)
  ans(4) = q1(1)*q2(4) - q1(2)*q2(3) - q1(3)*q2(2) - q1(4)*q2(1)
end function quat_mult

function quat_base_to_dependent(vec, quat) result(ans)
  implicit none
  real, intent(in) :: vec(3), quat(4)
  real :: ans(3), qv(3), cross1(3), cross2(3)
  real :: q0
  
  q0 = quat(1)
  qv(1) = quat(2)
  qv(2) = quat(3)
  qv(3) = quat(4)
  
  ! cross1 = qv x v
  cross1(1) = qv(2)*vec(3) - qv(3)*vec(2)
  cross1(2) = qv(3)*vec(1) - qv(1)*vec(3)
  cross1(3) = qv(1)*vec(2) - qv(2)*vec(1)
  
  ! cross2 = qv x cross1
  cross2(1) = qv(2)*cross1(3) - qv(3)*cross1(2)
  cross2(2) = qv(3)*cross1(1) - qv(1)*cross1(3)
  cross2(3) = qv(1)*cross1(2) - qv(2)*cross1(1)
  
  ans = vec - 2.0*q0*cross1 + 2.0*cross2
end function quat_base_to_dependent

function quat_dependent_to_base (vec, quat) result(ans)
  implicit none
  real, intent(in) :: vec(3), quat(4)
  real :: temp(4), ans(3)
  temp(1) =  quat(1)
  temp(2) = -quat(2)
  temp(3) = -quat(3)
  temp(4) = -quat(4)
  ans = quat_base_to_dependent(vec,temp)
end function quat_dependent_to_base

subroutine quat_norm(quat)
  implicit none
  real, intent(inout) :: quat(4)
  real :: denom
  denom = sqrt(quat(1)**2 + quat(2)**2 + quat(3)**2 + quat(4)**2)
  quat = quat/denom
end subroutine quat_norm

function euler_to_quat(eul) result(quat)
  implicit none
  real, intent(in) :: eul(3)
  real :: quat(4)
  real :: c1,c2,c3,s1,s2,s3
  c1 = cos(eul(1)*0.5);  s1 = sin(eul(1)*0.5)
  c2 = cos(eul(2)*0.5);  s2 = sin(eul(2)*0.5)
  c3 = cos(eul(3)*0.5);  s3 = sin(eul(3)*0.5)
  quat(1) =  c1*c2*c3 + s1*s2*s3
  quat(2) =  s1*c2*c3 - c1*s2*s3
  quat(3) =  c1*s2*c3 + s1*c2*s3
  quat(4) =  c1*c2*s3 - s1*s2*c3
end function euler_to_quat

function quat_to_euler(quat) result(eul)
  implicit none
  real, intent(in) :: quat(4)
  real :: eul(3)
  real :: condition, e0, ex, ey, ez, e02, ex2,ey2,ez2
  real :: s2, arg
  e0 = quat(1); ex = quat(2); ey = quat(3); ez = quat(4)
  s2 = 1.4142135623730950
  condition = e0*ey - ex*ez
  if ( abs( abs(condition) - 0.5) < TOLERANCE) then
    arg = ex * s2
    if (arg > 1.0) arg = 1.0
    if (arg < -1.0) arg = -1.0
    if (condition > 0.0) then
      eul(1) = 2.0*asin(arg)
      eul(2) = 1.5707963269489
      eul(3) = 0.0
    else
      eul(1) = 2.0*asin(arg)
      eul(2) = -1.5707963269489
      eul(3) = 0.0
    end if
  else
    e02 = e0**2; ex2 = ex**2; ey2 = ey**2; ez2 = ez**2
    eul(1) = atan2( 2.0*(e0*ex + ey*ez), (e02 + ez2 - ey2 - ex2) )
    arg    = 2.0*condition
    if (arg > 1.0)  arg = 1.0
    if (arg < -1.0) arg = -1.0
    eul(2) = asin(arg)
    eul(3) = atan2( 2.0*(e0*ez + ey*ex), (e02 + ex2 - ey2 - ez2) )
  end if
end function quat_to_euler

function gravity(h) result(g)
  implicit none
  real, intent(in) :: h
  real :: g
  real, parameter :: g0 = 9.80665
  real, parameter :: re = 6356766.0
  g = g0 * (re / (re + h))**2
end function gravity

function gravity_f(h) result(g)
  implicit none
  real, intent(in) :: h
  real :: g
  real, parameter :: g0 = 32.17404855643044
  real, parameter :: re = 20855531.5
  g = g0 * (re / (re + h))**2
end function gravity_f

subroutine test_gravity_SI_table(io_unit)
  implicit none
  integer, intent(in) :: io_unit
  integer :: k
  real :: g
  write(io_unit,'(A)') 'Gravity vs Altitude (SI)'
  write(io_unit,'(A)') 'z_geometric [m]   g [m/s^2]'
  do k = 1, size(alt_list)
    g = gravity(alt_list(k))
    write(io_unit,'(F12.1,3X,F12.6)') alt_list(k), g
  end do
end subroutine test_gravity_SI_table

subroutine test_gravity_English_table(io_unit)
  implicit none
  integer, intent(in) :: io_unit
  integer :: k
  real :: g
  write(io_unit,'(A)') 'Gravity vs Altitude (English units)'
  write(io_unit,'(A)') 'z_geometric [ft]  g [ft/s^2]'
  do k = 1, size(alt_list)
    g = gravity_f(alt_list(k))
    write(io_unit,'(F12.1,3X,F12.6)') alt_list(k), g
  end do
end subroutine test_gravity_English_table

subroutine standard_atmosphere_SI(z, h, T, p, rho, a)
  implicit none
  real, intent(in)  :: z
  real, intent(out) :: h, T, p, rho, a
  real, parameter :: re    = 6356766.0
  real, parameter :: g0    = 9.80665
  real, parameter :: R     = 287.0528
  real, parameter :: gamma = 1.4

  ! Geopotential layers & tables (8 layers)
  real, parameter :: height_tab(8) = (/ &
      0.0, 11000.0, 20000.0, 32000.0, 47000.0, 52000.0, 61000.0, 79000.0 /)
  real, parameter :: temp_tab(8) = (/ &
      288.150, 216.650, 216.650, 228.650, 270.650, 270.650, 252.650, 180.650 /)
  real, parameter :: lapse_tab(8) = (/ &
     -0.0065, 0.0, 0.0010, 0.0028, 0.0, -0.0020, -0.0040, 0.0 /)
  real, parameter :: pressure_tab(8) = (/ &
      101325.0, 22632.018222212, 5474.87352827083, 868.014769086723, &
      110.905588989225, 59.0005242789244, 18.2099249050177, 1.03770045489203 /)

  integer :: i
  ! write(*,*)'z =', z

  ! geometric->geopotential
  h = re * z / (re + z)

  ! pick layer

  i = 1
  do while (i < size(height_tab) .and. h >= height_tab(i+1))
    i = i + 1
  end do

  ! Temperature
  T = temp_tab(i) + lapse_tab(i) * (h - height_tab(i))
  ! write(*,*)'TEMP TAB=', temp_tab(i)
  ! write(*,*)'Lapse TAB=', lapse_tab(i)
  ! write(*,*)'Height TAB=', height_tab(i)
  ! write(*,*)'h =', h


  ! Pressure (gradient vs isothermal)
  if (abs(lapse_tab(i)) > 1.0e-12) then
    p = pressure_tab(i) * ( T / temp_tab(i) ) ** ( -g0 / (R * lapse_tab(i)) )
  else
    p = pressure_tab(i) * exp( -g0 * (h - height_tab(i)) / (R * temp_tab(i)) )
  end if

  ! Density & speed of sound
  rho = p / (R * T)
  ! write(*,*)'Pre SI=', p
  ! write(*,*)'R SI=', R
  ! write(*,*)'Temp SI=', T
  ! write(*,*)'DEnsity SI=', rho
  a   = sqrt(gamma * R * T)
end subroutine standard_atmosphere_SI

subroutine standard_atmospheric_SI_table(io_unit)
  implicit none
  integer, intent(in) :: io_unit
  integer :: z_m
  real :: z, h, T, p, rho, a
  write(io_unit,'(A)') 'Properties vs. Altitude (SI)'
  write(io_unit,'(A)') 'Columns: z[m] h[m] T[K] p[Pa] rho[kg/m^3] a[m/s]'
  write(io_unit,*)
  do z_m = 0, 90000, 5000
    z = real(z_m)
    call standard_atmosphere_SI(z, h, T, p, rho, a)
    write(io_unit,'(F10.1,1X,F10.1,1X,F9.2,1X,ES12.5,1X,ES12.5,1X,F9.3)') &
         z, h, T, p, rho, a
  end do
end subroutine standard_atmospheric_SI_table

subroutine standard_atmospheric_English(z_ft, h_ft, T_R, p_psf, rho_slug, a_fts)
  implicit none
  real, intent(in)  :: z_ft
  real, intent(out) :: h_ft, T_R, p_psf, rho_slug, a_fts
  real :: z_m, h_m, T_K, p_Pa, rho_si, a_si
  real, parameter :: ft_per_m = 3.28084
  real, parameter :: R_per_K  = 9.0/5.0
  real, parameter :: psf_per_Pa = 0.0208854342
  real, parameter :: slugft3_per_kgm3 = 0.00194032
  real, parameter :: fts_per_ms = 3.28084

  !z_m = z_ft / ft_per_m
  z_m = z_ft * 0.3048
  call standard_atmosphere_SI(z_m, h_m, T_K, p_Pa, rho_si, a_si)
  h_ft    = h_m / 0.3048
  T_R     = T_K * R_per_K
  p_psf   = p_Pa * psf_per_Pa
  !rho_slug= rho_si * slugft3_per_kgm3 
  rho_slug = rho_si/515.379
  a_fts   = a_si * fts_per_ms
end subroutine standard_atmospheric_English

subroutine standard_atmospheric_English_table(io_unit)
  implicit none
  integer, intent(in) :: io_unit
  integer :: zft
  real :: z_ft, h_ft, T_R, p_psf, rho_slug, a_fts
  write(io_unit,'(A)') 'Standard Atmosphere (English units)'
  write(io_unit,'(A)') 'Columns: z_geometric[ft] h_geopotential[ft] T[deg R] p[lbf/ft^2] rho[slug/ft^3] a[ft/s]'
  write(io_unit,*)
  do zft = 0, 200000, 10000
    z_ft = real(zft)
    call standard_atmospheric_English(z_ft, h_ft, T_R, p_psf, rho_slug, a_fts)
    write(io_unit,'(F10.1,1X,F12.1,1X,F9.2,1X,ES12.5,1X,ES12.5,1X,F9.2)') &
         z_ft, h_ft, T_R, p_psf, rho_slug, a_fts
  end do
end subroutine standard_atmospheric_English_table

end module babu_m