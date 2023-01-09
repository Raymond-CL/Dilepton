module photonff
  
  use nrtype
  use comvar

contains

  ! photon distribution
  function xf(k,kp)
  implicit none
  real(sp) :: xf
  real(sp), intent(in) :: k,kp
  xf = atomZ*atomZ*alphae/PI/PI
  xf = xf * ffac(k)/k/k * ffac(kp)/kp/kp
  return
  end function xf

  ! photon form factor
  function ffac(k)
  implicit none
  real(sp) :: ffac
  real(sp), intent(in) :: k
  real(sp), parameter :: Qs2 = 0.08d0
  if(formfac) then
    ffac = 4d0*PI*rho0/atomA/k**3
    ffac = ffac / (1d0+a0**2*k**2)
    ffac = ffac * (sin(k*RA) - k*RA * cos(k*RA))
  else
    ffac = exp(-k*k/Qs2)
  endif
  return
  end function ffac

end module photonff
