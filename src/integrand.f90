module integrand

  use nrtype
  use comvar

contains

  subroutine setkinematics
  implicit none
  real(sp) :: p1x,p1y,p2x,p2y
  real(sp) :: qtx,qty,ptx,pty
  real(sp) :: k1x,k1y,k2x,k2y
  real(sp) :: kax,kay,kbx,kby

  ! set lepton qT and PT
  thp1 = 0d0  ! fix lepton 1 azimuth angle at 0
  p1x = p1t*cos(thp1);  p1y = p1t*sin(thp1)
  p2x = p2t*cos(thp2);  p2y = p2t*sin(thp2)
  qtx = p1x + p2x;      qty = p1y + p2y
  ptx = p1x - p2x;      pty = p1y - p2y
  qt = sqrt(qtx*qtx + qty*qty)
  thqt = atan2(qty,qtx)
  pt = sqrt(ptx*ptx + pty*pty) / 2d0
  thpt = atan2(pty,ptx)

  ! set lepton transverse mass
  m1t = sqrt(m1*m1 + p1t*p1t)
  m2t = sqrt(m2*m2 + p2t*p2t)
  
  ! set momentum fraction
  xa = (m1t*exp(+y1) + m2t*exp(+y2))/CME
  xb = (m1t*exp(-y1) + m2t*exp(-y2))/CME

  ! set Mandelstam variable
  mans = +xa*xb*CME2
  mant = m1*m1 - xa*CME*m1t*exp(-y1)
  manu = m2*m2 - xa*CME*m2t*exp(-y2)
  Mll = sqrt(mans)

  ! set photon momentum
  k1x = k1t*cos(thk1);  k1y = k1t*sin(thk1)
  kax = kat*cos(thka);  kay = kat*sin(thka)
  k2x = qtx - k1x;      k2y = qty - k1y
  kbx = qtx - kax;      kby = qty - kay
  k2t = sqrt(k2x*k2x + k2y*k2y)
  thk2 = atan2(k2y,k2x)
  kbt = sqrt(kbx*kbx + kby*kby)
  thkb = atan2(kby,kbx)
  k1sq = xa*xa*M_pro*M_pro + k1t*k1t
  kasq = xa*xa*M_pro*M_pro + kat*kat
  k2sq = xb*xb*M_pro*M_pro + k2t*k2t
  kbsq = xb*xb*M_pro*M_pro + kbt*kbt

  ! photon momentum difference (should be the same)
  k1ka = sqrt((k1x-kax)**2 + (k1y-kay)**2)
  k2kb = sqrt((k2x-kbx)**2 + (k2y-kby)**2)

  ! acoplanarity
  aco = abs(thp1-thp2)
  aco = merge(twoPI-aco, aco, aco.gt.PI)
  aco = 1d0 - aco/PI

  ! asymmetry
  asym = abs(p1t-p2t)/(p1t+p2t)

  ! kperp
  kt = PI * aco * pt
  
  ! set dot-cross term
  dotcross = (k1x*kax + k1y*kay) * (k2x*kbx + k2y*kby) &
            -(k1x*k2x + k1y*k2y) * (kax*kbx + kay*kby) &
            +(k1x*kbx + k1y*kby) * (k2x*kax + k2y*kay) 
  end subroutine setkinematics

  function diff_xsec() result(dxs)
  use photonff
  real(sp) :: dxs
  real(sp) :: bessel,sigma0
  bessel = bp * bessel_j0(k1ka*bp) !/ twoPI
  sigma0 = 2d0 * alphae**2 / mans / mans * (mant/manu + manu/mant)
  dxs = p1t * p2t * k1t * kat
  dxs = dxs * xf(sqrt(k1sq) , sqrt(kasq))
  dxs = dxs * xf(sqrt(k2sq) , sqrt(kbsq))
  dxs = dxs * bessel * sigma0 * dotcross
  dxs = dxs * (gevfm*1e2)**2    ! GeV^-2 -> fm^2 -> micro-barn
  return
  end function diff_xsec

end module integrand
