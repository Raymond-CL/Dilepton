module initial
  
  use comvar

contains

  subroutine initialize
  CME2 = CME*CME
  nprn = -1
  avgi = 0d0; sd = 0d0; chi2a = 0d0
  if(Mli.eq.1) then       ! set lepton masses
    m1 = M_ele; m2 = M_ele
  elseif(Mli.eq.2) then
    m1 = M_mu;  m2 = M_mu
  elseif(Mli.eq.3) then
    m1 = M_tau; m2 = M_tau
  else
    call exit(13)   !m1 = 0d0;   m2 = 0d0
  endif
  a0 = 1d0/sqrt(2d0) /gevfm
  RA = 1.13d0 * atomA**(1d0/3d0) /gevfm
  rho0 = 3d0*atomA/4d0/PI/(RA**3)
  bpmin = 15d0/gevfm;    bpmax = 100d0/gevfm
  thpmin = (1d0-acocut)*PI; thpmax = (1d0+acocut)*PI
  end subroutine initialize 
  
  subroutine set_limits
  ndim = 10; 
  region(1) = ptmin;      region(11) = ptmax        ! dp1t
  region(2) = ptmin;      region(12) = ptmax        ! dp2t
  region(3) = thpmin;     region(13) = thpmax       ! dthp2
  region(4) = ylmin;      region(14) = ylmax        ! dy1
  region(5) = ylmin;      region(15) = ylmax        ! dy2
  region(6) = ktmin;      region(16) = ktmax        ! dk1t
  region(7) = 0d0;        region(17) = twoPI        ! dthk1
  region(8) = ktmin;      region(18) = ktmax        ! dkat
  region(9) = 0d0;        region(19) = twoPI        ! dthka
  region(10) = bpmin;     region(20) = bpmax        ! dbp
  end subroutine set_limits

end module initial
