module comvar

  use nrtype
  implicit none
  public
 
  ! conversion units
  real(sp), parameter :: planck = 6.62607015d-34    ! Planck constant (J.s)
  real(sp), parameter :: planckbar = planck/twoPI   ! reduced Planck (J.s)
  real(sp), parameter :: charge = 1.602176634d-19   ! chrg mag (J/eV)
  real(sp), parameter :: hbar = planckbar/charge    ! hbar (eV.s)
  real(sp), parameter :: light = 299792458d0        ! light speed (m/s)
  real(sp), parameter :: hbarc = hbar*light         ! conversion unit (eV.m)
  real(sp), parameter :: giga = 1d+9
  real(sp), parameter :: femto = 1d-15
  real(sp), parameter :: gevfm = hbarc/giga/femto   ! common conversion (GeV.fm)

  ! temporary test variables
  real    :: tempr1,tempr2
  integer :: tempi1,tempi2
   
  ! center of mass energy
  real(sp) :: CME,CME2

  ! lepton transverse momentum, azi_angle, rapidity, mass
  real(sp) :: ptmin,ptmax
  real(sp) :: p1t,p2t
  real(sp) :: thpmin,thpmax
  real(sp) :: thp1,thp2
  real(sp) :: ylmin,ylmax
  real(sp) :: y1,y2
  real(sp) :: m1,m2
  real(sp) :: m1tsq,m2tsq
  real(sp) :: m1t,m2t

  ! photon transverse momentum, angle
  real(sp) :: ktmin,ktmax
  real(sp) :: k1t,kat
  real(sp) :: k2t,kbt
  real(sp) :: thk1,thka
  real(sp) :: thk2,thkb
  real(sp) :: k1sq,kasq,k2sq,kbsq
  real(sp) :: k1ka,k2kb
  real(sp) :: dotcross

  ! imbalances, qT, pT, angles
  real(sp) :: qt,thqt
  real(sp) :: pt,thpt

  ! momentum fraction
  real(sp) :: xa,xb
  ! Mandelstam
  real(sp) :: mans,mant,manu

  ! nucleus info
  integer(i2b) :: atomA,atomZ
  real(sp) :: RA,rho0
  real(sp), parameter :: a0 = 1d0/sqrt(2d0) /gevfm !~0.71 fm

  ! impact parameter
  real(sp) :: bpmin,bp,bpmax

  ! acoplanarity cut
  real(sp) :: aco,acocut
  ! asymmetry cut
  real(sp) :: asym,asymcut
  ! kt cut
  real(sp) :: kt,ktcut

  ! masses (GeV)
  integer(i1b) :: Mli
  real(sp), parameter :: M_pro = 0.93827208816d0
  real(sp), parameter :: M_neu = 0.93956542052d0
  real(sp), parameter :: M_ele = 0.00051d0
  real(sp), parameter :: M_mu  = 0.10566d0
  real(sp), parameter :: M_tau = 1.77686d0
  real(sp) :: Mllmin,Mll,Mllmax

  ! coupling 
  real(sp), parameter :: alphae = 1d0/137d0
  real(sp) :: alphas

  ! vegas
  integer(i4b) :: ncall1,itmax1
  integer(i4b) :: ncall2,itmax2
  integer(i4b) :: init,ndim,nprn,dind
  real(sp) :: avgi,chi2a,sd
  real(sp), dimension(30) :: region

  ! event counters
  integer(8) :: totevnt,accevnt

  ! io-file name
  character(len=*), parameter :: ifile='input.dat'
  character(len=*), parameter :: ofile='output.dat'

  ! flags
  logical :: fill_hist
  logical :: formfac
  logical :: debug

end module comvar
