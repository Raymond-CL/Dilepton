program main
  use, intrinsic :: iso_fortran_env, only:stdout=>output_unit
  use nrtype
  use nr
  use comvar
  use sysio
  use initial
  use bookkeep
  implicit none
  real :: tstart,tstop
  interface
    function func(dx,wgt)
    use nrtype
    implicit none
    real(sp), dimension(:), intent(in) :: dx
    real(sp), intent(in) :: wgt
    real(sp) :: func
    end function func
  end interface
  call cpu_time(tstart)

  ! read from input file
  call readin()

  ! initialize fixed values
  call initialize()

  ! initialize histogram 
  call initbook()
  call newbook(1,.false.,'dsig/dqt',40,0.0,0.2)
  !call newbook(2,.false.,'dsig/dpt',20,4.0,14.0)
  call newbook(3,.false.,'dsig/dkt',30,0.0,0.15)
  !call newbook(4,.false.,'dsig/dMll',20,4.0,45.0)
  call newbook(5,.false.,'dsig/daco',40,0.0,0.01)
  !call newbook(6,.true.,'dsig/dxa',24,0.0001,0.1)
  !call newbook(7,.true.,'dsig/dxb',24,0.0001,0.1)
  !call newbook(8,.false.,'dsig/dy1',20,-2.5,2.5)
  !call newbook(9,.false.,'dsig/dy2',20,-2.5,2.5)

  ! set integration limits
  call set_limits()

  ! print program info
  call print_prog_info(stdout)

  ! vegas integration (main evaluation procedure)
  init = -1
  fill_hist = .false.
  write(*,*) 'warm up ... ...'
  call vegas(region(1:2*ndim),func,init,ncall1,itmax1,nprn,avgi,sd,chi2a)

  init = +1
  fill_hist = .true.
  write(*,*) 'end run ... ...'
  call vegas(region(1:2*ndim),func,init,ncall2,itmax2,nprn,avgi,sd,chi2a)

  ! print evaluation results and histogram
  call print_results()

  ! end of program
  call cpu_time(tstop)
  write(*,*) 'time elapsed: ',tstop-tstart,'seconds.'
end program main



function func(dx,wgt)
  use nrtype
  use comvar
  use sysio
  use integrand
  use photonff
  use bookkeep
  implicit none
  real(sp), dimension(:), intent(in) :: dx
  real(sp), intent(in) :: wgt
  real(sp) :: func

  if(fill_hist) totevnt = totevnt + 1
  func = 0d0

  p1t = dx(1);  p2t = dx(2)
  thp2 = dx(3)
  y1 = dx(4);   y2 = dx(5)
  k1t = dx(6);  thk1 = dx(7)
  kat = dx(8);  thka = dx(9)
  bp = dx(10)

  ! set kinematics
  call setkinematics()

  ! set cuts
  !if(qt .ge. pt) return
  !if(aco .ge. acocut) return
  !if(asym .ge. asymcut) return
  !if(kt .ge. ktcut) return
  !if((Mll .le. Mllmin) .or. (Mll .ge. Mllmax)) return

  ! calculate diff. X-sectn
  func = diff_xsec()

  ! debug
  if((debug .or. isnan(func)) .and. fill_hist) then
    call print_debug_info()
  endif

  ! fill histogram
  if(fill_hist) then
    accevnt = accevnt + 1
    call fillbook(1, qt, func*wgt/itmax2)
    !call fillbook(2, pt, func*wgt/itmax2)
    call fillbook(3, kt, func*wgt/itmax2)
    !call fillbook(4, Mll, func*wgt/itmax2)
    call fillbook(5, aco, func*wgt/itmax2)
    !call fillbook(6, xa, func*wgt/itmax2)
    !call fillbook(7, xb, func*wgt/itmax2)
    !call fillbook(8, y1, func*wgt/itmax2)
    !call fillbook(9, y2, func*wgt/itmax2)
  endif

  return
end function func
