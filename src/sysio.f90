module sysio

  use, intrinsic :: iso_fortran_env, &
    only:compiler_version, compiler_options, stdout=>output_unit
  use nrtype
  use comvar
  use bookkeep

contains

  function getu() result(u)
  implicit none
  integer :: u
  logical :: ex
  do u=10,100
    inquire(unit=u,opened=ex)
    if(.not. ex) return
  enddo
  end function getu

  subroutine readin()
  implicit none
  integer :: u,stat
  u = getu()
  open(u,iostat=stat,file=ifile,status='old')
  if(stat.ne.0) call exit(11)
  read(u,*) CME
  read(u,*) ptmin,ptmax
  read(u,*) ylmin,ylmax
  read(u,*) ktmin,ktmax
  read(u,*) Mllmin,Mllmax
  read(u,*) acocut
  read(u,*) asymcut
  read(u,*) ktcut
  read(u,*) ncall1,itmax1
  read(u,*) ncall2,itmax2
  read(u,*) formfac
  read(u,*) debug
  close(u)
  end subroutine

  subroutine print_prog_info(u)
  implicit none
  integer :: u
  write(u,*) '***********************************'
  write(u,*) '     dilepton photo-production     '
  write(u,*) 'Phys. Rev. Lett. 122, 132301 (2019)'
  write(u,*) 'Phys. Rev. D 102, 094013 (2020)    '
  write(u,*) '***********************************'
  write(u,*) 'program compiled by:',compiler_version()
  !write(u,*) 'with compiler options:',compiler_options()
  write(u,*) 'kinematics:'
  write(u,*) '- center of mass energy per nucleon:',CME
  write(u,*) '- lepton pt range:',ptmin, '< pt(GeV)   <',ptmax
  write(u,*) '- lepton yl range:',ylmin, '< yl(GeV)   <',ylmax
  write(u,*) '- photon kt range:',ktmin, '< kt(GeV)   <',ktmax
  write(u,*) '- inv. mass range:',Mllmin,'< Mll(GeV)  <',Mllmax
  write(u,*) '- impact bp range:',bpmin, '< bp(1/GeV) <',bpmax
  write(u,*) '- lepton masses(GeV):',m1,m2
  write(u,*) 'nucleus info:'
  write(u,*) '- mass:',int(atomA,2),', charge:',int(atomZ,2)
  write(u,*) '- Yukawa pot. rng (GeV^-1):',a0
  write(u,*) '- nucleus radius  (GeV^-1):',RA
  write(u,*) '- nucleus density (GeV^+3):',rho0
  write(u,*) 'kinematics cuts:'
  write(u,*) '- acoplanarity cut:',acocut
  write(u,*) '- asymmetry cut   :',asymcut
  write(u,*) '- ATLAS kperp cut :',ktcut
  write(u,*) 'vegas:'
  write(u,*) '- warm-up:',ncall1,' X ',int(itmax1,1)
  write(u,*) '- end-run:',ncall2,' X ',int(itmax2,1)
  if(debug) write(u,*) 'running debug mode'
  write(u,*) '***********************************'
  end subroutine print_prog_info
  
  subroutine print_debug_info()
  implicit none
  integer :: u
  u = stdout
  write(u,*) '=== debugging for event:',totevnt
  write(u,*) '- lepton info:'
  write(u,*) 'lepton pt:',p1t,p2t
  write(u,*) 'lepton th:',thp1,thp2
  write(u,*) 'lepton y :',y1,y2
  write(u,*) 'lepton mt:',m1t,m2t
  write(u,*) 'inv. mass:',mll
  write(u,*) '- momentum fraction and Mandelstam:'
  write(u,*) 'mom frac:',xa,xb
  write(u,*) 'Man var :',mans,mant,manu
  write(u,*) 'stu sum check :',mans+mant+manu
  write(u,*) 'mass sum check:',m1*m1+m2*m2
  write(u,*) '- photon info:'
  write(u,*) 'pt 1a     :',k1t,kat
  write(u,*) 'pt 2b     :',k2t,kbt
  write(u,*) 'angle 1a  :',thk1,thka
  write(u,*) 'angle 2b  :',thk2,thkb
  write(u,*) 'masses 1a :',k1sq,kasq
  write(u,*) 'masses 2b :',k2sq,kbsq
  write(u,*) 'difference:',k1ka,k2kb
  write(u,*) '- imbalance info:'
  write(u,*) 'pair qT  :',qt,thqt
  write(u,*) 'pair PT  :',pt,thpt
  write(u,*) 'kt       :',kt
  write(u,*) 'acoplanar:',aco
  write(u,*) 'asymmetry:',asym
  write(u,*) '- overall factors:'
  write(u,*) 'dcterm:',dotcross
  write(u,*) '=== end of debug for this event',new_line('a')
  !call sleep(1)
  end subroutine print_debug_info

  subroutine print_results()
  implicit none
  integer :: u
  u = stdout  ! print to terminal
  write(u,*) '***********************************'
  write(u,*) 'total generated events:',totevnt
  write(u,*) 'total accepted events :',accevnt
  write(u,*) 'efficiency after cuts :',real(dble(accevnt)/dble(totevnt)*100d0),'%'
  write(u,*) 'total cross-section:',avgi
  write(u,*) 'printing histograms to:',ofile
  write(u,*) '***********************************'
  u = getu()  ! print to file
  open(u,file=ofile,status='replace')
  write(u,*) 'total generated events:',totevnt
  write(u,*) 'total accepted events :',accevnt
  write(u,*) 'efficiency after cuts :',real(dble(accevnt)/dble(totevnt)*100d0),'%'
  write(u,*) 'total cross-section:',avgi
  write(u,*) 'printing histograms:',new_line('a')
  call printbook(u)
  close(u)
  end subroutine print_results

end module sysio
