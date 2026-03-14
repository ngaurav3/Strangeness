program pythia6truth
  implicit double precision(a-h, o-z)
  integer pychge
  external pydata

  integer :: n, npad
  integer :: k(4000,5)
  double precision :: p(4000,5), v(4000,5)
  common /pyjets/ n, npad, k, p, v
  integer :: mstu(200), mstj(200)
  double precision :: paru(200), parj(200)
  common /pydat1/ mstu, paru, mstj, parj
  integer :: kchg(500,4)
  double precision :: pmas(500,4), parf(2000), vckm(4,4)
  common /pydat2/ kchg, pmas, parf, vckm
  integer :: mdcy(500,3), mdme(8000,2), kfdp(8000,5)
  double precision :: brat(8000)
  common /pydat3/ mdcy, mdme, brat, kfdp
  integer :: msel, mselpd, msub(500), kfin(2,-40:40)
  double precision :: ckin(200)
  common /pysubs/ msel, mselpd, msub, kfin, ckin
  integer :: mstp(200), msti(200)
  double precision :: parp(200), pari(200)
  common /pypars/ mstp, parp, msti, pari

  character(len=256) :: arg, outfile, particlefile
  integer :: nargs, nev, iev, i, idc, ks, kf, iq3, apkf, nfinal
  logical :: write_particles
  integer :: nacc, nch, ncheta05, npipt0405, nkpt0405
  double precision :: ecm, px, py, pz, pt, eta, pabs, kpi, q, ev, mv

  nargs = command_argument_count()
  if (nargs < 2) then
    write(*,*) 'Usage: pythia6_generate_truth_kpi_dndeta <nEvents> <summary.txt> [final_state.dat]'
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg,*) nev
  call get_command_argument(2, outfile)
  particlefile = ''
  write_particles = .false.
  if (nargs >= 3) then
    call get_command_argument(3, particlefile)
    write_particles = .true.
  end if

  open(unit=20, file=trim(outfile), status='replace', action='write')
  write(20,'(A)') '# event nChEta05 nPiPt0405 nKPt0405 kPiPt0405'
  if (write_particles) then
    open(unit=21, file=trim(particlefile), status='replace', action='write')
    write(21,'(A)') '# PYTHIA6_FINAL_STATE v1 | idx pid charge E Px Py Pz M'
  end if

  ecm = 91.2d0
  msel = 0
  msub(1) = 1

  do idc = mdcy(23,2), mdcy(23,2) + mdcy(23,3) - 1
    if (abs(kfdp(idc,1)) >= 6) mdme(idc,1) = min(0, mdme(idc,1))
  end do

  call pyinit('CMS', 'e+', 'e-', ecm)

  nacc = 0
  do iev = 1, nev
    call pyevnt
    nacc = nacc + 1

    nch = 0
    ncheta05 = 0
    npipt0405 = 0
    nkpt0405 = 0
    kpi = -1d0
    nfinal = 0

    do i = 1, n
      ks = k(i,1)
      if (ks /= 1) cycle
      nfinal = nfinal + 1
      kf = k(i,2)
      iq3 = pychge(kf)

      px = p(i,1)
      py = p(i,2)
      pz = p(i,3)
      pt = dsqrt(px*px + py*py)
      if (pt <= 0d0) then
        eta = 0d0
      else
        pabs = dsqrt(pt*pt + pz*pz)
        if (pabs > dabs(pz)) then
          eta = 0.5d0 * dlog((pabs + pz) / (pabs - pz))
        else if (pz >= 0d0) then
          eta = 1d6
        else
          eta = -1d6
        end if
      end if

      if (iq3 /= 0) then
        nch = nch + 1
        if (dabs(eta) < 0.5d0) ncheta05 = ncheta05 + 1
      end if

      apkf = abs(kf)
      if (iq3 /= 0 .and. pt >= 0.4d0 .and. pt < 5.0d0) then
        if (apkf == 211) npipt0405 = npipt0405 + 1
        if (apkf == 321) nkpt0405 = nkpt0405 + 1
      end if
    end do

    if (write_particles) then
      write(21,'(A,I10,1X,A,I10)') '# Event ', nacc, ' N_particles ', nfinal
      do i = 1, n
        ks = k(i,1)
        if (ks /= 1) cycle
        kf = k(i,2)
        iq3 = pychge(kf)
        q = dble(iq3) / 3.d0
        ev = p(i,4)
        mv = p(i,5)
        write(21,'(I10,1X,I10,1X,F8.3,1X,ES20.10,1X,ES20.10,1X,ES20.10,1X,ES20.10,1X,ES20.10)') &
          i, kf, q, ev, p(i,1), p(i,2), p(i,3), mv
      end do
    end if

    if (npipt0405 > 0) kpi = dble(nkpt0405) / dble(npipt0405)
    write(20,'(I10,1X,I10,1X,I10,1X,I10,1X,F12.6)') &
      nacc, ncheta05, npipt0405, nkpt0405, kpi
  end do

  call pystat(1)
  close(20)
  if (write_particles) close(21)
  write(*,*) 'Wrote truth-level sample to ', trim(outfile)
  if (write_particles) write(*,*) 'Wrote particle-level sample to ', trim(particlefile)
  write(*,*) 'Accepted events: ', nacc
end program pythia6truth
