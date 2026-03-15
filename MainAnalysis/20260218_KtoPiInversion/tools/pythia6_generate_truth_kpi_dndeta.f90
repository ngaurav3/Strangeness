program pythia6truth
  implicit double precision(a-h, o-z)
  integer pychge, pycomp
  external pydata, pychge, pycomp

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
  integer :: nargs, nev, iev, i, idc, ks, kf, iq3, apkf, nfinal, kc
  logical :: write_particles, weakdaughter
  integer :: nacc, nch, ncheta05, npipt0405, nkpt0405
  integer :: nch_incl, ncheta05_incl, npipt0405_incl, nkpt0405_incl
  double precision :: ecm, px, py, pz, pt, eta, pabs, kpi, q, ev, mv, kpi_incl

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
  write(20,'(A)') '# event nChEta05 nPiPt0405 nKPt0405 kPiPt0405 ' // &
                  'nChEta05Inclusive nPiPt0405Inclusive nKPt0405Inclusive kPiPt0405Inclusive'
  if (write_particles) then
    open(unit=21, file=trim(particlefile), status='replace', action='write')
    write(21,'(A)') '# PYTHIA6_FINAL_STATE v2 | idx pid charge E Px Py Pz M weakDaughter'
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
    nch_incl = 0
    ncheta05_incl = 0
    npipt0405_incl = 0
    nkpt0405_incl = 0
    kpi_incl = -1d0
    nfinal = 0

    do i = 1, n
      ks = k(i,1)
      if (ks /= 1) cycle
      nfinal = nfinal + 1
      kf = k(i,2)
      iq3 = pychge(kf)
      weakdaughter = has_long_lived_ancestor(i)

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

      if (is_counted_charged_for_activity_pdg(kf)) then
        nch_incl = nch_incl + 1
        if (.not. weakdaughter) nch = nch + 1
        if (dabs(eta) < 0.5d0) then
          ncheta05_incl = ncheta05_incl + 1
          if (.not. weakdaughter) ncheta05 = ncheta05 + 1
        end if
      end if

      apkf = abs(kf)
      if (is_counted_pion_for_ratio_pdg(kf, px, py, pz)) then
        npipt0405_incl = npipt0405_incl + 1
        if (.not. weakdaughter) npipt0405 = npipt0405 + 1
      end if
      if (is_counted_kaon_for_ratio_pdg(kf, px, py, pz)) then
        nkpt0405_incl = nkpt0405_incl + 1
        if (.not. weakdaughter) nkpt0405 = nkpt0405 + 1
      end if
    end do

    if (write_particles) then
      write(21,'(A,I10,1X,A,I10)') '# Event ', nacc, ' N_particles ', nfinal
      do i = 1, n
        ks = k(i,1)
        if (ks /= 1) cycle
        kf = k(i,2)
        iq3 = pychge(kf)
        q = counted_charge_from_pdg(kf)
        ev = p(i,4)
        mv = p(i,5)
        weakdaughter = has_long_lived_ancestor(i)
        write(21,'(I10,1X,I10,1X,F8.3,1X,ES20.10,1X,ES20.10,1X,ES20.10,1X,ES20.10,1X,ES20.10,1X,I2)') &
          i, kf, q, ev, p(i,1), p(i,2), p(i,3), mv, merge(1, 0, weakdaughter)
      end do
    end if

    if (npipt0405 > 0) kpi = dble(nkpt0405) / dble(npipt0405)
    if (npipt0405_incl > 0) kpi_incl = dble(nkpt0405_incl) / dble(npipt0405_incl)
    write(20,'(I10,1X,I10,1X,I10,1X,I10,1X,F12.6,1X,I10,1X,I10,1X,I10,1X,F12.6)') &
      nacc, ncheta05, npipt0405, nkpt0405, kpi, ncheta05_incl, npipt0405_incl, nkpt0405_incl, kpi_incl
  end do

  call pystat(1)
  close(20)
  if (write_particles) close(21)
  write(*,*) 'Wrote truth-level sample to ', trim(outfile)
  if (write_particles) write(*,*) 'Wrote particle-level sample to ', trim(particlefile)
  write(*,*) 'Accepted events: ', nacc
contains
  logical function is_counted_charged_for_activity_pdg(kf) result(flag)
    integer, intent(in) :: kf
    integer :: apkf
    apkf = abs(kf)
    flag = (apkf == 11 .or. apkf == 13 .or. apkf == 15 .or. &
            apkf == 211 .or. apkf == 321 .or. apkf == 2212 .or. &
            apkf == 3112 .or. apkf == 3222 .or. apkf == 3312 .or. &
            apkf == 3334 .or. apkf == 411 .or. apkf == 431 .or. &
            apkf == 521 .or. apkf == 541 .or. apkf == 24)
  end function is_counted_charged_for_activity_pdg

  logical function pass_pid_fiducial_from_mom(px, py, pz) result(flag)
    double precision, intent(in) :: px, py, pz
    double precision :: p2, abscos
    p2 = px*px + py*py + pz*pz
    if (p2 <= 0.d0) then
      flag = .false.
      return
    end if
    abscos = dabs(pz / dsqrt(p2))
    flag = (abscos > 0.15d0 .and. abscos < 0.675d0)
  end function pass_pid_fiducial_from_mom

  double precision function counted_charge_from_pdg(kf) result(q)
    integer, intent(in) :: kf
    q = 0.d0
    select case (kf)
    case (11, 13, 15, 3112, 3312, 3334)
      q = -1.d0
    case (-11, -13, -15, -3112, -3312, -3334)
      q = +1.d0
    case (211, 321, 2212, 3222, 411, 431, 521, 541, 24)
      q = +1.d0
    case (-211, -321, -2212, -3222, -411, -431, -521, -541, -24)
      q = -1.d0
    case default
      q = 0.d0
    end select
  end function counted_charge_from_pdg

  logical function is_counted_pion_for_ratio_pdg(kf, px, py, pz) result(flag)
    integer, intent(in) :: kf
    double precision, intent(in) :: px, py, pz
    double precision :: pt
    pt = dsqrt(px*px + py*py)
    flag = (abs(kf) == 211 .and. pt >= 0.4d0 .and. pt < 5.0d0 .and. &
            pass_pid_fiducial_from_mom(px, py, pz))
  end function is_counted_pion_for_ratio_pdg

  logical function is_counted_kaon_for_ratio_pdg(kf, px, py, pz) result(flag)
    integer, intent(in) :: kf
    double precision, intent(in) :: px, py, pz
    double precision :: pt
    pt = dsqrt(px*px + py*py)
    flag = (abs(kf) == 321 .and. pt >= 0.4d0 .and. pt < 5.0d0 .and. &
            pass_pid_fiducial_from_mom(px, py, pz))
  end function is_counted_kaon_for_ratio_pdg

  recursive logical function has_long_lived_ancestor(idx) result(flag)
    integer, intent(in) :: idx
    integer :: mother, parentkf, parentkc

    if (idx <= 0 .or. idx > n) then
      flag = .false.
      return
    end if

    mother = k(idx,3)
    if (mother <= 0 .or. mother > n) then
      flag = .false.
      return
    end if

    parentkf = k(mother,2)
    parentkc = pycomp(parentkf)
    if (parentkc > 0) then
      if (pmas(parentkc,4) > 10.d0) then
        flag = .true.
        return
      end if
    end if

    flag = has_long_lived_ancestor(mother)
  end function has_long_lived_ancestor
end program pythia6truth
