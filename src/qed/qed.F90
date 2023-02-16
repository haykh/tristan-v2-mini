module m_qedphysics

  use m_aux
  use m_globalnamespace
  use m_qednamespace
  use m_particles
  use m_readinput, only: getInput
  use m_bwpairproduction
  use m_compton
  use m_annihilation
#ifdef QED

  implicit none

  !--- PRIVATE variables/functions -------------------------------!
  !...............................................................!
contains
  subroutine initializeQED()
    implicit none

    call getInput('qed', 'tau0', QED_tau0, 0.1)

#ifdef RADIATION
#ifdef SYNCHROTRON
    QED_tau0 = rad_beta_rec / (cool_gamma_syn**2 * sqrt(sigma) * c_omp)
    if (mpi_rank .eq. 0) then
      print *, 'WARNING: `QED_tau0` is defined via `cool_gamma_syn`; input value is ignored.'
    end if
#endif
#endif

#ifdef BWPAIRPRODUCTION
    call initializeBWPairProduction()
#endif

#ifdef COMPTONSCATTERING
    call initializeComptonScattering()
#endif

#ifdef PAIRANNIHILATION
    call initializePairAnnihilation()
#endif

    call printDiag("initializeQED()", 1)
  end subroutine initializeQED

#ifdef BWPAIRPRODUCTION
  subroutine initializeBWPairProduction()
    implicit none
    integer :: s
    character(len=STR_MAX) :: var_name

    call getInput('bw_pp', 'interval', BW_interval, 1)
    call getInput('bw_pp', 'algorithm', BW_algorithm, 2)
    call getInput('bw_pp', 'electron_sp', BW_electron_sp, 1)
    call getInput('bw_pp', 'positron_sp', BW_positron_sp, 2)
    if ((nspec .lt. BW_electron_sp) .or. (nspec .lt. BW_positron_sp)) then
      call throwError('Number of specified species does not match with `BW_electron_sp` and/or `BW_positron_sp`.')
    end if
    if ((species(BW_electron_sp) % m_sp * species(BW_positron_sp) % m_sp .eq. 0) .or. &
        (species(BW_electron_sp) % ch_sp * species(BW_positron_sp) % ch_sp .eq. 0)) then
      call throwError('BW only produces massive and charged particles.')
    end if
    if ((species(BW_electron_sp) % m_sp .ne. species(BW_positron_sp) % m_sp) .or. &
        (abs(species(BW_electron_sp) % ch_sp) .ne. abs(species(BW_positron_sp) % ch_sp))) then
      call throwError('Masses and charges of `BW_electron_sp` and `BW_positron_sp` have to match.')
    end if

    do s = 1, nspec
      write (var_name, "(A2,I1)") "bw", s
      call getInput('particles', var_name, species(s) % bw_sp, 0)
      if ((species(s) % bw_sp .ne. 0) .and. &
          ((species(s) % ch_sp .ne. 0) .or. (species(s) % m_sp .ne. 0))) then
        call throwError('`ch != 0` or `m != 0` particles cannot pair-produce via BW.')
      end if
      if ((species(s) % bw_sp .gt. 2)) then
        call throwError('Only two BW photon populations are allowed.')
      end if
    end do
  end subroutine initializeBWPairProduction
#endif

#ifdef COMPTONSCATTERING
  subroutine initializeComptonScattering()
    implicit none
    integer :: s
    character(len=STR_MAX) :: var_name
    real :: Pmax

    call getInput('compton', 'interval', Compton_interval, 1)
    call getInput('compton', 'algorithm', Compton_algorithm, 2)
    if (Compton_algorithm .ne. 2) then
      call throwError('Compton scattering currently only supports the MC algorithm.')
    end if
    call getInput('compton', 'el_recoil', Compton_el_recoil, .true.)
    call getInput('compton', 'Thomson_lim', Thomson_lim, 1d-6)
    call getInput('compton', 'nph_over_ne', Compton_nph_over_ne, 1.0)

    do s = 1, nspec
      write (var_name, "(A7,I1)") "compton", s
      call getInput('particles', var_name, species(s) % compton_sp, .true.)
      if (species(s) % compton_sp) then
        if (.not. (((species(s) % m_sp .eq. 0) .and. (species(s) % ch_sp .eq. 0)) .or. &
                   ((species(s) % m_sp .eq. 1.0) .and. (abs(species(s) % ch_sp) .eq. 1.0)))) then
          call throwError('Only electron/positron and photon species can Compton scatter.')
        end if
      end if
    end do
    if (mpi_rank .eq. 0) then
      Pmax = 2.0 * QED_tau0 * REAL(Compton_interval) * CC * max(Compton_nph_over_ne, 1.0)
      print *, '  Reference nph/ne =', Compton_nph_over_ne
      print *, '  Reference max. probability for Compton scattering (MC algorithm) =', Pmax
    end if
  end subroutine initializeComptonScattering
#endif

#ifdef PAIRANNIHILATION
  subroutine initializePairAnnihilation()
    implicit none
    integer :: s
    character(len=STR_MAX) :: var_name

    call getInput('annihilation', 'interval', Annihilation_interval, 1)
    call getInput('annihilation', 'photon_sp', Annihilation_photon_sp, 3)
    call getInput('annihilation', 'algorithm', Annihilation_algorithm, 2)
    call getInput('annihilation', 'sporadic', Annihilation_sporadic, .false.)

    if (Annihilation_algorithm .ne. 2) then
      call throwError('Annihilation currently only supports the MC algorithm.')
    end if

    do s = 1, nspec
      write (var_name, "(A12,I1)") "annihilation", s
      call getInput('particles', var_name, species(s) % annihilation_sp, .false.)
      if (species(s) % annihilation_sp) then
        if ((species(s) % m_sp .eq. 0) .or. (species(s) % ch_sp .eq. 0)) then
          call throwError('Massless/chargeless particles cannot participate in pair annihilation.')
        end if
      end if
    end do
  end subroutine initializePairAnnihilation
#endif

  subroutine QEDstep(timestep)
    implicit none
    integer, intent(in) :: timestep

#ifdef COMPTONSCATTERING
    if (modulo(timestep, Compton_interval) .eq. 0) then
      call comptonScattering(timestep)
    end if
#endif

#ifdef BWPAIRPRODUCTION
    if (modulo(timestep, BW_interval) .eq. 0) then
      call bwPairProduction()
    end if
#endif

#ifdef PAIRANNIHILATION
    if (((.not. Annihilation_sporadic) .and. (modulo(timestep, Annihilation_interval) .eq. 0)) .or. &
        (Annihilation_sporadic)) then
      call pairAnnihilation()
    end if
#endif

    call printDiag("QEDstep()", 2)
  end subroutine QEDstep

#endif
end module m_qedphysics
