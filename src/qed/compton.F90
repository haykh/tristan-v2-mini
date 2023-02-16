module m_compton

  use m_globalnamespace
  use m_qednamespace
  use m_aux
  use m_errors
  use m_bincoupling
  use m_particlelogistics
#ifdef COMPTONSCATTERING
  implicit none

  real(kind=8), parameter :: low_eph_lim = 2d-3
  real, parameter :: p_lim = 0.1
  integer, private :: current_step

  !--- PRIVATE variables/functions -------------------------------!
  private :: comptonOnTile_bin, comptonOnTile_mc, scatterPhoton
  private :: generateRandomCosTheta, du_KN_Newt, boostPhoton
  private :: low_eph_lim
  !...............................................................!
contains
  subroutine comptonScattering(timestep)
    implicit none
    integer, intent(in) :: timestep
    integer :: compton_species_1(20), compton_species_2(20)
    integer :: s, si_1, si_2
    integer :: ti, tj, tk

    current_step = timestep

    ! find all the species that participate in Compton scattering
    si_1 = 0; si_2 = 0
    do s = 1, nspec
      ! go through the list and assign species either to
      ! group1 (electrons/positrons) or to group2 (photons):
      if (species(s) % compton_sp) then
        ! electrons or positrons:
        if ((species(s) % m_sp .eq. 1.0) .and. (abs(species(s) % ch_sp) .eq. 1.0)) then
          compton_species_1(si_1 + 1) = s
          si_1 = si_1 + 1
          ! the photons:
        else if ((species(s) % m_sp .eq. 0) .and. (species(s) % ch_sp .eq. 0)) then
          compton_species_2(si_2 + 1) = s
          si_2 = si_2 + 1
        end if
      end if
    end do

    ! proceed if scattering is turned on at least for one electron/positron
    ! species and one photon species:
    if ((si_1 .ne. 0) .and. (si_2 .ne. 0)) then
      ! assuming all tiles are equal
      s = compton_species_1(1)
      ! loop through all the tiles
      do tk = 1, species(s) % tile_nz
        do tj = 1, species(s) % tile_ny
          do ti = 1, species(s) % tile_nx
            call comptonOnTile_mc(ti, tj, tk, &
                                  compton_species_1(1:si_1), si_1, &
                                  compton_species_2(1:si_2), si_2)
          end do
        end do
      end do
    end if
  end subroutine comptonScattering

  subroutine comptonOnTile_bin(ti, tj, tk, &
                               sp_arr_1, n_sp_1, &
                               sp_arr_2, n_sp_2)
    implicit none
    integer, intent(in) :: ti, tj, tk
    integer, intent(in) :: n_sp_1 ! # of species in set
    integer, intent(in) :: sp_arr_1(n_sp_1)  ! el/positrons
    integer, intent(in) :: n_sp_2 ! # of species in set
    integer, intent(in) :: sp_arr_2(n_sp_2) ! photons

    !TODO
  end subroutine comptonOnTile_bin

  subroutine comptonOnTile_mc(ti, tj, tk, &
                              sp_arr_1, n_sp_1, &
                              sp_arr_2, n_sp_2)
    implicit none
    integer, intent(in) :: ti, tj, tk
    integer, intent(in) :: n_sp_1 ! # of species in set
    integer, intent(in) :: sp_arr_1(n_sp_1)  ! el/positrons
    integer, intent(in) :: n_sp_2 ! # of species in set
    integer, intent(in) :: sp_arr_2(n_sp_2) ! photons
    type(couple), allocatable :: el_photon_pairs(:)
    integer :: num_pairs, el_ph, s1, s2, p1, p2
    real :: rnd, P_12
    integer :: tile_x, tile_y, tile_z
    real :: P0, P_corr_el, P_corr_ph, P_el, P_ph, ppt0
    real :: P_max_el, P_max_ph
    logical :: KleinNishina
    integer :: num_1, num_2, num_iter, iter, num_scatter_max
    real(kind=8) :: el_gamma, pel_x, pel_y, pel_z
    real(kind=8) :: eph, kph_x, kph_y, kph_z
    real(kind=8) :: eph_RF, kph_RF_x, kph_RF_y, kph_RF_z
    real, pointer :: u_el, v_el, w_el, wei_el
    real, pointer :: u_ph, v_ph, w_ph, wei_ph
    real :: u_el_new, v_el_new, w_el_new
    real :: u_ph_new, v_ph_new, w_ph_new
    real :: wei_split, wei_1, wei_2, wei_split_tot
    real :: wei_split_el, wei_split_ph
    real, parameter :: p_max = 1.0 + TINYREAL

    ! couple the electrons/positrons (group1) and photons (group2):
    call coupleParticlesOnTile(ti, tj, tk, sp_arr_1, n_sp_1, sp_arr_2, n_sp_2, &
                               el_photon_pairs, num_pairs, num_1, num_2, wei_1, wei_2)
    if (num_pairs .eq. 0) return
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! calculate prob. correction factor:
    tile_x = species(1) % prtl_tile(ti, tj, tk) % x2 - &
             species(1) % prtl_tile(ti, tj, tk) % x1
    tile_y = species(1) % prtl_tile(ti, tj, tk) % y2 - &
             species(1) % prtl_tile(ti, tj, tk) % y1
    tile_z = species(1) % prtl_tile(ti, tj, tk) % z2 - &
             species(1) % prtl_tile(ti, tj, tk) % z1
    ! reference # pairs on a tile:
    ppt0 = ppc0 * REAL(tile_x * tile_y * tile_z)
    ! make it independent of ppt0 & qed step:
    P0 = QED_tau0 * REAL(Compton_interval) * CC / ppt0
    wei_split_tot = 0.0
    num_scatter_max = max(num_1, num_2)
    do el_ph = 1, num_scatter_max
      wei_split_tot = wei_split_tot + min(el_photon_pairs(el_ph) % part_1 % wei, &
                                          el_photon_pairs(el_ph) % part_2 % wei)
    end do
    ! correction to match binary pairing, correct for non-integer weights, etc:
    P0 = P0 * wei_1 * wei_2 / wei_split_tot
    ! the correction factors for a paired photon and electron/positron:
    P_corr_ph = P0
    P_corr_el = P0 * Compton_nph_over_ne
    if (.not. Compton_el_recoil) P_corr_el = P_corr_ph
    if (2.0 * max(P_corr_ph, P_corr_el) .le. p_lim) then
      ! at most one candidate scattering per particle for the larger set:
      !num_scatter_max = CEILING(2.0 * REAL(num_scatter_max) * max(P_corr_ph, P_corr_el) + TINYREAL)
      !num_scatter_max = min(num_scatter_max, max(num_1, num_2))
      num_iter = 1
    else
      ! multiple candidate scatterings per particle:
      num_iter = CEILING(2.0 * max(P_corr_ph, P_corr_el) / p_lim + TINYREAL)
      num_iter = min(num_iter, min(num_1, num_2))
    end if
    P0 = wei_split_tot
    wei_split_tot = 0.0
    ! additional pairings are generated by cyclically shifting the particle list of the smaller set:
    if (num_1 .le. num_2) then
      do iter = 0, num_iter - 1
        do el_ph = 1, num_scatter_max
          wei_split_tot = wei_split_tot + min(el_photon_pairs(modulo(el_ph - 1 + iter, num_1) + 1) % part_1 % wei, &
                                              el_photon_pairs(el_ph) % part_2 % wei)
        end do
      end do
    else
      do iter = 0, num_iter - 1
        do el_ph = 1, num_scatter_max
          wei_split_tot = wei_split_tot + min(el_photon_pairs(el_ph) % part_1 % wei, &
                                              el_photon_pairs(modulo(el_ph - 1 + iter, num_2) + 1) % part_2 % wei)
        end do
      end do
    end if
    P0 = P0 / wei_split_tot
    P_corr_ph = P_corr_ph * P0
    P_corr_el = P_corr_el * P0
    P_max_ph = 2.0 * P_corr_ph
    P_max_el = 2.0 * P_corr_el
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#ifdef DEBUG
    if ((P_max_ph .lt. 0) .or. (P_max_ph .gt. p_max)) then
      print *, 'max(P_ph) = ', P_max_ph
      call throwError('Maximum photon Compton scattering probability out of bounds!')
    end if
    if (((P_max_el .lt. 0) .or. (P_max_el .gt. p_max)) .and. Compton_el_recoil) then
      print *, 'max(P_el) = ', P_max_el
      call throwError('Maximum electron Compton scattering probability out of bounds!')
    end if
#else
    if ((P_max_ph .gt. p_max) .or. ((P_max_el .gt. p_max) .and. Compton_el_recoil)) then
      call addWarning(2)
    end if
#endif

    do iter = 0, num_iter - 1
      do el_ph = 1, num_scatter_max

        ! "extract" the el-photon pair:
        if (num_1 .le. num_2) then
          s1 = el_photon_pairs(modulo(el_ph - 1 + iter, num_1) + 1) % part_1 % spec
          p1 = el_photon_pairs(modulo(el_ph - 1 + iter, num_1) + 1) % part_1 % index
          wei_split_el = el_photon_pairs(modulo(el_ph - 1 + iter, num_1) + 1) % part_1 % wei
          s2 = el_photon_pairs(el_ph) % part_2 % spec
          p2 = el_photon_pairs(el_ph) % part_2 % index
          wei_split_ph = el_photon_pairs(el_ph) % part_2 % wei
        else
          s1 = el_photon_pairs(el_ph) % part_1 % spec
          p1 = el_photon_pairs(el_ph) % part_1 % index
          wei_split_el = el_photon_pairs(el_ph) % part_1 % wei
          s2 = el_photon_pairs(modulo(el_ph - 1 + iter, num_2) + 1) % part_2 % spec
          p2 = el_photon_pairs(modulo(el_ph - 1 + iter, num_2) + 1) % part_2 % index
          wei_split_ph = el_photon_pairs(modulo(el_ph - 1 + iter, num_2) + 1) % part_2 % wei
        end if

#ifdef DEBUG
        if ((species(s1) % m_sp .ne. 1) .or. (abs(species(s1) % ch_sp) .ne. 1)) then
          call throwError('Wrong particle assigned to electron/positron group in Compton particle coupling!')
        end if
        if ((species(s2) % m_sp .ne. 0) .or. (species(s2) % ch_sp .ne. 0)) then
          call throwError('Wrong particle assigned to photon group in Compton particle coupling!')
        end if
#endif

        u_el => species(s1) % prtl_tile(ti, tj, tk) % u(p1)
        v_el => species(s1) % prtl_tile(ti, tj, tk) % v(p1)
        w_el => species(s1) % prtl_tile(ti, tj, tk) % w(p1)
        wei_el => species(s1) % prtl_tile(ti, tj, tk) % weight(p1)
        u_ph => species(s2) % prtl_tile(ti, tj, tk) % u(p2)
        v_ph => species(s2) % prtl_tile(ti, tj, tk) % v(p2)
        w_ph => species(s2) % prtl_tile(ti, tj, tk) % w(p2)
        wei_ph => species(s2) % prtl_tile(ti, tj, tk) % weight(p2)

        pel_x = REAL(u_el, 8); pel_y = REAL(v_el, 8); pel_z = REAL(w_el, 8)
        kph_x = REAL(u_ph, 8); kph_y = REAL(v_ph, 8); kph_z = REAL(w_ph, 8)

        el_gamma = sqrt(1.0d0 + pel_x**2 + pel_y**2 + pel_z**2)
        eph = sqrt(kph_x**2 + kph_y**2 + kph_z**2)

        ! boost photon momentum into electron frame:
        call boostPhoton(el_gamma, pel_x, pel_y, pel_z, &
                         eph, kph_x, kph_y, kph_z, &
                         eph_RF, kph_RF_x, kph_RF_y, kph_RF_z)

        ! compute cross section:
        call computeComptonCrossSection(eph, eph_RF, el_gamma, P_12, KleinNishina)
#ifdef DEBUG
        if ((P_12 .lt. 0.0) .or. (P_12 .gt. 2.0)) then
          print *, 'P_12 = ', P_12
          print *, eph, wei_ph, kph_x, kph_y, kph_z
          call throwError('Compton cross section P_12 out of bounds!')
        end if
#endif

        rnd = random(dseed)
        P_ph = P_corr_ph * P_12
        P_el = P_corr_el * P_12
        if (rnd .le. max(P_ph, P_el)) then
          ! scatter the photon in the electron rest frame:
          call scatterPhoton(KleinNishina, eph_RF, kph_RF_x, kph_RF_y, kph_RF_z)

          ! boost back into lab frame:
          pel_x = -pel_x; pel_y = -pel_y; pel_z = -pel_z
          call boostPhoton(el_gamma, pel_x, pel_y, pel_z, &
                           eph_RF, kph_RF_x, kph_RF_y, kph_RF_z, &
                           eph, kph_x, kph_y, kph_z)
          u_ph_new = REAL(kph_x)
          v_ph_new = REAL(kph_y)
          w_ph_new = REAL(kph_z)
          ! obtain the recoil on the electron via momentum conservation:
          u_el_new = u_el + u_ph - u_ph_new
          v_el_new = v_el + v_ph - v_ph_new
          w_el_new = w_el + w_ph - w_ph_new

          ! within some numeric tolerance, check if the split weight ...
          ! ... matches the original particle weight:
          if (abs(wei_el - wei_split_el) .le. TINYWEI) wei_split_el = wei_el
          if (abs(wei_ph - wei_split_ph) .le. TINYWEI) wei_split_ph = wei_ph

          ! take the smaller of the two weights for the scattering:
          wei_split = min(wei_split_el, wei_split_ph)
#ifdef DEBUG
          if (wei_split .le. TINYWEI) then
            call throwError('ERROR: Weight of to-be splitted particle in Compton <= 0!')
          end if
          if ((wei_split - wei_el .gt. TINYWEI) .or. &
              (wei_split - wei_ph .gt. TINYWEI)) then
            print *, wei_split_el, wei_split_ph, wei_split
            print *, wei_el, wei_ph
            call throwError('ERROR: Weight of to-be splitted particle in Compton exceeds initial particle weight!')
          end if
#endif
          ! el update:
          if (Compton_el_recoil .and. (rnd .le. P_el)) then
            if (abs(wei_el - wei_split) .le. TINYWEI) then
              u_el = u_el_new
              v_el = v_el_new
              w_el = w_el_new
            else ! split electron:
              wei_el = wei_el - wei_split
              ! this is done for safety but is not supposed to happen:
              if (wei_el .le. TINYWEI) then
#ifdef DEBUG
                call throwError('ERROR: Electron weight after splitting in Compton <= 0!')
#endif
                species(s1) % prtl_tile(ti, tj, tk) % proc(p1) = -1
              end if
              call createParticle(s1, species(s1) % prtl_tile(ti, tj, tk) % xi(p1), &
                                  species(s1) % prtl_tile(ti, tj, tk) % yi(p1), &
                                  species(s1) % prtl_tile(ti, tj, tk) % zi(p1), &
                                  species(s1) % prtl_tile(ti, tj, tk) % dx(p1), &
                                  species(s1) % prtl_tile(ti, tj, tk) % dy(p1), &
                                  species(s1) % prtl_tile(ti, tj, tk) % dz(p1), &
                                  u_el_new, v_el_new, w_el_new, weight=wei_split)
            end if
          end if
          ! photon update:
          if (rnd .le. P_ph) then
#ifdef PRTLPAYLOADS
            species(s2) % prtl_tile(ti, tj, tk) % payload2(p2) = REAL(current_step)
            species(s2) % prtl_tile(ti, tj, tk) % payload3(p2) = &
              species(s2) % prtl_tile(ti, tj, tk) % payload3(p2) + 1.0
#endif
            if (abs(wei_ph - wei_split) .le. TINYWEI) then
              u_ph = u_ph_new
              v_ph = v_ph_new
              w_ph = w_ph_new
            else ! split photon:
              wei_ph = wei_ph - wei_split
              if (wei_ph .le. TINYWEI) then
#ifdef DEBUG
                call throwError('ERROR: Photon weight after splitting in Compton <= 0!')
#endif
                species(s2) % prtl_tile(ti, tj, tk) % proc(p2) = -1
              end if
              call createParticle(s2, species(s2) % prtl_tile(ti, tj, tk) % xi(p2), &
                                  species(s2) % prtl_tile(ti, tj, tk) % yi(p2), &
                                  species(s2) % prtl_tile(ti, tj, tk) % zi(p2), &
                                  species(s2) % prtl_tile(ti, tj, tk) % dx(p2), &
                                  species(s2) % prtl_tile(ti, tj, tk) % dy(p2), &
                                  species(s2) % prtl_tile(ti, tj, tk) % dz(p2), &
                                  u_ph_new, v_ph_new, w_ph_new, weight=wei_split)
            end if
          end if
        end if
        u_el => null(); v_el => null(); w_el => null(); wei_el => null()
        u_ph => null(); v_ph => null(); w_ph => null(); wei_ph => null()
      end do
    end do
  end subroutine comptonOnTile_mc

  subroutine computeComptonCrossSection(eph, eph_RF, el_gamma, P_12, KleinNishina)
    implicit none
    real(kind=8), intent(in) :: eph, eph_RF, el_gamma
    real, intent(out) :: P_12
    logical, intent(out) :: KleinNishina
    real(kind=8) :: over_eph_RF, f_KN

    ! note: f_KN is normalized to sigma_T
    if (eph_RF .lt. Thomson_lim) then
      KleinNishina = .false.  ! use classical Thomson cross-section
      f_KN = 1.0d0
    else if (eph_RF .lt. low_eph_lim) then
      ! correctly handle the eph_RF << 1 limit using 2nd order expansion of f_KN:
      KleinNishina = .true.  ! Klein-Nishina
      f_KN = 1.0d0 - 2.0d0 * eph_RF + 5.2d0 * eph_RF**2
    else
      KleinNishina = .true.
      over_eph_RF = 1.0d0 / eph_RF
      f_KN = 0.375d0 * over_eph_RF * ((1.0d0 - 2.0d0 * over_eph_RF - 2.0d0 * over_eph_RF**2) * &
                            log(1.0d0 + 2.0d0 * eph_RF) + 0.5d0 + &
                            4.0d0 * over_eph_RF - 0.5d0 / (1.0d0 + 2.0d0 * eph_RF)**2)
    end if
    ! Cross section in the *lab* frame:
    P_12 = REAL(f_KN * eph_RF / (el_gamma * eph))
  end subroutine computeComptonCrossSection

  subroutine boostPhoton(gam, p_x, p_y, p_z, &
                         eph, k_x, k_y, k_z, &
                         eph1, k1_x, k1_y, k1_z)
    implicit none
    real(kind=8), intent(in) :: gam, p_x, p_y, p_z
    ! the input momentum:
    real(kind=8), intent(in) :: eph, k_x, k_y, k_z
    ! the transformed momentum:
    real(kind=8), intent(out) :: eph1, k1_x, k1_y, k1_z
    real(kind=8) :: p_dot_k, check

    p_dot_k = p_x * k_x + p_y * k_y + p_z * k_z
    ! The transformed photon momentum:
    eph1 = gam * eph - p_dot_k
    k1_x = k_x + (p_dot_k / (1.0d0 + gam) - eph) * p_x
    k1_y = k_y + (p_dot_k / (1.0d0 + gam) - eph) * p_y
    k1_z = k_z + (p_dot_k / (1.0d0 + gam) - eph) * p_z
#ifdef DEBUG
    check = eph1 / (gam * eph)
    if ((check .le. 0.0d0) .or. (check .ge. 2.0d0)) then
      print *, 'Ratio Eph1 / (gamma * Eph)  = ', check
      call throwError('Error in boostPhoton(). Transformed energy out of bounds!')
    end if
#endif
  end subroutine boostPhoton

  subroutine scatterPhoton(KleinNishina, eph_RF, kph_RF_x, kph_RF_y, kph_RF_z)
    implicit none
    logical, intent(in) :: KleinNishina
    real(kind=8), intent(inout) :: eph_RF, kph_RF_x, kph_RF_y, kph_RF_z
    real(kind=8) :: a_RF_x, a_RF_y, a_RF_z
    real(kind=8) :: b_RF_x, b_RF_y, b_RF_z
    real(kind=8) :: c_RF_x, c_RF_y, c_RF_z
    real(kind=8) :: rand_costheta_RF, rand_sintheta_RF
    real(kind=8) :: rand_phi_RF, rand_cosphi_RF, rand_sinphi_RF, norm

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Define a basis in the electron frame:
    norm = 1.0d0 / eph_RF
    a_RF_x = kph_RF_x * norm; a_RF_y = kph_RF_y * norm; a_RF_z = kph_RF_z * norm
    if (a_RF_x .ne. 0.0d0) then
      b_RF_x = -a_RF_y / a_RF_x; b_RF_y = 1.0d0; b_RF_z = 0.0d0
      norm = 1.0d0 / sqrt(b_RF_x**2 + b_RF_y**2)
      b_RF_x = b_RF_x * norm; b_RF_y = b_RF_y * norm
    else
      b_RF_x = 1.0d0; b_RF_y = 0.0d0; b_RF_z = 0.0d0
    end if
    c_RF_x = b_RF_z * a_RF_y - b_RF_y * a_RF_z
    c_RF_y = b_RF_x * a_RF_z - b_RF_z * a_RF_x
    c_RF_z = b_RF_y * a_RF_x - b_RF_x * a_RF_y
    ! note: c_RF is normalized by construction
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    ! Generate random vector in the electron frame ...
    ! ... respecting the differential (Klein-Nishina) cross section ...
    call generateRandomCosTheta(eph_RF, rand_costheta_RF, KleinNishina)
    rand_phi_RF = REAL(2.0 * M_PI * random(dseed), 8)
    rand_sintheta_RF = sqrt(1.0d0 - rand_costheta_RF**2)
    rand_cosphi_RF = cos(rand_phi_RF)
    rand_sinphi_RF = sin(rand_phi_RF)

    ! update the momentum:
    eph_RF = eph_RF / (1.0d0 + eph_RF * (1.0d0 - rand_costheta_RF))
    kph_RF_x = eph_RF * (rand_costheta_RF * a_RF_x + &
                         rand_sintheta_RF * rand_cosphi_RF * b_RF_x + &
                         rand_sintheta_RF * rand_sinphi_RF * c_RF_x)
    kph_RF_y = eph_RF * (rand_costheta_RF * a_RF_y + &
                         rand_sintheta_RF * rand_cosphi_RF * b_RF_y + &
                         rand_sintheta_RF * rand_sinphi_RF * c_RF_y)
    kph_RF_z = eph_RF * (rand_costheta_RF * a_RF_z + &
                         rand_sintheta_RF * rand_cosphi_RF * b_RF_z + &
                         rand_sintheta_RF * rand_sinphi_RF * c_RF_z)
  end subroutine scatterPhoton

  ! sample random u=cos(theta) in the electron frame for the Thomson or
  ! Klein-Nishina differential cross-section. The cross-sections are
  ! for un-polarized photons.
  subroutine generateRandomCosTheta(eph_RF, costheta_RF, KleinNishina)
    implicit none
    real(kind=8), intent(in) :: eph_RF
    real(kind=8), intent(out) :: costheta_RF
    logical, intent(in) :: KleinNishina
    real(kind=8) :: rnd, u, du
    real(kind=8) :: c0, c1, c2, c3, c4
    integer :: iter
    real(kind=8), parameter :: thresh = 1d-7
    integer, parameter :: max_iter = 30
    logical :: converged

    rnd = REAL(random(dseed), 8)
    if (.not. KleinNishina) then
      ! `u` for Thomson can be sampled by transforming `rnd` with
      ! an analytic formula, which is obtained by inverting the
      ! cumulative distribution:
      u = (4.0d0 * rnd - 2.0d0 + sqrt(5.0d0 + 16.0d0 * rnd * (rnd - 1.0d0)))**(1.0d0 / 3.0d0)
      u = u - 1.0d0 / u
    else
      ! generate random costheta for Klein-Nishina
      ! by solving iteratively (via Newton method) for F(u=cos(theta)) = rnd \in [0,1]
      iter = 0
      converged = .false.
      c0 = 1.0d0 + 2.0d0 * eph_RF
      c1 = eph_RF / c0
      c2 = eph_RF**2 - 2.0d0 * eph_RF - 2.0d0
      c3 = eph_RF - 1.0d0 - 0.5d0 * c1**2
      c4 = 1.0d0 / (4.0d0 * eph_RF + 2.0d0 * eph_RF * (1.0d0 + eph_RF) * c1**2 + c2 * log(c0))
      u = 2.0d0 * rnd - 1.0d0
      do while (iter .lt. max_iter)
        iter = iter + 1
        du = du_KN_Newt(eph_RF, u, rnd, c0, c1, c2, c3, c4)
        u = u + du
        if (u .gt. 1.0d0) u = 1.0d0
        if (u .lt. -1.0d0) u = -1.0d0
        if (abs(du) .lt. thresh) then
          converged = .true.
          exit
        end if
      end do
      if (.not. converged) then
        print '(1X,A,F6.3,A,ES10.3,A,F6.3,A)', 'Warning: Cos(theta) = ', u, ' did not converge (eph_RF = ', eph_RF, ', rnd = ', rnd, ')'
#ifdef DEBUG
        call throwError('Random value for cos(theta) in Compton scattering failed to converge.')
#endif
      end if
    end if
    costheta_RF = u
#ifdef DEBUG
    if ((costheta_RF .lt. -1.0d0) .or. (costheta_RF .gt. 1.0d0)) then
      print *, 'Cos(theta) = ', costheta_RF
      call throwError('Cos(theta) for Compton scattering went out of bounds.')
    end if
#endif
  end subroutine generateRandomCosTheta

  ! increment for a Newton root find of F_KN(u=cos(theta)) = rnd
  real(8) function du_KN_Newt(eph_RF, u, rnd, c0, c1, c2, c3, c4)
    implicit none
    real(kind=8), intent(in) :: eph_RF, u, rnd, c0, c1, c2, c3, c4
    real(kind=8) :: Fu, dFdu, g, eg, c0g

    if (eph_RF .lt. low_eph_lim) then
      ! use a 2nd order expansion in eph_RF to avoid numerical issues
      ! when eph_RF << 1:
      Fu = 0.125d0 * (u**3 + 3.0d0 * u + 4.0d0) + &
           0.1875d0 * (u**4 + 2.0d0 * u**2 - 3.0d0) * eph_RF + &
           0.0375d0 * (6.0d0 * u**5 - 5.0d0 * u**4 + 6.0d0 * u**3 - &
                       20.0d0 * u**2 - 12.0d0 * u + 25.0d0) * eph_RF**2
      dFdu = 0.375d0 * (u**2 + 1.0d0) + &
             0.75d0 * (u**3 + u) * eph_RF + &
             0.075d0 * (15.0d0 * u**4 - 10.0d0 * u**3 + 9.0d0 * u**2 - &
                        20.0d0 * u - 6.0d0) * eph_RF**2
    else
      g = 1.0d0 / (1.0d0 + eph_RF * (1.0d0 - u))
      eg = eph_RF * g
      c0g = c0 * g
      Fu = c4 * (c3 + eph_RF * u + 0.5d0 * eg**2 + c0g + c2 * log(c0g))
      dFdu = c4 * (eph_RF + eg**3 + c1 * c0g**2 + c2 * eg)
    end if
    du_KN_Newt = (rnd - Fu) / dFdu
  end function du_KN_Newt

#endif
end module m_compton
