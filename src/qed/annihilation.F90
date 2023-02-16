module m_annihilation

  use m_globalnamespace
  use m_qednamespace
  use m_aux
  use m_errors
  use m_particlelogistics
#ifdef PAIRANNIHILATION

  implicit none

  type, private :: particleID
    ! particle is uniquely identified by its:
    ! ... `s`, `ti, tj, tk` and `p`
    integer :: s, ti, tj, tk, p
    real :: wei = -1
  end type particleID

  type, private :: particleGroup
    integer :: npart
    type(particleID), allocatable :: prtls(:)
  end type particleGroup

  type, private :: particlePair
    type(particleID) :: prtl1, prtl2
  end type particlePair

  !--- PRIVATE variables/functions -------------------------------!
  private :: pairAnnihilationWithGroups, pairAnnihilationWithGroups_mc, &
             breakDownParticles, shuffleGroup, computeAnnihilationCrossSection, &
             annihilatePairs, generateRandomThetaAnn
  !...............................................................!
contains
  subroutine pairAnnihilation()
    implicit none
    type(particleGroup), allocatable :: lec_cell_bins(:, :, :), pos_cell_bins(:, :, :)
    type(particle_tile) :: electrons_tile, positrons_tile
    integer :: ann_electrons(20), ann_positrons(20)
    integer :: s, s_lec, s_pos, npart_lec, npart_pos, s0, si, nn
    integer :: ti, tj, tk, nx_bin, ny_bin, nz_bin, pi, pj, pk
    integer :: p, p_ind, s_ind, x1_tile, y1_tile, z1_tile

    ! find all the species that participate in the annihilation process
    s_lec = 0; s_pos = 0
    do s = 1, nspec
      if (species(s) % annihilation_sp) then
        if (species(s) % ch_sp .lt. 0) then
          ann_electrons(s_lec + 1) = s
          s_lec = s_lec + 1
        else if (species(s) % ch_sp .gt. 0) then
          ann_positrons(s_pos + 1) = s
          s_pos = s_pos + 1
        else
          call throwError("Only charged particles can participate in annihilation.")
        end if
      end if
    end do

    if ((s_lec * s_pos .eq. 0) .and. (s_lec + s_pos .ne. 0)) then
      call throwError("Pair annihilation requires at least 2 species with opposite charges to participate.")
    end if

    if (s_lec * s_pos .ne. 0) then
      ! this assumes that all the species have the same tiles
      ! ... so picking just a random `s0`
      s0 = ann_electrons(1)
      ! loop over all tiles
      do tk = 1, species(s0) % tile_nz
        do tj = 1, species(s0) % tile_ny
          do ti = 1, species(s0) % tile_nx
            if ((Annihilation_sporadic) .and. (random(dseed) * Annihilation_interval .gt. 1.0)) then
              cycle
            end if

            ! number of cells on a tile
            nx_bin = species(s0) % prtl_tile(ti, tj, tk) % x2 - species(s0) % prtl_tile(ti, tj, tk) % x1
            ny_bin = species(s0) % prtl_tile(ti, tj, tk) % y2 - species(s0) % prtl_tile(ti, tj, tk) % y1
            nz_bin = species(s0) % prtl_tile(ti, tj, tk) % z2 - species(s0) % prtl_tile(ti, tj, tk) % z1

            ! count total # of electrons in the tile
            npart_lec = 0
            do si = 1, s_lec
              npart_lec = npart_lec + species(ann_electrons(si)) % prtl_tile(ti, tj, tk) % npart_sp
            end do
            if (allocated(lec_cell_bins)) deallocate (lec_cell_bins)
            allocate (lec_cell_bins(nx_bin, ny_bin, nz_bin))

            if (npart_lec .eq. 0) then
              cycle
            end if

            ! count total # of positrons in the tile
            npart_pos = 0
            do si = 1, s_pos
              npart_pos = npart_pos + species(ann_positrons(si)) % prtl_tile(ti, tj, tk) % npart_sp
            end do
            if (allocated(pos_cell_bins)) deallocate (pos_cell_bins)
            allocate (pos_cell_bins(nx_bin, ny_bin, nz_bin))

            if (npart_pos .eq. 0) then
              cycle
            end if

            ! initialize cell-based bins
            do pk = 1, nz_bin
              do pj = 1, ny_bin
                do pi = 1, nx_bin
                  lec_cell_bins(pi, pj, pk) % npart = 0
                  allocate (lec_cell_bins(pi, pj, pk) % prtls(npart_lec))
                  pos_cell_bins(pi, pj, pk) % npart = 0
                  allocate (pos_cell_bins(pi, pj, pk) % prtls(npart_pos))
                end do
              end do
            end do

            ! put particles into correct groups according to their cells
            x1_tile = species(s0) % prtl_tile(ti, tj, tk) % x1
            y1_tile = species(s0) % prtl_tile(ti, tj, tk) % y1
            z1_tile = species(s0) % prtl_tile(ti, tj, tk) % z1
            do si = 1, s_lec
              do p = 1, species(ann_electrons(si)) % prtl_tile(ti, tj, tk) % npart_sp
                pi = species(ann_electrons(si)) % prtl_tile(ti, tj, tk) % xi(p) - x1_tile + 1
                pj = species(ann_electrons(si)) % prtl_tile(ti, tj, tk) % yi(p) - y1_tile + 1
                pk = species(ann_electrons(si)) % prtl_tile(ti, tj, tk) % zi(p) - z1_tile + 1

                lec_cell_bins(pi, pj, pk) % npart = lec_cell_bins(pi, pj, pk) % npart + 1
                nn = lec_cell_bins(pi, pj, pk) % npart
                lec_cell_bins(pi, pj, pk) % prtls(nn) % s = ann_electrons(si)
                lec_cell_bins(pi, pj, pk) % prtls(nn) % ti = ti
                lec_cell_bins(pi, pj, pk) % prtls(nn) % tj = tj
                lec_cell_bins(pi, pj, pk) % prtls(nn) % tk = tk
                lec_cell_bins(pi, pj, pk) % prtls(nn) % p = p
              end do
            end do

            do si = 1, s_pos
              do p = 1, species(ann_positrons(si)) % prtl_tile(ti, tj, tk) % npart_sp
                pi = species(ann_positrons(si)) % prtl_tile(ti, tj, tk) % xi(p) - x1_tile + 1
                pj = species(ann_positrons(si)) % prtl_tile(ti, tj, tk) % yi(p) - y1_tile + 1
                pk = species(ann_positrons(si)) % prtl_tile(ti, tj, tk) % zi(p) - z1_tile + 1

                pos_cell_bins(pi, pj, pk) % npart = pos_cell_bins(pi, pj, pk) % npart + 1
                nn = pos_cell_bins(pi, pj, pk) % npart
                pos_cell_bins(pi, pj, pk) % prtls(nn) % s = ann_positrons(si)
                pos_cell_bins(pi, pj, pk) % prtls(nn) % ti = ti
                pos_cell_bins(pi, pj, pk) % prtls(nn) % tj = tj
                pos_cell_bins(pi, pj, pk) % prtls(nn) % tk = tk
                pos_cell_bins(pi, pj, pk) % prtls(nn) % p = p
              end do
            end do

            ! at this point positrons and electrons on a tile are distributed ...
            ! ... into groups based on their cells
            ! loop over all cells on a tile
            do pk = 1, nz_bin
              do pj = 1, ny_bin
                do pi = 1, nx_bin
                  if ((lec_cell_bins(pi, pj, pk) % npart .ge. 1) .and. (pos_cell_bins(pi, pj, pk) % npart .ge. 1)) then
                    call pairAnnihilationWithGroups(lec_cell_bins(pi, pj, pk), pos_cell_bins(pi, pj, pk))
                  end if
                end do
              end do
            end do
          end do
        end do
      end do
    end if
  end subroutine pairAnnihilation

  subroutine pairAnnihilationWithGroups(electron_group, positron_group)
    implicit none
    type(particleGroup), intent(in) :: electron_group, positron_group
    if (Annihilation_algorithm .eq. 1) then
      call throwError('Annihilation algorithm currently only supports MC pairing.')
    else if (Annihilation_algorithm .eq. 2) then
      call pairAnnihilationWithGroups_mc(electron_group, positron_group)
    else
      call throwError('Unrecognized annihilation algorithm.')
    end if
  end subroutine pairAnnihilationWithGroups

  subroutine pairAnnihilationWithGroups_mc(electron_group, positron_group)
    implicit none
    type(particleGroup), intent(in) :: electron_group, positron_group
    type(particleGroup) :: splitted_lec_group, splitted_pos_group
    real :: lec_weight, pos_weight
    type(particlePair), allocatable :: ep_pairs(:)
    integer :: s1, s2, p1, p2
    integer :: num_couples, n, ti1, tj1, tk1, ti2, tj2, tk2
    real :: P_12, wei1, wei2, P_corr, wei_split_tot, rnd

    call breakDownParticles(electron_group, splitted_lec_group, lec_weight)
    call breakDownParticles(positron_group, splitted_pos_group, pos_weight)
    call shuffleGroup(splitted_lec_group)
    call shuffleGroup(splitted_pos_group)

    num_couples = min(splitted_lec_group % npart, splitted_pos_group % npart)
    allocate (ep_pairs(num_couples))

    ! couple particles
    wei_split_tot = 0.0
    do n = 1, num_couples
      ep_pairs(n) % prtl1 = splitted_lec_group % prtls(n)
      ep_pairs(n) % prtl2 = splitted_pos_group % prtls(n)
      wei_split_tot = wei_split_tot + min(ep_pairs(n) % prtl1 % wei, ep_pairs(n) % prtl2 % wei)
    end do

    ! make it independent of ppc0 & qed step:
    P_corr = (3.0 / 8.0) * QED_tau0 * REAL(Annihilation_interval) * CC / ppc0
    ! match with binary pairing:
    P_corr = P_corr * max(lec_weight, pos_weight)
    ! `min(wei_1, wei_2)` is what could be scattered in an ideal pairing world, ...
    ! ... `wei_split_tot` is what is actually available to scatter due to ...
    ! ... non-ideal pairing:
    P_corr = P_corr * min(lec_weight, pos_weight) / wei_split_tot

    do n = 1, num_couples
      ! compute P_12 for each pair of e+e-
      call computeAnnihilationCrossSection(ep_pairs(n), P_12)
      ! to match the optical depth with the binary pairing case:
      P_12 = P_12 * P_corr
#ifdef DEBUG
      if ((P_12 .lt. 0.0) .or. (P_12 .gt. 1.0)) then
        print *, 'P_12 = ', P_12
        call throwError('Annihilation cross section P_12 out of bounds!')
      end if
#else
      if ((P_12 .gt. 1.0)) then
        call addWarning(3)
      end if
#endif
      rnd = random(dseed)
      if (rnd .le. P_12) then
        ! pair produce
        call annihilatePairs(ep_pairs(n))
        ! schedule particles for deletion
        s1 = ep_pairs(n) % prtl1 % s
        p1 = ep_pairs(n) % prtl1 % p
        ti1 = ep_pairs(n) % prtl1 % ti
        tj1 = ep_pairs(n) % prtl1 % tj
        tk1 = ep_pairs(n) % prtl1 % tk
        wei1 = ep_pairs(n) % prtl1 % wei
        s2 = ep_pairs(n) % prtl2 % s
        p2 = ep_pairs(n) % prtl2 % p
        ti2 = ep_pairs(n) % prtl2 % ti
        tj2 = ep_pairs(n) % prtl2 % tj
        tk2 = ep_pairs(n) % prtl2 % tk
        wei2 = ep_pairs(n) % prtl2 % wei
#ifdef DEBUG
        if ((wei1 .ne. 1) .or. (wei2 .ne. 1)) then
          call throwError("Only particles in groups with weight 1 can annihilate.")
        end if
#endif
        species(s1) % prtl_tile(ti1, tj1, tk1) % weight(p1) = species(s1) % prtl_tile(ti1, tj1, tk1) % weight(p1) - wei1
        species(s2) % prtl_tile(ti2, tj2, tk2) % weight(p2) = species(s2) % prtl_tile(ti2, tj2, tk2) % weight(p2) - wei2
        if (species(s1) % prtl_tile(ti1, tj1, tk1) % weight(p1) .lt. 1e-6) then
          species(s1) % prtl_tile(ti1, tj1, tk1) % proc(p1) = -1
        end if
        if (species(s2) % prtl_tile(ti2, tj2, tk2) % weight(p2) .lt. 1e-6) then
          species(s2) % prtl_tile(ti2, tj2, tk2) % proc(p2) = -1
        end if
      end if
    end do
  end subroutine pairAnnihilationWithGroups_mc

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! . . . . Physics functions . . . .
  subroutine computeAnnihilationCrossSection(ep_pair, P_12)
    implicit none
    type(particlePair), intent(in) :: ep_pair
    real, intent(out) :: P_12
    integer :: s1, s2, p1, p2
    integer :: ti1, tj1, tk1, ti2, tj2, tk2
    real(kind=8) :: sigma_ann
    real(kind=8) :: lec_u, lec_v, lec_w, pos_u, pos_v, pos_w, lec_gamma, pos_gamma
    real(kind=8) :: COM_beta_u, COM_beta_v, COM_beta_w
    real(kind=8) :: COM_u, COM_v, COM_w, COM_gamma
    real(kind=8) :: gamma_2, gamma_2_sqr, u_2
    real :: wei1, wei2

    s1 = ep_pair % prtl1 % s
    p1 = ep_pair % prtl1 % p
    ti1 = ep_pair % prtl1 % ti
    tj1 = ep_pair % prtl1 % tj
    tk1 = ep_pair % prtl1 % tk
    wei1 = ep_pair % prtl1 % wei
    s2 = ep_pair % prtl2 % s
    p2 = ep_pair % prtl2 % p
    ti2 = ep_pair % prtl2 % ti
    tj2 = ep_pair % prtl2 % tj
    tk2 = ep_pair % prtl2 % tk
    wei2 = ep_pair % prtl2 % wei

#ifdef DEBUG
    if ((wei1 .ne. wei2) .or. (wei1 .ne. 1.0)) then
      call throwError('Unequal weights in `computeAnnihilationCrossSection`: '//STR(wei1)//':'//STR(wei2))
    end if
#endif

    lec_u = REAL(species(s1) % prtl_tile(ti1, tj1, tk1) % u(p1), 8)
    lec_v = REAL(species(s1) % prtl_tile(ti1, tj1, tk1) % v(p1), 8)
    lec_w = REAL(species(s1) % prtl_tile(ti1, tj1, tk1) % w(p1), 8)
    lec_gamma = sqrt(1d0 + lec_u**2 + lec_v**2 + lec_w**2)
    pos_u = REAL(species(s2) % prtl_tile(ti2, tj2, tk2) % u(p2), 8)
    pos_v = REAL(species(s2) % prtl_tile(ti2, tj2, tk2) % v(p2), 8)
    pos_w = REAL(species(s2) % prtl_tile(ti2, tj2, tk2) % w(p2), 8)
    pos_gamma = sqrt(1d0 + pos_u**2 + pos_v**2 + pos_w**2)

    gamma_2 = lec_gamma * pos_gamma - (lec_u * pos_u + lec_v * pos_v + lec_w * pos_w)
    gamma_2_sqr = gamma_2 * gamma_2

    ! using assymptotic relations for `gamma_2 >> 1` and `gamma_2 ~ 1` ...
    ! ... with an error of <0.01%
    if (gamma_2 .lt. 1.01) then
      sigma_ann = 1d0 / sqrt(1d0 - 1d0 / gamma_2_sqr)
    else if (gamma_2 .gt. 100) then
      sigma_ann = ((log(2d0 * gamma_2) - 1d0) / gamma_2) + ((3d0 * log(2d0 * gamma_2) - 2d0) / gamma_2_sqr)
    else
      u_2 = sqrt(gamma_2_sqr - 1d0)
      sigma_ann = ((gamma_2_sqr + 4d0 * gamma_2 + 1d0) * log(gamma_2 + u_2) / (gamma_2_sqr - 1d0) - &
                   (gamma_2 + 3d0) / u_2) / (gamma_2 + 1d0)
    end if

    ! take into account the relative velocity
    P_12 = sigma_ann * REAL(sqrt(gamma_2**2 - 1d0) / (lec_gamma * pos_gamma))
  end subroutine computeAnnihilationCrossSection

  subroutine annihilatePairs(ep_pair)
    implicit none
    type(particlePair), intent(in) :: ep_pair
    integer :: s1, s2, p1, p2
    integer :: ti1, tj1, tk1, ti2, tj2, tk2
    real(kind=8) :: dummy
    real(kind=8) :: x_ph, y_ph, z_ph, lec_x, lec_y, lec_z, pos_x, pos_y, pos_z
    real(kind=8) :: lec_Ux, lec_Uy, lec_Uz, pos_Ux, pos_Uy, pos_Uz, lec_gamma, pos_gamma
    real(kind=8) :: beta_CM_x, beta_CM_y, beta_CM_z, beta_CM_sq, gamma_inCM
    real(kind=8) :: a_x, a_y, a_z                     ! CoM basis vector along momentum
    real(kind=8) :: b_x, b_y, b_z, c_x, c_y, c_z      ! CoM basis vector perp to momentum
    real(kind=8) :: rand_theta_CM, rand_phi_CM
    real(kind=8) :: cos_rand_theta_CM, cos_rand_phi_CM, sin_rand_theta_CM, sin_rand_phi_CM
    real(kind=8) :: kprime_x, kprime_y, kprime_z
    real(kind=8) :: eph_prime, k_x, k_y, k_z, eph

    s1 = ep_pair % prtl1 % s
    p1 = ep_pair % prtl1 % p
    ti1 = ep_pair % prtl1 % ti
    tj1 = ep_pair % prtl1 % tj
    tk1 = ep_pair % prtl1 % tk
    s2 = ep_pair % prtl2 % s
    p2 = ep_pair % prtl2 % p
    ti2 = ep_pair % prtl2 % ti
    tj2 = ep_pair % prtl2 % tj
    tk2 = ep_pair % prtl2 % tk

    ! Photons will be put in the center of mass of two particles
    lec_x = REAL(species(s1) % prtl_tile(ti1, tj1, tk1) % dx(p1), 8) + REAL(species(s1) % prtl_tile(ti1, tj1, tk1) % xi(p1), 8)
    lec_y = REAL(species(s1) % prtl_tile(ti1, tj1, tk1) % dy(p1), 8) + REAL(species(s1) % prtl_tile(ti1, tj1, tk1) % yi(p1), 8)
    lec_z = REAL(species(s1) % prtl_tile(ti1, tj1, tk1) % dz(p1), 8) + REAL(species(s1) % prtl_tile(ti1, tj1, tk1) % zi(p1), 8)
    pos_x = REAL(species(s2) % prtl_tile(ti2, tj2, tk2) % dx(p2), 8) + REAL(species(s2) % prtl_tile(ti2, tj2, tk2) % xi(p2), 8)
    pos_y = REAL(species(s2) % prtl_tile(ti2, tj2, tk2) % dy(p2), 8) + REAL(species(s2) % prtl_tile(ti2, tj2, tk2) % yi(p2), 8)
    pos_z = REAL(species(s2) % prtl_tile(ti2, tj2, tk2) % dz(p2), 8) + REAL(species(s2) % prtl_tile(ti2, tj2, tk2) % zi(p2), 8)

#if defined (oneD) || defined (twoD) || defined(threeD)
    x_ph = 0.5 * (lec_x + pos_x)
#else
    call throwError('ERROR: No dimension specified.')
#endif

#if defined (twoD) || defined(threeD)
    y_ph = 0.5 * (lec_y + pos_y)
#else
    y_ph = 0.5
#endif

#if defined(threeD)
    z_ph = 0.5 * (lec_z + pos_z)
#else
    z_ph = 0.5
#endif

    ! Compensate for currents
    dummy = 1.0d0 / REAL(species(s1) % prtl_tile(ti1, tj1, tk1) % weight(p1), 8)
    call depositCurrentsFromSingleParticle(s1, species(s1) % prtl_tile(ti1, tj1, tk1), p1, &
                                           REAL(lec_x), REAL(lec_y), REAL(lec_z), &
                                           REAL(x_ph), REAL(y_ph), REAL(z_ph), &
                                           multiplier=REAL(dummy))
    dummy = 1.0d0 / REAL(species(s2) % prtl_tile(ti2, tj2, tk2) % weight(p2), 8)
    call depositCurrentsFromSingleParticle(s2, species(s2) % prtl_tile(ti2, tj2, tk2), p2, &
                                           REAL(pos_x), REAL(pos_y), REAL(pos_z), &
                                           REAL(x_ph), REAL(y_ph), REAL(z_ph), &
                                           multiplier=REAL(dummy))

    lec_Ux = REAL(species(s1) % prtl_tile(ti1, tj1, tk1) % u(p1), 8)
    lec_Uy = REAL(species(s1) % prtl_tile(ti1, tj1, tk1) % v(p1), 8)
    lec_Uz = REAL(species(s1) % prtl_tile(ti1, tj1, tk1) % w(p1), 8)
    lec_gamma = sqrt(1d0 + lec_Ux**2 + lec_Uy**2 + lec_Uz**2)
    pos_Ux = REAL(species(s2) % prtl_tile(ti2, tj2, tk2) % u(p2), 8)
    pos_Uy = REAL(species(s2) % prtl_tile(ti2, tj2, tk2) % v(p2), 8)
    pos_Uz = REAL(species(s2) % prtl_tile(ti2, tj2, tk2) % w(p2), 8)
    pos_gamma = sqrt(1d0 + pos_Ux**2 + pos_Uy**2 + pos_Uz**2)

    dummy = lec_gamma + pos_gamma
    beta_CM_x = (lec_Ux + pos_Ux) / dummy
    beta_CM_y = (lec_Uy + pos_Uy) / dummy
    beta_CM_z = (lec_Uz + pos_Uz) / dummy
    call LorentzBoost(beta_CM_x, beta_CM_y, beta_CM_z, &
                      lec_gamma, &
                      lec_Ux, lec_Uy, lec_Uz, &
                      a_x, a_y, a_z)
    dummy = a_x**2 + a_y**2 + a_z**2
    gamma_inCM = sqrt(1.0d0 + dummy) ! electron/positron energy in CoM frame

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Pick 3 basis vectors in the CoM frame
    ! normalize 1st basis vector
    dummy = sqrt(dummy)
    a_x = a_x / dummy
    a_y = a_y / dummy
    a_z = a_z / dummy
    ! choose 2nd basis vector
    if (a_x .ne. 0.0d0) then
      b_x = -a_y / a_x; b_y = 1.0d0; b_z = 0.0d0
      dummy = sqrt(b_x**2 + b_y**2)
      b_x = b_x / dummy; b_y = b_y / dummy
    else
      b_x = 1.0d0; b_y = 0.0d0; b_z = 0.0d0
    end if
    ! ... 3rd basis vector
    c_x = a_y * b_z - a_z * b_y
    c_y = -a_x * b_z + a_z * b_x
    c_z = a_x * b_y - a_y * b_x
    dummy = sqrt(c_x**2 + c_y**2 + c_z**2)
    c_x = c_x / dummy; c_y = c_y / dummy; c_z = c_z / dummy

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Generate random vector in the CoM frame ...
    ! ... respecting the differential cross section ...
    ! ... at angle `theta` w.r.t. `k1_CM`
    call generateRandomThetaAnn(gamma_inCM, rand_theta_CM)
    rand_phi_CM = 2.0d0 * REAL(M_PI * random(dseed), 8)
    cos_rand_theta_CM = cos(rand_theta_CM)
    sin_rand_theta_CM = sin(rand_theta_CM)
    cos_rand_phi_CM = cos(rand_phi_CM)
    sin_rand_phi_CM = sin(rand_phi_CM)

    ! direction of the newly created photon in CoM frame
    kprime_x = a_x * cos_rand_theta_CM + &
               b_x * sin_rand_theta_CM * cos_rand_phi_CM + &
               c_x * sin_rand_theta_CM * sin_rand_phi_CM
    kprime_y = a_y * cos_rand_theta_CM + &
               b_y * sin_rand_theta_CM * cos_rand_phi_CM + &
               c_y * sin_rand_theta_CM * sin_rand_phi_CM
    kprime_z = a_z * cos_rand_theta_CM + &
               b_z * sin_rand_theta_CM * cos_rand_phi_CM + &
               c_z * sin_rand_theta_CM * sin_rand_phi_CM
#ifdef DEBUG
    if (abs(sqrt(kprime_x**2 + kprime_y**2 + kprime_z**2) - 1.0d0) .gt. 1e-6) then
      call throwError('ERROR: `|k_prime|` is not 1 in annihilatePairs.')
    end if
#endif
    eph_prime = gamma_inCM

    kprime_x = kprime_x * eph_prime
    kprime_y = kprime_y * eph_prime
    kprime_z = kprime_z * eph_prime

    call LorentzBoost(-beta_CM_x, -beta_CM_y, -beta_CM_z, &
                      eph_prime, &
                      kprime_x, kprime_y, kprime_z, &
                      k_x, k_y, k_z)
    call injectParticleLocally(Annihilation_photon_sp, REAL(x_ph), REAL(y_ph), REAL(z_ph), &
                               REAL(k_x), REAL(k_y), REAL(k_z))
    call LorentzBoost(-beta_CM_x, -beta_CM_y, -beta_CM_z, &
                      eph_prime, &
                      -kprime_x, -kprime_y, -kprime_z, &
                      k_x, k_y, k_z)
    call injectParticleLocally(Annihilation_photon_sp, REAL(x_ph), REAL(y_ph), REAL(z_ph), &
                               REAL(k_x), REAL(k_y), REAL(k_z))
  end subroutine annihilatePairs

  subroutine generateRandomThetaAnn(gamma, theta_final)
    implicit none
    real(kind=8), intent(in) :: gamma ! lorentz factor of e+/e- in CoM frame
    real(kind=8), intent(out) :: theta_final
    real(kind=8) :: beta2, beta4, dSigma_dO
    real(kind=8) :: rand_theta, rand_prob, dummy0
    integer :: iter
    iter = 0
    if (gamma .gt. 1.05d0) then
      ! use full expression for dSigma/dO
      beta2 = 1.0d0 - 1.0d0 / gamma**2
      beta4 = beta2**2
      dummy0 = 1.0d0 / (4.0d0 * sqrt(beta2) * gamma**2)
      do while (.true.)
        rand_prob = random(dseed)
        rand_theta = REAL(M_PI * random(dseed), 8)
        dSigma_dO = REAL(dummy0 * (1.0d0 + beta2 * sin(2.0 * rand_theta)**2 - &
                                   beta4 * (1.0d0 - sin(rand_theta)**4)) * sin(rand_theta) / &
                         (1.0d0 - beta2 * cos(rand_theta)**2)**2)
        if (rand_prob .le. dSigma_dO) then
          exit
        end if
        iter = iter + 1
#ifdef DEBUG
        if (iter .gt. 10000) then
          call throwError('Too many iterations in `generateRandomThetaAnn()`.')
        end if
#endif
      end do
    else
      ! use normalized Taylor expansion
      dummy0 = gamma**2 - 1.0d0
      do while (.true.)
        rand_prob = random(dseed)
        rand_theta = REAL(M_PI * random(dseed), 8)
        dSigma_dO = (2.0d0 * (1.0d0 + dummy0 + dummy0 * cos(2.0d0 * rand_theta)) - &
                     dummy0 * cos(4.0d0 * rand_theta)) * &
                    sin(rand_theta) / (2.0d0 - dummy0)
        if (rand_prob .le. dSigma_dO) then
          exit
        end if
        iter = iter + 1
#ifdef DEBUG
        if (iter .gt. 10000) then
          call throwError('Too many iterations in `generateRandomThetaAnn()`.')
        end if
#endif
      end do
    end if
    theta_final = rand_theta
  end subroutine generateRandomThetaAnn

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! . . . . Technical functions . . . .
  subroutine LorentzBoost(L_vx, L_vy, L_vz, &
                          U0_old, &
                          Ux_old, Uy_old, Uz_old, &
                          Ux_new, Uy_new, Uz_new)
    implicit none
    real(kind=8), intent(in) :: L_vx, L_vy, L_vz
    real(kind=8), intent(in) :: U0_old, Ux_old, Uy_old, Uz_old
    real(kind=8), intent(out) :: Ux_new, Uy_new, Uz_new
    real(kind=8) :: L_gamma, dummy, L_vSQR

    L_vSQR = L_vx**2 + L_vy**2 + L_vz**2
    if (L_vSQR .gt. 0.0d0) then
      L_gamma = 1.0d0 / sqrt(1.0d0 - L_vSQR)
      dummy = ((L_gamma - 1.0d0) * (L_vx * Ux_old + L_vy * Uy_old + L_vz * Uz_old) / L_vSQR) - (L_gamma * U0_old)
      Ux_new = Ux_old + dummy * L_vx
      Uy_new = Uy_old + dummy * L_vy
      Uz_new = Uz_old + dummy * L_vz
    else
      Ux_new = Ux_old; Uy_new = Uy_old; Uz_new = Uz_old
    end if
  end subroutine LorentzBoost

  subroutine breakDownParticles(group, set, set_weight)
    implicit none
    type(particleGroup) :: group
    real, intent(out) :: set_weight
    type(particleGroup), intent(out) :: set
    integer :: set_size_, set_size
    type(particleGroup) :: set_
    integer :: s, si, i, p, q, pi, ti, tj, tk
    real :: wei

    ! computing number of particles in the set
    set_size_ = 0
    do pi = 1, group % npart
      s = group % prtls(pi) % s
      ti = group % prtls(pi) % ti
      tj = group % prtls(pi) % tj
      tk = group % prtls(pi) % tk
      p = group % prtls(pi) % p
      set_size_ = set_size_ + CEILING(species(s) % prtl_tile(ti, tj, tk) % weight(p))
    end do
    set_ % npart = 0
    allocate (set_ % prtls(set_size_))
    ! assigning particles in the set
    i = 1
    set_weight = 0.0
    do pi = 1, group % npart
      s = group % prtls(pi) % s
      ti = group % prtls(pi) % ti
      tj = group % prtls(pi) % tj
      tk = group % prtls(pi) % tk
      p = group % prtls(pi) % p
      wei = species(s) % prtl_tile(ti, tj, tk) % weight(p)
      ! distribute the weight so that all weights are 1
      set_weight = set_weight + FLOOR(wei)
      do q = 1, FLOOR(wei)
        set_ % prtls(i) % s = s
        set_ % prtls(i) % ti = group % prtls(pi) % ti
        set_ % prtls(i) % tj = group % prtls(pi) % tj
        set_ % prtls(i) % tk = group % prtls(pi) % tk
        set_ % prtls(i) % p = p
        set_ % prtls(i) % wei = 1.0
        set_ % npart = set_ % npart + 1
        i = i + 1
      end do
    end do
    set_size = i - 1
    allocate (set % prtls(set_size))
    set % npart = set_size
    set % prtls(1:set_size) = set_ % prtls(1:set_size)
  end subroutine breakDownParticles

  ! Knuth's algorithm to randomly shuffle a group
  subroutine shuffleGroup(group)
    implicit none
    type(particleGroup), intent(inout) :: group
    type(particleID) :: temp
    integer :: i, j
    do i = 1, group % npart - 1
      j = randomInt(dseed, i, group % npart + 1)
      temp = group % prtls(i)
      group % prtls(i) = group % prtls(j)
      group % prtls(j) = temp
    end do
  end subroutine shuffleGroup
#endif
end module m_annihilation
