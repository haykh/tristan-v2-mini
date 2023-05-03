module m_userfile
  use m_globalnamespace
  use m_aux
  use m_readinput
  use m_domain
  use m_particles
  use m_fields
  use m_thermalplasma
  use m_particlelogistics
  implicit none

  !--- PRIVATE variables -----------------------------------------!
  real, private :: Te, Tph0, ph_fraction, t_inject, t_escape
  !...............................................................!

contains
  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
    call getInput('problem', 'Te', Te)
    call getInput('problem', 'Tph0', Tph0)
    call getInput('problem', 'ph_fraction', ph_fraction, 1.0)
    call getInput('problem', 't_inject', t_inject, 0.0)
    call getInput('problem', 't_escape', t_escape, 0.0)
  end subroutine userReadInput

  real function planckSample()
    implicit none
    real :: prob, n
    real :: rnd
    rnd = random(dseed)
    n = 0.0
    prob = 0.0
    do while ((prob .lt. rnd) .and. (n .lt. 40.0))
      n = n + 1.0
      prob = prob + 1.0 / (1.20206 * n**3)
    end do
    planckSample = -log(random(dseed) * random(dseed) * random(dseed) + 1e-16) / n
    return
  end function planckSample

  subroutine injectPhotons(nphotons)
    implicit none
    integer, intent(in) :: nphotons
    integer :: n
    real :: xg, yg, zg, eph, kx, ky, kz
    real :: U_, TH_
    do n = 1, INT(REAL(nphotons) * ph_fraction)
      xg = random(dseed) * (global_mesh % sx)
      yg = random(dseed) * (global_mesh % sy)
      ! zg = random(dseed) * (global_mesh%sz)
      zg = 0.5
      U_ = 2 * (random(dseed) - 0.5)
      TH_ = 2 * M_PI * random(dseed)
      eph = planckSample() * Tph0
      kx = eph * sqrt(1 - U_**2) * cos(TH_)
      ky = eph * sqrt(1 - U_**2) * sin(TH_)
      kz = eph * U_
      call injectParticleGlobally(3, xg, yg, zg, kx, ky, kz, weight=1.0, payload1=0.0)
    end do
  end subroutine injectPhotons

  subroutine removePhotons(tmax)
    implicit none
    real, intent(in) :: tmax
    integer :: s, ti, tj, tk, p
    s = 3

    do ti = 1, species(s) % tile_nx
      do tj = 1, species(s) % tile_ny
        do tk = 1, species(s) % tile_nz
          do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
            if (species(s) % prtl_tile(ti, tj, tk) % payload1(p) .ge. tmax) then
              species(s) % prtl_tile(ti, tj, tk) % proc(p) = -1
            end if
          end do
        end do
      end do
    end do
  end subroutine removePhotons

  subroutine userInitParticles()
    implicit none
    real :: dens
    integer :: ntot
    type(region) :: back_region

    ! the electron and positron thermal background:
    dens = 0.5 * ppc0
    back_region % x_min = 0.0
    back_region % x_max = REAL(global_mesh % sx)
#if defined(twoD) || defined (threeD)
    back_region % y_min = 0.0
    back_region % y_max = REAL(global_mesh % sy)
#endif
#if defined(threeD)
    back_region % z_min = 0.0
    back_region % z_max = REAL(global_mesh % sz)
#endif
    call fillRegionWithThermalPlasma(back_region, (/1, 2/), 2, dens, Te)

    ! init the isotropic photon field:
    if (t_inject .eq. 0.0) then
      ntot = global_mesh % sx * global_mesh % sy * global_mesh % sz * ppc0
      call injectPhotons(ntot)
    end if
  end subroutine userInitParticles

  subroutine userInitFields()
    implicit none
    integer :: i, j, k
    integer :: i_glob, j_glob, k_glob
    ex(:, :, :) = 0; ey(:, :, :) = 0; ez(:, :, :) = 0
    bx(:, :, :) = 0; by(:, :, :) = 0; bz(:, :, :) = 0
    jx(:, :, :) = 0; jy(:, :, :) = 0; jz(:, :, :) = 0
  end subroutine userInitFields
  !............................................................!

  !--- driving ------------------------------------------------!
  subroutine userCurrentDeposit(step)
    implicit none
    integer, optional, intent(in) :: step
  end subroutine userCurrentDeposit

  subroutine userDriveParticles(step)
    implicit none
    integer, optional, intent(in) :: step
  end subroutine userDriveParticles

  subroutine userExternalFields(xp, yp, zp, &
                                ex_ext, ey_ext, ez_ext, &
                                bx_ext, by_ext, bz_ext)
    implicit none
    real, intent(in) :: xp, yp, zp
    real, intent(out) :: ex_ext, ey_ext, ez_ext
    real, intent(out) :: bx_ext, by_ext, bz_ext
    ! some functions of xp, yp, zp
    ex_ext = 0.0; ey_ext = 0.0; ez_ext = 0.0
    bx_ext = 0.0; by_ext = 0.0; bz_ext = 0.0
  end subroutine userExternalFields
  !............................................................!

  !--- boundaries ---------------------------------------------!
  subroutine userParticleBoundaryConditions(step)
    implicit none
    integer, optional, intent(in) :: step
    integer :: ntot

    ! init the isotropic photon field
    if (t_inject .gt. 0.0) then
      ntot = global_mesh % sx * global_mesh % sy * global_mesh % sz * ppc0 / t_inject
      call injectPhotons(ntot)
    end if

    if (t_escape .gt. 0.0) then
      call removePhotons(t_escape)
    end if
  end subroutine userParticleBoundaryConditions

  subroutine userFieldBoundaryConditions(step, updateE, updateB)
    implicit none
    integer, optional, intent(in) :: step
    logical, optional, intent(in) :: updateE, updateB
    logical :: updateE_, updateB_

    if (present(updateE)) then
      updateE_ = updateE
    else
      updateE_ = .true.
    end if

    if (present(updateB)) then
      updateB_ = updateB
    else
      updateB_ = .true.
    end if
  end subroutine userFieldBoundaryConditions
  !............................................................!

  elemental subroutine usrSetPhPld(u0, v0, w0, over_e_temp, &
                                   incr_pld1, incr_pld2, incr_pld3)
    !$omp declare simd(usrSetPhPld)
    real, intent(in) :: u0, v0, w0, over_e_temp
    real, intent(out) :: incr_pld1, incr_pld2, incr_pld3
    incr_pld1 = 1.0; incr_pld2 = 0.0; incr_pld3 = 0.0
  end subroutine

  ! DUMMY
  subroutine writeUsrRestart(rst_file)
    implicit none
    integer, intent(in) :: rst_file
  end subroutine writeUsrRestart

  subroutine readUsrRestart(rst_file)
    implicit none
    integer, intent(in) :: rst_file
  end subroutine readUsrRestart

  subroutine userDeallocate()
    implicit none
  end subroutine userDeallocate

  elemental subroutine usrSetElPld(q_over_m, u0, v0, w0, over_e_temp, &
                                   ex0, ey0, ez0, bx0, by0, bz0, &
                                   incr_pld1, incr_pld2, incr_pld3)
    !$omp declare simd(usrSetElPld)
    real, intent(in) :: q_over_m, u0, v0, w0, over_e_temp, &
                        ex0, ey0, ez0, &
                        bx0, by0, bz0
    real, intent(out) :: incr_pld1, incr_pld2, incr_pld3
    incr_pld1 = 0.0; incr_pld2 = 0.0; incr_pld3 = 0.0
  end subroutine
end module m_userfile
