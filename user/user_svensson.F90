! Configuration for this userfile:
! ```
!   $ python configure.py -2d --user=user_svensson -qed -ann -bw_pp -compton
! ```

module m_userfile
  use m_globalnamespace
  use m_aux
  use m_helpers
  use m_readinput
  use m_domain
  use m_particles
  use m_fields
  use m_thermalplasma
  use m_particlelogistics
  implicit none

  !--- PRIVATE variables -----------------------------------------!
  real, private :: T_init
  logical, private :: inj_mono, inj_photons
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userSpatialDistribution
  !...............................................................!
contains
  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
    call getInput('problem', 'T_init', T_init)
    call getInput('problem', 'inj_mono', inj_mono)
    call getInput('problem', 'inj_photons', inj_photons)
  end subroutine userReadInput

  function userSpatialDistribution(x_glob, y_glob, z_glob, &
                                   dummy1, dummy2, dummy3)
    real :: userSpatialDistribution
    real, intent(in), optional :: x_glob, y_glob, z_glob
    real, intent(in), optional :: dummy1, dummy2, dummy3

    return
  end function

  function userSLBload(x_glob, y_glob, z_glob, &
                       dummy1, dummy2, dummy3)
    real :: userSLBload
    ! global coordinates
    real, intent(in), optional :: x_glob, y_glob, z_glob
    ! global box dimensions
    real, intent(in), optional :: dummy1, dummy2, dummy3
    return
  end function

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

  subroutine generateRandomDirection(kx, ky, kz)
    implicit none
    real, intent(out) :: kx, ky, kz
    real :: rand_costh, rand_phi
    rand_costh = 2.0 * random(dseed) - 1.0
    rand_phi = 2.0 * M_PI * random(dseed)
    kx = sqrt(1.0 - rand_costh**2) * cos(rand_phi)
    ky = sqrt(1.0 - rand_costh**2) * sin(rand_phi)
    kz = rand_costh
  end subroutine generateRandomDirection

  subroutine injectPhoton(step)
    implicit none
    integer, optional, intent(in) :: step
    real :: xg, yg, zg, kx, ky, kz, energy
    real :: ph_temperature
    ph_temperature = T_init

    xg = random(dseed) * REAL(global_mesh % sx)
    yg = random(dseed) * REAL(global_mesh % sy)
    zg = 0.5

    if (inj_mono) then
      energy = ph_temperature
    else
      energy = ph_temperature * planckSample()
    end if

    call generateRandomDirection(kx, ky, kz)
    call injectParticleGlobally(3, xg, yg, zg, energy * kx, energy * ky, energy * kz, 1.0)
  end subroutine injectPhoton

  subroutine userInitParticles()
    implicit none
    type(region) :: back_region
    integer :: s, ti, tj, tk, p
    real :: vel
    integer :: ncells, nphotons, n
    procedure(spatialDistribution), pointer :: spat_distr_ptr => null()
    spat_distr_ptr => userSpatialDistribution

    if (inj_photons) then
      ! injecting photons
      ncells = global_mesh % sx * global_mesh % sy * global_mesh % sz
      nphotons = INT(REAL(ppc0, 8) * REAL(ncells, 8))
      do n = 1, nphotons
        call injectPhoton(0)
      end do
    else
      ! injecting pairs
      back_region % x_min = 0.0
      back_region % x_max = global_mesh % sx
      back_region % y_min = 0.0
      back_region % y_max = global_mesh % sy
      call fillRegionWithThermalPlasma(back_region, (/1, 2/), 2, ppc0, T_init)
      
      if (inj_mono) then
        ! reset pair energies
        do s = 1, 2
          do ti = 1, species(s) % tile_nx
            do tj = 1, species(s) % tile_ny
              do tk = 1, species(s) % tile_nz
                do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
                  vel = sqrt(species(s) % prtl_tile(ti, tj, tk) % u(p)**2 + &
                            species(s) % prtl_tile(ti, tj, tk) % v(p)**2 + &
                            species(s) % prtl_tile(ti, tj, tk) % w(p)**2)
                  species(s) % prtl_tile(ti, tj, tk)%u(p) = T_init * species(s) % prtl_tile(ti, tj, tk)%u(p) / vel
                  species(s) % prtl_tile(ti, tj, tk)%v(p) = T_init * species(s) % prtl_tile(ti, tj, tk)%v(p) / vel
                  species(s) % prtl_tile(ti, tj, tk)%w(p) = T_init * species(s) % prtl_tile(ti, tj, tk)%w(p) / vel
                end do
              end do
            end do
          end do
        end do
      end if
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
  end subroutine userParticleBoundaryConditions

  subroutine userFieldBoundaryConditions(step, updateE, updateB)
    implicit none
    integer, optional, intent(in) :: step
    logical, optional, intent(in) :: updateE, updateB
  end subroutine userFieldBoundaryConditions
  !............................................................!

#include "optional.F"
end module m_userfile
