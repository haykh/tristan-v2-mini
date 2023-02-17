! Configuration for this userfile:
! ```
!   $ python configure.py -2d --user=user_cooling --radiation=sync -payload
! ```

module m_userfile
  use m_globalnamespace
  use m_aux
  use m_helpers
  use m_readinput
  use m_domain
  use m_particles
  use m_fields
  use m_powerlawplasma
  use m_particlelogistics
  implicit none

  !--- PRIVATE variables -----------------------------------------!
  real, private :: plaw_gmin, plaw_gmax, plaw_ind, t_esc
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userSpatialDistribution
  !...............................................................!
contains
  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
    call getInput('problem', 'plaw_gmin', plaw_gmin)
    call getInput('problem', 'plaw_gmax', plaw_gmax)
    call getInput('problem', 'plaw_ind', plaw_ind)
    call getInput('problem', 't_esc', t_esc)
    if (plaw_gmin .lt. 1.0) then
      call throwError("plaw_gmin must be >= 1.0")
    else if (plaw_gmin .ge. plaw_gmax) then
      call throwError("plaw_gmin must be < plaw_gmax")
    end if
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

  subroutine userInitParticles()
    implicit none
    procedure(spatialDistribution), pointer :: spat_distr_ptr => null()
    spat_distr_ptr => userSpatialDistribution
  end subroutine userInitParticles

  subroutine userInitFields()
    implicit none
    integer :: i, j, k
    integer :: i_glob, j_glob, k_glob
    ! initialize uniform B-field in z
    ex(:, :, :) = 0; ey(:, :, :) = 0; ez(:, :, :) = 0
    bx(:, :, :) = 0; by(:, :, :) = 0; bz(:, :, :) = 1.0
    jx(:, :, :) = 0; jy(:, :, :) = 0; jz(:, :, :) = 0
  end subroutine userInitFields
  !............................................................!

  !--- driving ------------------------------------------------!
  subroutine userCurrentDeposit(step)
    implicit none
    integer, optional, intent(in) :: step
    ! called after particles move and deposit ...
    ! ... and before the currents are added to the electric field
  end subroutine userCurrentDeposit

  subroutine userDriveParticles(step)
    implicit none
    integer, optional, intent(in) :: step
    type(region) :: back_region
    integer, dimension(1000) :: nbefore
    integer :: s, ti, tj, tk, p, cntr
    nbefore(:) = -1

    ! record last particle indices in each tile before injection
    cntr = 0
    do s = 1, 2
      do ti = 1, species(s)%tile_nx
        do tj = 1, species(s)%tile_ny
          do tk = 1, species(s)%tile_nz
            nbefore(cntr) = species(s)%prtl_tile(ti, tj, tk)%npart_sp
            cntr = cntr + 1
            if (cntr .gt. 1000) then
              call throwError("ERROR: nbefore array too small")
            end if
          end do
        end do
      end do
    end do

    back_region % x_min = 0.0
    back_region % x_max = global_mesh % sx
    back_region % y_min = 0.0
    back_region % y_max = global_mesh % sy

    call fillRegionWithPowerlawPlasma(back_region, (/1, 2/), 2, REAL(ppc0 / t_esc), &
                                      plaw_gmin, plaw_gmax, plaw_ind, .true.)
    !                                                                   ^
    !                                                                   |
    !                                                       this flag ensures particle
    !                                                       velocities are in XY plane
    
    ! assign injection time to new particles
    cntr = 0
    do s = 1, 2
      do ti = 1, species(s)%tile_nx
        do tj = 1, species(s)%tile_ny
          do tk = 1, species(s)%tile_nz

            if (nbefore(cntr) .lt. 0) then
              call throwError("ERROR: invalid operation in userDriveParticles")
            end if

            do p = nbefore(cntr) + 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
              ! record the injection time for the new particles
              species(s) % prtl_tile(ti, tj, tk) % payload1(p) = step
            end do
            cntr = cntr + 1
          end do
        end do
      end do
    end do

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
    integer :: s, ti, tj, tk, p
    do s = 1, 2
      do ti = 1, species(s)%tile_nx
        do tj = 1, species(s)%tile_ny
          do tk = 1, species(s)%tile_nz
            do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
              ! schedule particle for deletion once its age is > t_esc
              if (step - species(s) % prtl_tile(ti, tj, tk) % payload1(p) .gt. t_esc) then
                species(s) % prtl_tile(ti, tj, tk) % proc(p) = -1
              end if
            end do
          end do
        end do
      end do
    end do
  end subroutine userParticleBoundaryConditions

  subroutine userFieldBoundaryConditions(step, updateE, updateB)
    implicit none
    integer, optional, intent(in) :: step
    logical, optional, intent(in) :: updateE, updateB
  end subroutine userFieldBoundaryConditions
  !............................................................!

#include "optional.F"
end module m_userfile
