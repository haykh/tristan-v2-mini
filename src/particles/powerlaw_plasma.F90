module m_powerlawplasma
  use m_globalnamespace
  use m_aux
  use m_helpers
  use m_errors
  use m_domain
  use m_particles
  use m_fields
  use m_particlelogistics
  implicit none

  !--- PRIVATE functions -----------------------------------------!
  !...............................................................!
contains

  subroutine fillRegionWithPowerlawPlasma(fill_region, fill_species, num_species, ndens_sp, &
                                          plaw_gmin, plaw_gmax, plaw_ind, &
                                          init_2dQ, weights, &
                                          spat_distr_ptr, &
                                          dummy1, dummy2, dummy3)
    implicit none
    ! assuming that the charges of all species given in `fill_species` add up to `0`
    type(region), intent(in) :: fill_region
    integer, intent(in) :: num_species
    integer, intent(in) :: fill_species(num_species)
    real, intent(in) :: ndens_sp, plaw_gmin, plaw_gmax, plaw_ind
    integer :: num_part, n, s, spec_
    integer(kind=2) :: xi_, yi_, zi_
    real :: fill_xmin, fill_xmax, &
            fill_ymin, fill_ymax, &
            fill_zmin, fill_zmax
    real :: u_, v_, w_, dx_, dy_, dz_
    real :: x_, y_, z_, gam_, bet_, TH, ZT, rnd
    real(kind=8) :: num_part_r
    real :: x_glob, y_glob, z_glob
    logical, optional :: init_2dQ
    logical :: init_2dQ_
    real, optional :: weights
    real :: weights_

    procedure(spatialDistribution), pointer, intent(in), optional :: spat_distr_ptr
    real, intent(in), optional :: dummy1, dummy2, dummy3
    real :: dummy1_, dummy2_, dummy3_

    if (present(dummy1)) then
      dummy1_ = dummy1
    else
      dummy1_ = 0.0
    end if
    if (present(dummy2)) then
      dummy2_ = dummy2
    else
      dummy2_ = 0.0
    end if
    if (present(dummy3)) then
      dummy3_ = dummy3
    else
      dummy3_ = 0.0
    end if

    if (.not. present(weights)) then
      weights_ = 1.0
    else
      weights_ = weights
    end if

    if (present(init_2dQ)) then
      init_2dQ_ = init_2dQ
    else
      init_2dQ_ = .false.
    end if

    ! global to local coordinates
#ifdef oneD
    call globalToLocalCoords(fill_region % x_min, 0.0, 0.0, &
                             fill_xmin, fill_ymin, fill_zmin, adjustQ=.true.)
    call globalToLocalCoords(fill_region % x_max, 0.0, 0.0, &
                             fill_xmax, fill_ymax, fill_zmax, adjustQ=.true.)
    num_part_r = REAL(ndens_sp, 8) * REAL(fill_xmax - fill_xmin, 8)
#elif defined(twoD)
    call globalToLocalCoords(fill_region % x_min, fill_region % y_min, 0.0, &
                             fill_xmin, fill_ymin, fill_zmin, adjustQ=.true.)
    call globalToLocalCoords(fill_region % x_max, fill_region % y_max, 0.0, &
                             fill_xmax, fill_ymax, fill_zmax, adjustQ=.true.)
    num_part_r = REAL(ndens_sp, 8) * REAL(fill_xmax - fill_xmin, 8) &
                 * REAL(fill_ymax - fill_ymin, 8)
#elif defined(threeD)
    call globalToLocalCoords(fill_region % x_min, fill_region % y_min, fill_region % z_min, &
                             fill_xmin, fill_ymin, fill_zmin, adjustQ=.true.)
    call globalToLocalCoords(fill_region % x_max, fill_region % y_max, fill_region % z_max, &
                             fill_xmax, fill_ymax, fill_zmax, adjustQ=.true.)
    num_part_r = REAL(ndens_sp, 8) * REAL(fill_xmax - fill_xmin, 8) &
                 * REAL(fill_ymax - fill_ymin, 8) &
                 * REAL(fill_zmax - fill_zmin, 8)
#endif

    num_part = INT(num_part_r) + 1

    do n = 1, num_part
      if (n .eq. num_part) then
        ! last particle
        num_part_r = num_part_r - REAL(num_part - 1, 8)
        if (random(dseed) .gt. num_part_r) then
          exit
        end if
      end if
      ! generate coords for all species
      call generateCoordInRegion(fill_xmin, fill_xmax, fill_ymin, fill_ymax, fill_zmin, fill_zmax, &
                                 x_, y_, z_, xi_, yi_, zi_, dx_, dy_, dz_)

      ! if spatial distribution function is present, compute it
      !   otherwise use uniform distribution
      if (present(spat_distr_ptr)) then
        x_glob = REAL(this_meshblock % ptr % x0) + x_
        y_glob = REAL(this_meshblock % ptr % y0) + y_
        z_glob = REAL(this_meshblock % ptr % z0) + z_
        rnd = spat_distr_ptr(x_glob=x_glob, y_glob=y_glob, z_glob=z_glob, &
                             dummy1=dummy1_, dummy2=dummy2_, dummy3=dummy3_)
      else
        rnd = 1.0
      end if
      if (random(dseed) .lt. rnd) then
        do s = 1, num_species
          ! generate momenta for every species individually
          spec_ = fill_species(s)
          if ((spec_ .le. 0) .or. (spec_ .gt. nspec)) then
            call throwError('Wrong species specified in fillRegionWithPowerlawPlasma.')
          end if
          !   generate powerlaw gamma:
          rnd = random(dseed)
          gam_ = ((plaw_gmax**(plaw_ind + 1) - plaw_gmin**(plaw_ind + 1)) * rnd + &
                  plaw_gmin**(plaw_ind + 1))**(1.0 / (plaw_ind + 1.0))
          bet_ = sqrt(1.0 - 1.0 / gam_**2)
          !   generate random direction
          if (init_2dQ_) then
            ! initialize momentum on 2d plane
            TH = random(dseed) * 2.0 * M_PI
            u_ = gam_ * bet_ * cos(TH)
            v_ = gam_ * bet_ * sin(TH)
            w_ = 0.0
          else
            ! initialize momentum in 3d
            TH = random(dseed) * 2.0 * M_PI; ZT = random(dseed) * 2.0 - 1.0
            u_ = gam_ * bet_ * sqrt(1.0 - ZT**2) * cos(TH)
            v_ = gam_ * bet_ * sqrt(1.0 - ZT**2) * sin(TH)
            w_ = gam_ * bet_ * ZT
          end if
          call createParticle(spec_, xi_, yi_, zi_, dx_, dy_, dz_, u_, v_, w_, weight=weights_)
        end do
      end if
    end do
  end subroutine fillRegionWithPowerlawPlasma
end module m_powerlawplasma
