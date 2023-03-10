module m_mover
  use m_globalnamespace
  use m_aux
  use m_helpers
  use m_domain
  use m_particles
  use m_fields
  use m_userfile
  use m_exchangearray, only: exchangeArray

  ! extra physics
#ifdef RADIATION
  use m_radiation
#endif

  implicit none
contains
  subroutine moveParticles(timestep)
    ! DEP_PRT [particle-dependent]
    implicit none
    integer, intent(in) :: timestep
    integer :: s, p, ti, tj, tk
    integer(kind=2) :: temp_i
    real :: g_temp, over_e_temp, temp_r
    integer(kind=2), pointer, contiguous :: pt_xi(:), pt_yi(:), pt_zi(:)
    real, pointer, contiguous :: pt_dx(:), pt_dy(:), pt_dz(:), &
                                 pt_u(:), pt_v(:), pt_w(:), pt_wei(:)
    real :: ex0, ey0, ez0, bx0, by0, bz0, q_over_m
    real :: u0, v0, w0, u1, v1, w1, dummy_, dx, dy, dz
#if defined (EXTERNALFIELDS)
    real :: ex_ext, ey_ext, ez_ext
    real :: bx_ext, by_ext, bz_ext
#endif
    real :: c000, c100, c001, c101, c010, c110, c011, c111, &
            c00, c01, c10, c11, c0, c1
    integer :: iy, iz, lind

#ifdef GCA
    ! variables with `_` at the end are temporary (or correspond to `t = n-1/2` or `t = n+1/2`) ...
    ! ... variables with `_n` at the end correspond to `t = n` ...
    ! ... variables with `_n1` at the end correspond to `t = n+1` ...
    integer(kind=2), pointer, contiguous :: pt_xi_past(:), pt_yi_past(:), pt_zi_past(:)
    real, pointer, contiguous :: pt_dx_past(:), pt_dy_past(:), pt_dz_past(:)
    real, pointer, contiguous :: pt_u_eff(:), pt_v_eff(:), pt_w_eff(:), pt_u_par(:), pt_u_perp(:)
    integer(kind=2) :: xi_, yi_, zi_
    real :: x_, y_, z_, dx_, dy_, dz_, Gamma_, mu_
    real :: x_n1, y_n1, z_n1, x_n, y_n, z_n
    real :: Gamma_n, Gamma_n1, e0_SQR, b0_SQR, ee0, bb0
    real :: ex0_n, ey0_n, ez0_n, bx0_n, by0_n, bz0_n, b0_SQR_n, bb0_n
    real :: ex0_n1, ey0_n1, ez0_n1, bx0_n1, by0_n1, bz0_n1, b0_SQR_n1, bb0_n1
    real :: u_par, epar0_n, bpar_x, bpar_y, bpar_z, bperp_x, bperp_y, bperp_z
    real :: uperp_g_x, uperp_g_y, uperp_g_z, uperp_g
    real :: vE_x, vE_y, vE_z, gammaE, wE_x, wE_y, wE_z, wE_SQR
    real :: vE_x_n, vE_y_n, vE_z_n, gammaE_n, wE_x_n, wE_y_n, wE_z_n, wE_SQR_n
    real :: vE_x_n1, vE_y_n1, vE_z_n1, gammaE_n1, wE_x_n1, wE_y_n1, wE_z_n1, wE_SQR_n1
    logical :: doNormalPushQ
    integer :: iter
#endif

#ifdef PRTLPAYLOADS
    real, pointer, contiguous :: pt_pld1(:), pt_pld2(:), pt_pld3(:)
    real :: incr_pld1, incr_pld2, incr_pld3
#endif

#if defined (RADIATION) || defined (GCA)
    integer, pointer, contiguous :: pt_proc(:)
#endif

#ifdef RADIATION
    logical :: dummy_flag
    integer(kind=2) :: xi_rad, yi_rad, zi_rad
    real :: ex_rad, ey_rad, ez_rad, bx_rad, by_rad, bz_rad
    real :: u_init, v_init, w_init, dx_rad, dy_rad, dz_rad
    integer, pointer, contiguous :: pt_ind(:)
#endif

    if (.false.) print *, timestep

    iy = this_meshblock % ptr % i2 - this_meshblock % ptr % i1 + 1
    iz = iy * (this_meshblock % ptr % j2 - this_meshblock % ptr % j1 + 1)

#ifdef RADIATION
    if (rad_dens_lim .gt. 0) then
      ! if density limit is enabled
      dummy_flag = .true.
      do s = 1, nspec
        if (species(s) % m_sp .ne. 0) then
          call computeDensity(s, reset=dummy_flag)
          dummy_flag = .false.
        end if
      end do
      call exchangeArray()
    end if
#endif

    do s = 1, nspec
      if (.not. species(s) % move_sp) cycle
      if (species(s) % m_sp .eq. 0) then
        do tk = 1, species(s) % tile_nz
          do tj = 1, species(s) % tile_ny
            do ti = 1, species(s) % tile_nx
              pt_xi => species(s) % prtl_tile(ti, tj, tk) % xi
              pt_yi => species(s) % prtl_tile(ti, tj, tk) % yi
              pt_zi => species(s) % prtl_tile(ti, tj, tk) % zi

              pt_dx => species(s) % prtl_tile(ti, tj, tk) % dx
              pt_dy => species(s) % prtl_tile(ti, tj, tk) % dy
              pt_dz => species(s) % prtl_tile(ti, tj, tk) % dz

              pt_u => species(s) % prtl_tile(ti, tj, tk) % u
              pt_v => species(s) % prtl_tile(ti, tj, tk) % v
              pt_w => species(s) % prtl_tile(ti, tj, tk) % w

#ifdef PRTLPAYLOADS
              pt_pld1 => species(s) % prtl_tile(ti, tj, tk) % payload1
              pt_pld2 => species(s) % prtl_tile(ti, tj, tk) % payload2
              pt_pld3 => species(s) % prtl_tile(ti, tj, tk) % payload3
#endif

              ! routine for massless particles
#ifdef PRTLPAYLOADS
              !$omp simd private(over_e_temp, temp_r, temp_i, incr_pld1, incr_pld2, incr_pld3, u0, v0, w0)
#else
              !$omp simd private(over_e_temp, temp_r, temp_i)
#endif
              !dir$ vector aligned
              do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
                ! move particle
                over_e_temp = 1.0 / sqrt(pt_u(p)**2 + pt_v(p)**2 + pt_w(p)**2)
#ifdef PRTLPAYLOADS
                u0 = pt_u(p); v0 = pt_v(p); w0 = pt_w(p)
                call usrSetPhPld(u0, v0, w0, over_e_temp, incr_pld1, incr_pld2, incr_pld3)
                pt_pld1(p) = pt_pld1(p) + incr_pld1
                pt_pld2(p) = pt_pld2(p) + incr_pld2
                pt_pld3(p) = pt_pld3(p) + incr_pld3
#endif
                ! this "function" takes
                ! ... inverse energy: `over_e_temp` ...
                ! ... reads the velocities from: `pt_*(p)` ...
                ! ... and updates the particle position `pt_*(p)`
#include "position_update.F08"
              end do ! p
              pt_xi => null(); pt_yi => null(); pt_zi => null()
              pt_dx => null(); pt_dy => null(); pt_dz => null()
              pt_u => null(); pt_v => null(); pt_w => null()
#ifdef PRTLPAYLOADS
              pt_pld1 => null(); pt_pld2 => null(); pt_pld3 => null()
#endif
            end do ! ti
          end do ! tj
        end do ! tk
      else ! massive particles
        do tk = 1, species(s) % tile_nz
          do tj = 1, species(s) % tile_ny
            do ti = 1, species(s) % tile_nx
              pt_xi => species(s) % prtl_tile(ti, tj, tk) % xi
              pt_yi => species(s) % prtl_tile(ti, tj, tk) % yi
              pt_zi => species(s) % prtl_tile(ti, tj, tk) % zi

              pt_dx => species(s) % prtl_tile(ti, tj, tk) % dx
              pt_dy => species(s) % prtl_tile(ti, tj, tk) % dy
              pt_dz => species(s) % prtl_tile(ti, tj, tk) % dz

              pt_u => species(s) % prtl_tile(ti, tj, tk) % u
              pt_v => species(s) % prtl_tile(ti, tj, tk) % v
              pt_w => species(s) % prtl_tile(ti, tj, tk) % w

              pt_wei => species(s) % prtl_tile(ti, tj, tk) % weight

#ifdef PRTLPAYLOADS
              pt_pld1 => species(s) % prtl_tile(ti, tj, tk) % payload1
              pt_pld2 => species(s) % prtl_tile(ti, tj, tk) % payload2
              pt_pld3 => species(s) % prtl_tile(ti, tj, tk) % payload3
#endif

#ifdef GCA
              pt_xi_past => species(s) % prtl_tile(ti, tj, tk) % xi_past
              pt_yi_past => species(s) % prtl_tile(ti, tj, tk) % yi_past
              pt_zi_past => species(s) % prtl_tile(ti, tj, tk) % zi_past

              pt_dx_past => species(s) % prtl_tile(ti, tj, tk) % dx_past
              pt_dy_past => species(s) % prtl_tile(ti, tj, tk) % dy_past
              pt_dz_past => species(s) % prtl_tile(ti, tj, tk) % dz_past

              pt_u_eff => species(s) % prtl_tile(ti, tj, tk) % u_eff
              pt_v_eff => species(s) % prtl_tile(ti, tj, tk) % v_eff
              pt_w_eff => species(s) % prtl_tile(ti, tj, tk) % w_eff

              pt_u_par => species(s) % prtl_tile(ti, tj, tk) % u_par
              pt_u_perp => species(s) % prtl_tile(ti, tj, tk) % u_perp
#endif

#if defined(RADIATION) || defined(GCA)
              pt_proc => species(s) % prtl_tile(ti, tj, tk) % proc
#endif

#if defined(RADIATION)
              pt_ind => species(s) % prtl_tile(ti, tj, tk) % ind
#endif

              ! routine for massive particles
              q_over_m = species(s) % ch_sp / species(s) % m_sp
#if !defined(RADIATION) && !defined(EXTERNALFIELDS) && !defined(GCA)
              !$omp simd private(lind, dummy_, g_temp, over_e_temp,&
              !$omp  temp_r, temp_i, u0, v0, w0, u1, v1, w1,&
              !$omp  ex0, ey0, ez0, bx0, by0, bz0,&
#ifdef PRTLPAYLOADS
              !$omp  incr_pld1, incr_pld2, incr_pld3,&
#endif
              !$omp  c000, c100, c001, c101, c010, c110, c011, c111,&
              !$omp  c00, c01, c10, c11, c0, c1)
              !dir$ vector aligned
#endif
#ifdef GCA
              !$omp simd private(lind, dummy_, g_temp, over_e_temp,&
              !$omp  temp_r, temp_i, u0, v0, w0, u1, v1, w1,&
              !$omp  ex0, ey0, ez0, bx0, by0, bz0,&
              !$omp  c000, c100, c001, c101, c010, c110, c011, c111,&
              !$omp  c00, c01, c10, c11, c0, c1,&
              !$omp  xi_, yi_, zi_, x_, y_, z_, dx_, dy_, dz_, Gamma_, mu_,&
              !$omp  x_n1, y_n1, z_n1, x_n, y_n, z_n,&
              !$omp  Gamma_n, Gamma_n1, e0_SQR, b0_SQR, ee0, bb0,&
              !$omp  ex0_n, ey0_n, ez0_n, bx0_n, by0_n, bz0_n, b0_SQR_n, bb0_n,&
              !$omp  ex0_n1, ey0_n1, ez0_n1, bx0_n1, by0_n1, bz0_n1, b0_SQR_n1, bb0_n1,&
              !$omp  u_par, epar0_n, bpar_x, bpar_y, bpar_z, bperp_x, bperp_y, bperp_z,&
              !$omp  uperp_g_x, uperp_g_y, uperp_g_z, uperp_g,&
              !$omp  vE_x, vE_y, vE_z, gammaE, wE_x, wE_y, wE_z, wE_SQR,&
              !$omp  vE_x_n, vE_y_n, vE_z_n, gammaE_n, wE_x_n, wE_y_n, wE_z_n, wE_SQR_n,&
              !$omp  vE_x_n1, vE_y_n1, vE_z_n1, gammaE_n1, wE_x_n1, wE_y_n1, wE_z_n1, wE_SQR_n1,&
              !$omp  iter, doNormalPushQ)
              !dir$ vector aligned
#endif
              do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
#ifndef DEBUG
                ! use "fast" interpolation
#ifdef oneD
                lind = pt_xi(p)
#elif defined(twoD)
                lind = pt_xi(p) + (NGHOST + pt_yi(p)) * iy
#elif defined(threeD)
                lind = pt_xi(p) + (NGHOST + pt_yi(p)) * iy + (NGHOST + pt_zi(p)) * iz
#endif
                dx = pt_dx(p); dy = pt_dy(p); dz = pt_dz(p)
                ! these "functions" take
                ! ... coordinates and linear index: `dx`, `dy`, `dz` and `lind` ...
                ! ... and "return" `bx0`, `by0`, `bz0`, `ex0`, `ey0`, `ez0`
#include "interp_efield.F08"
#include "interp_bfield.F08"
#else
                dx = pt_dx(p); dy = pt_dy(p); dz = pt_dz(p)
                call interpFromEdges(dx, dy, dz, pt_xi(p), pt_yi(p), pt_zi(p), &
                                     ex, ey, ez, ex0, ey0, ez0)
                call interpFromFaces(dx, dy, dz, pt_xi(p), pt_yi(p), pt_zi(p), &
                                     bx, by, bz, bx0, by0, bz0)
#endif

#if defined (EXTERNALFIELDS) && !defined (GCA)
                call userExternalFields(REAL(pt_xi(p)) + pt_dx(p), &
                                        REAL(pt_yi(p)) + pt_dy(p), &
                                        REAL(pt_zi(p)) + pt_dz(p), &
                                        ex_ext, ey_ext, ez_ext, &
                                        bx_ext, by_ext, bz_ext)
                ex0 = ex0 + ex_ext; ey0 = ey0 + ey_ext; ez0 = ez0 + ez_ext
                bx0 = bx0 + bx_ext; by0 = by0 + by_ext; bz0 = bz0 + bz_ext
#endif

#ifdef RADIATION
                ! save fields at time `t = n`
                ex_rad = ex0; ey_rad = ey0; ez_rad = ez0
                bx_rad = bx0; by_rad = by0; bz_rad = bz0

                ! save velocities before the push
                u_init = pt_u(p); v_init = pt_v(p); w_init = pt_w(p)

                ! save coordinates before the push
                dx_rad = pt_dx(p); dy_rad = pt_dy(p); dz_rad = pt_dz(p)
                xi_rad = pt_xi(p); yi_rad = pt_yi(p); zi_rad = pt_zi(p)
#endif

#ifndef GCA
                ! . . . . simple Boris/Vay pusher . . . .
                u0 = pt_u(p); v0 = pt_v(p); w0 = pt_w(p)
                ! this "function" takes
                ! ... the field quantities: `bx0`, `by0`, `bz0`, `ex0`, `ey0`, `ez0` ...
                ! ... and the velocities: `u0`, `v0`, `w0` ...
                ! ... and returns the updated velocities `u0`, `v0`, `w0`
#ifndef VAY
#include "boris_push.F08"
#else
#include "vay_push.F08"
#endif
                pt_u(p) = u0; pt_v(p) = v0; pt_w(p) = w0
                over_e_temp = 1.0 / sqrt(1.0 + pt_u(p)**2 + pt_v(p)**2 + pt_w(p)**2)
#ifdef PRTLPAYLOADS
                call usrSetElPld(q_over_m, u0, v0, w0, over_e_temp, ex0, ey0, ez0, bx0, by0, bz0, incr_pld1, incr_pld2, incr_pld3)
                pt_pld1(p) = pt_pld1(p) + incr_pld1
                pt_pld2(p) = pt_pld2(p) + incr_pld2
                pt_pld3(p) = pt_pld3(p) + incr_pld3
#endif
                ! this "function" takes
                ! ... inverse energy: `over_e_temp` ...
                ! ... reads the velocities from: `pt_*(p)` ...
                ! ... and updates the particle position `pt_*(p)`
#include "position_update.F08"
#else
                ! . . . . hybrid Boris/GCA pusher . . . .
                ! this "function"
                ! ... reads velocities: `pt_*(p)`
                ! ... and coordinates: `pt_*(p)`
                ! ... then updates the particle position, velocity, and past coordinates
#include "gca_routine.F08"
#endif

                ! RADIATION >
#ifdef RADIATION
#ifdef SYNCHROTRON
                if (species(s) % cool_sp) then
                  call particleRadiateSync(timestep, s, &
                                           pt_u(p), pt_v(p), pt_w(p), u_init, v_init, w_init, &
                                           dx_rad, dy_rad, dz_rad, xi_rad, yi_rad, zi_rad, pt_wei(p), &
                                           bx_rad, by_rad, bz_rad, ex_rad, ey_rad, ez_rad, &
                                           index=pt_ind(p), proc=pt_proc(p))
                end if
#endif
#ifdef INVERSECOMPTON
                if (species(s) % cool_sp) then
                  call particleRadiateIC(timestep, s, &
                                         pt_u(p), pt_v(p), pt_w(p), u_init, v_init, w_init, &
                                         dx_rad, dy_rad, dz_rad, xi_rad, yi_rad, zi_rad, pt_wei(p), &
                                         bx_rad, by_rad, bz_rad, ex_rad, ey_rad, ez_rad, &
                                         index=pt_ind(p), proc=pt_proc(p))
                end if
#endif
#endif
                ! </ RADIATION
              end do
              pt_xi => null(); pt_yi => null(); pt_zi => null()
              pt_dx => null(); pt_dy => null(); pt_dz => null()
              pt_u => null(); pt_v => null(); pt_w => null()

              pt_wei => null()

#ifdef PRTLPAYLOADS
              pt_pld1 => null(); pt_pld2 => null(); pt_pld3 => null()
#endif

#ifdef GCA
              pt_xi_past => null(); pt_yi_past => null(); pt_zi_past => null()
              pt_dx_past => null(); pt_dy_past => null(); pt_dz_past => null()
              pt_u_eff => null(); pt_v_eff => null(); pt_w_eff => null()

              pt_u_par => null(); pt_u_perp => null()
#endif

#if defined (RADIATION) || defined (GCA)
              pt_proc => null(); 
#endif

#if defined(RADIATION)
              pt_ind => null()
#endif

            end do ! ti
          end do ! tj
        end do ! tk
      end if
    end do ! species
    call printDiag("moveParticles()", 2)

  end subroutine moveParticles
end module m_mover

