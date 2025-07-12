module m_userfile
#ifdef HDF5
  use hdf5
  use m_globalnamespace, only: h5comm, h5info
#endif
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
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: ReadRealScalar, userSpatialDistribution
  !...............................................................!
contains
  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
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

  subroutine ReadRealScalar(filename, datasetname, buffer)
    character(len=*), intent(in) :: filename 
    character(len=*), intent(in) :: datasetname
    real, intent(out) :: buffer
    character(len=1), parameter :: groupname = "/"
    integer(HSIZE_T), allocatable, dimension(:) :: dims
    integer(HID_T) :: file_id, group_id, dset_id, dspace_id
    integer :: error, status
    integer :: h5_error_id

    h5_error_id = -1

    call h5open_f(error)
    if (error .ne. h5_error_id) then
      call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, error)
      if (error .ne. h5_error_id) then
        call h5gopen_f(file_id, groupname, group_id, error)
        if (error .ne. h5_error_id) then
          call h5dopen_f(group_id, trim(datasetname), dset_id, error)
          if (error .ne. h5_error_id) then
            call h5dget_space_f(dset_id, dspace_id, error)
            if (error .ne. h5_error_id) then
              call h5dread_f(dset_id, H5T_NATIVE_REAL, buffer, dims, error)
              if (error .ne. h5_error_id) then
                print *, 'successfully read dataset value: ', buffer
              else
                print *, 'error: failed to read dataset value.'
              end if
              call h5sclose_f(dspace_id, error)
            else
              print *, 'error: could not get dataspace.'
            end if
            call h5dclose_f(dset_id, error)
          else
            print *, 'error: could not open dataset ', trim(datasetname)
          end if
          call h5gclose_f(group_id, error)
        else
          print *, 'error: could not open group ', trim(groupname)
        end if
        call h5fclose_f(file_id, error)
      else
        print *, 'error: could not open file ', trim(filename)
      end if
      call h5close_f(error)
    else
      print *, 'error: hdf5 interface initialization failed.'
    end if
  end subroutine ReadRealScalar

  subroutine ReadRealArray(filename, datasetname, buffer)
    character(len=*), intent(in) :: filename 
    character(len=*), intent(in) :: datasetname
    real, allocatable, dimension(:,:,:), intent(out) :: buffer
    character(len=1), parameter :: groupname = "/"
    integer(HSIZE_T), dimension(3) :: dims, maxdims
    integer(HID_T) :: file_id, group_id, dset_id, dspace_id
    integer :: error, status
    integer :: h5_error_id

    h5_error_id = -1

    call h5open_f(error)
    if (error .ne. h5_error_id) then
      call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, error)
      if (error .ne. h5_error_id) then
        call h5gopen_f(file_id, groupname, group_id, error)
        if (error .ne. h5_error_id) then
          call h5dopen_f(group_id, trim(datasetname), dset_id, error)
          if (error .ne. h5_error_id) then
            call h5dget_space_f(dset_id, dspace_id, error)
            if (error .ne. h5_error_id) then
              call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, error)
              if (error .ne. h5_error_id) then
                print *, 'Dataset dimensions: ', dims(1), dims(2), dims(3)
                allocate (buffer(dims(1), dims(2), dims(3)))
                call h5dread_f(dset_id, H5T_NATIVE_REAL, buffer, dims, error)
                if (error .ne. h5_error_id) then
                  print *, 'successfully read 3d dataset.'
                else
                  print *, 'error: failed to read dataset.'
                end if
              else
                print *, 'error: could not get dataspace dimensions.'
              end if
              call h5sclose_f(dspace_id, error)
            else
              print *, 'error: could not get dataspace.'
            end if
            call h5dclose_f(dset_id, error)
          else
            print *, 'error: could not open dataset ', trim(datasetname)
          end if
          call h5gclose_f(group_id, error)
        else
          print *, 'error: could not open group ', trim(groupname)
        end if
        call h5fclose_f(file_id, error)
      else
        print *, 'error: could not open file ', trim(filename)
      end if
      call h5close_f(error)
    else
      print *, 'error: hdf5 interface initialization failed.'
    end if
  end subroutine ReadRealArray

  subroutine userInitFields()
    implicit none
    real :: c_value, mx0_value
    real, allocatable, dimension(:,:,:) :: bx_G, by_G, bz_G, ex_G, ey_G, ez_G
    integer :: xmin, xmax, ymin, ymax, zmin, zmax
    integer :: Gxmin, Gxmax, Gymin, Gymax, Gzmin, Gzmax

    character(len=STR_MAX), parameter :: params_filename = "/mnt/home/vrohoza/shared-valeriia/runs/ic_v2_2d_cf32/sig10.comp5.ppc4.bz1e-1/output/params.00000"
    character(len=STR_MAX), parameter :: fields_filename = "/mnt/home/vrohoza/shared-valeriia/runs/ic_v2_2d_cf32/sig10.comp5.ppc4.bz1e-1/output/flds/flds.tot.00000"
    call ReadRealScalar(params_filename, "algorithm:c", c_value)
    call ReadRealScalar(params_filename, "grid:mx0", mx0_value)
    call ReadRealArray(fields_filename, "bx", bx_G)
    call ReadRealArray(fields_filename, "by", by_G)
    call ReadRealArray(fields_filename, "bz", bz_G)
    call ReadRealArray(fields_filename, "ex", ex_G)
    call ReadRealArray(fields_filename, "ey", ey_G)
    call ReadRealArray(fields_filename, "ez", ez_G)

    xmin = 0; xmax = this_meshblock % ptr % sx
    ymin = 0; ymax = this_meshblock % ptr % sy
    zmin = 0; zmax = 1

    Gxmin = this_meshblock % ptr % x0 + 1; Gxmax = Gxmin + this_meshblock % ptr % sx
    Gymin = this_meshblock % ptr % y0 + 1; Gymax = Gymin + this_meshblock % ptr % sy
    Gzmin = 0; Gzmax = 1

    ex(xmin : xmax, ymin : ymax, zmin : zmax) = ex_G(Gxmin : Gxmax, Gymin : Gymax, Gzmin : Gzmax)
    ey(xmin : xmax, ymin : ymax, zmin : zmax) = ey_G(Gxmin : Gxmax, Gymin : Gymax, Gzmin : Gzmax)
    ez(xmin : xmax, ymin : ymax, zmin : zmax) = ez_G(Gxmin : Gxmax, Gymin : Gymax, Gzmin : Gzmax)
    by(xmin : xmax, ymin : ymax, zmin : zmax) = bx_G(Gxmin : Gxmax, Gymin : Gymax, Gzmin : Gzmax)
    by(xmin : xmax, ymin : ymax, zmin : zmax) = by_G(Gxmin : Gxmax, Gymin : Gymax, Gzmin : Gzmax)
    bz(xmin : xmax, ymin : ymax, zmin : zmax) = bz_G(Gxmin : Gxmax, Gymin : Gymax, Gzmin : Gzmax)
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
