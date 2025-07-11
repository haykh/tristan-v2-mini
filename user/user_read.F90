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
  private :: userSpatialDistribution
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

  subroutine userInitFields()
    implicit none
    CHARACTER(LEN=12), PARAMETER :: filename = "your_file.h5"
    CHARACTER(LEN=1), PARAMETER :: groupname = "/"
    CHARACTER(LEN=2), PARAMETER :: datasetname = "bx"

    ! HDF5 identifiers
    INTEGER(HID_T) :: file_id, group_id, dset_id, dspace_id
    INTEGER :: error, status

    ! Data and dimension variables
    INTEGER(HSIZE_T), DIMENSION(3) :: dims
    REAL, ALLOCATABLE, DIMENSION(:, :, :) :: data_out
    ! initialize uniform B-field in z

    ! 1. Initialize the HDF5 Fortran interface
    CALL h5open_f(h5_status)

    IF (h5_status .eq. 0) THEN
      ! 2. Open the HDF5 file for read-only access
      CALL h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)

      IF (error .eq. 0) THEN
        ! 3. Open the group
        CALL h5gopen_f(file_id, groupname, group_id, error)

        IF (error .eq. 0) THEN
          ! 4. Open the dataset
          CALL h5dopen_f(group_id, datasetname, dset_id, error)

          IF (error .eq. 0) THEN
            ! 5. Get the dataspace from the dataset
            CALL h5dget_space_f(dset_id, dspace_id, error)

            IF (error .eq. 0) THEN
              ! 6. Get the dimensions of the dataspace
              CALL h5sget_simple_extent_dims_f(dspace_id, dims, NULL, error)

              IF (error .eq. 0) THEN
                ! Fortran reads arrays in column-major order, while HDF5 uses row-major.
                ! The Fortran HDF5 wrapper handles this transposition automatically.
                PRINT *, 'Dataset dimensions: ', dims(1), dims(2), dims(3)

                ! 7. Allocate memory for the data array
                ALLOCATE (data_out(dims(1), dims(2), dims(3)), STAT=error)

                IF (error .eq. 0) THEN
                  ! 8. Read the data from the dataset
                  CALL h5dread_f(dset_id, H5T_NATIVE_REAL, data_out, dims, error)

                  IF (error .eq. 0) THEN
                    PRINT *, 'Successfully read 3D dataset "bx".'
                    !
                    ! You can now work with the data in the 'data_out' array
                    !
                  ELSE
                    PRINT *, 'Error: Failed to read dataset.'
                  END IF

                  ! Deallocate the array when done
                  DEALLOCATE (data_out)
                ELSE
                  PRINT *, 'Error: Could not allocate memory for data array.'
                END IF
              ELSE
                PRINT *, 'Error: Could not get dataspace dimensions.'
              END IF
              ! Close the dataspace
              CALL h5sclose_f(dspace_id, error)
            ELSE
              PRINT *, 'Error: Could not get dataspace.'
            END IF
            ! Close the dataset
            CALL h5dclose_f(dset_id, error)
          ELSE
            PRINT *, 'Error: Could not open dataset ', TRIM(datasetname)
          END IF
          ! Close the group
          CALL h5gclose_f(group_id, error)
        ELSE
          PRINT *, 'Error: Could not open group ', TRIM(groupname)
        END IF
        ! Close the file
        CALL h5fclose_f(file_id, error)
      ELSE
        PRINT *, 'Error: Could not open file ', TRIM(filename)
      END IF

      ! 10. Terminate the HDF5 Fortran interface
      CALL h5close_f(h5_status)
    ELSE
      PRINT *, 'Error: HDF5 interface initialization failed.'
    END IF
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
