!> Shears the mesh in the user mesh_setup routine, and logs the number of
!! elements flagged as deformed once the case is initialised.
module user
  use neko
  implicit none

  !> Shear of the mesh, the displacement in y per unit x
  real(kind=dp) :: shear = 0.0_dp

contains

  !> Register user defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%startup => startup
    user%mesh_setup => mesh_setup
    user%initialize => initialize
  end subroutine user_setup

  !> Read the shear from the case file
  subroutine startup(params)
    type(json_file), intent(inout) :: params

    call json_get(params, 'case.mesh_shear', shear)
  end subroutine startup

  !> Shear the mesh, y in proportion to x
  subroutine mesh_setup(msh, time)
    type(mesh_t), intent(inout) :: msh
    type(time_state_t), intent(in) :: time
    integer :: i

    do i = 1, size(msh%points)
       msh%points(i)%x(2) = msh%points(i)%x(2) + shear * msh%points(i)%x(1)
    end do
  end subroutine mesh_setup

  !> Log the number of elements flagged as deformed
  subroutine initialize(time)
    type(time_state_t), intent(in) :: time
    type(mesh_t), pointer :: msh
    character(len=LOG_SIZE) :: log_buf
    integer :: n_deformed

    msh => neko_user_access%case%msh
    n_deformed = count(msh%dfrmd_el)
    call MPI_Allreduce(MPI_IN_PLACE, n_deformed, 1, MPI_INTEGER, MPI_SUM, &
         NEKO_COMM)

    write(log_buf, '(A, I0, A, I0)') 'Deformed elements: ', n_deformed, &
         ' of ', msh%glb_nelv
    call neko_log%message(log_buf)
  end subroutine initialize

end module user
