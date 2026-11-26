!> A uniform state for the compressible solver: rho = 1, p = 1, u = (0.5, 0, 0).
module user
  use neko
  implicit none

contains

  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%initial_conditions => initial_conditions
  end subroutine user_setup

  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields

    type(field_t), pointer :: rho, u, v, w, p

    rho => fields%get_by_name("fluid_rho")
    u => fields%get_by_name("u")
    v => fields%get_by_name("v")
    w => fields%get_by_name("w")
    p => fields%get_by_name("p")

    rho%x = 1.0_rp
    p%x = 1.0_rp
    u%x = 0.5_rp
    v%x = 0.0_rp
    w%x = 0.0_rp

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call rho%copy_from(HOST_TO_DEVICE, sync = .false.)
       call u%copy_from(HOST_TO_DEVICE, sync = .false.)
       call v%copy_from(HOST_TO_DEVICE, sync = .false.)
       call w%copy_from(HOST_TO_DEVICE, sync = .false.)
       call p%copy_from(HOST_TO_DEVICE, sync = .true.)
    end if
  end subroutine initial_conditions

end module user
