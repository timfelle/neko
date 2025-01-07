! Copyright (c) 2024, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!>
module mask
  use, intrinsic :: iso_c_binding, only : c_ptr
  use neko_config, only: NEKO_BCKND_DEVICE
  use device, only: device_map, device_memcpy, HOST_TO_DEVICE, device_free

  implicit none
  private

  !>
  type, public :: mask_t
     private

     !> The actual mask array.
     integer, dimension(:), allocatable :: array
     !> The device mask array.
     type(c_ptr) :: array_d
     !> The number of elements in the mask.
     integer :: mask_size

   contains
     !> Constructor for mask_t.
     procedure, public, pass(this) :: init => mask_init
     !> Destructor for mask_t.
     procedure, public, pass(this) :: free => mask_free

     !> Returns the size of the mask.
     procedure, public, pass(this) :: get_size => mask_get_size

     !> Returns the mask.
     generic, public :: get_host => mask_host_array, mask_host_value
     !> Returns the mask array.
     procedure, pass(this) :: mask_host_array
     !> Returns the mask value for a given index.
     procedure, pass(this) :: mask_host_value

     !> Returns the mask array on devices.
     generic :: get_device => mask_device_array
     !> Returns the mask array on devices.
     procedure, pass(this) :: mask_device_array

     !> Sets the mask value for a given index.
     procedure, pass(this) :: set_mask_value
     !> Assign the entire mask array to an array.
     procedure, pass(this) :: set_mask_array
     !> Assign the entire mask to another mask.
     procedure, pass(this) :: set_mask_mask

     ! Operator overloads
     generic :: assignment(=) => set_mask_array, set_mask_mask

  end type mask_t

contains

  !> Constructor
  subroutine mask_init(this, n, array)
    class(mask_t), intent(inout) :: this
    integer, intent(in) :: n
    integer, intent(in), dimension(n) :: array

    this%mask_size = n
    allocate(this%array(n))
    this%array = array

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%array, this%array_d, n)
       call device_memcpy(this%array, this%array_d, n, &
            HOST_TO_DEVICE, sync = .false.)
    end if

  end subroutine mask_init

  !> Destructor
  subroutine mask_free(this)
    class(mask_t), intent(inout) :: this

    if (NEKO_BCKND_DEVICE .eq. 1) call device_free(this%array_d)
    deallocate(this%array)

  end subroutine mask_free

  !> Returns the size of the mask.
  pure function mask_get_size(this) result(size)
    class(mask_t), intent(in) :: this
    integer :: size
    size = this%mask_size
  end function mask_get_size

  !> Returns the mask array.
  function mask_host_array(this) result(mask_array)
    class(mask_t), intent(in), target :: this
    integer, dimension(:), pointer :: mask_array
    mask_array => this%array
  end function mask_host_array

  !> Returns the mask array.
  pure function mask_host_value(this, i) result(mask_value)
    class(mask_t), intent(in) :: this
    integer, intent(in) :: i
    integer :: mask_value
    mask_value = this%array(i)
  end function mask_host_value

  !> Returns the mask array.
  function mask_device_array(this) result(mask_array)
    class(mask_t), intent(in), target :: this
    type(c_ptr) :: mask_array
    mask_array = this%array_d
  end function mask_device_array

  ! -------------------------------------------------------------------------- !

  !> Sets the mask value for a given index.
  subroutine set_mask_value(this, i, input)
    class(mask_t), intent(inout) :: this
    integer, intent(in) :: i
    integer, intent(in) :: input

    this%array(i) = input

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%array, this%array_d, this%mask_size, &
            HOST_TO_DEVICE, sync = .false.)
    end if

  end subroutine set_mask_value

  !> Assign the entire mask array to an array.
  subroutine set_mask_array(this, input)
    class(mask_t), intent(inout) :: this
    integer, intent(in), dimension(:) :: input

    call this%free()
    call this%init(size(input), input)

  end subroutine set_mask_array

  !> Assign the entire mask to another mask.
  subroutine set_mask_mask(this, other)
    class(mask_t), intent(inout) :: this
    class(mask_t), intent(in) :: other

    call this%free()
    call this%init(other%get_size(), other%get_host())

  end subroutine set_mask_mask

end module mask
