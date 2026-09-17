! Copyright (c) 2020-2026, The Neko Authors
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
!> Defines a Matrix-vector product
module ax_product
  use num_types, only : rp
  use coefs, only : coef_t
  use space, only : space_t
  use mesh, only : mesh_t
  use utils, only : neko_error
  use, intrinsic :: iso_c_binding, only : c_ptr
  implicit none
  private

  !> Base type for a matrix-vector product providing \f$ Ax \f$
  type, public, abstract :: ax_t
   contains
     procedure(ax_compute), pass(this), deferred :: compute
     procedure(ax_compute_vector), pass(this), deferred :: compute_vector
     procedure, pass(this) :: compute_device => ax_compute_device
     procedure, pass(this) :: compute_vector_device => ax_compute_vector_device
     procedure, pass(this) :: free => ax_free
  end type ax_t

  interface
     !> Allocate a Helmholtz problem matrix-vector product.
     !! The implementation is selected by name and compute backend.
     !! @param object The matrix-vector product type to be allocated.
     !! @param type_name The name of the matrix-vector product type.
     module subroutine ax_helm_allocator(object, type_name)
       class(ax_t), allocatable, intent(inout) :: object
       character(len=*), intent(in) :: type_name
     end subroutine ax_helm_allocator
  end interface

  !
  ! Machinery for injecting user-defined types
  !

  !> Interface for a matrix-vector product allocator.
  !! Implemented in user modules, it should allocate `obj` to the custom user
  !! type.
  abstract interface
     subroutine ax_helm_allocate(obj)
       import ax_t
       class(ax_t), allocatable, intent(inout) :: obj
     end subroutine ax_helm_allocate
  end interface

  interface
     !> Called in user modules to add an allocator for custom types.
     module subroutine register_ax_helm(type_name, allocator)
       character(len=*), intent(in) :: type_name
       procedure(ax_helm_allocate), pointer, intent(in) :: allocator
     end subroutine register_ax_helm
  end interface

  !> A name-allocator pair for user-defined matrix-vector product types.
  type ax_helm_allocator_entry
     character(len=20) :: type_name
     procedure(ax_helm_allocate), pointer, nopass :: allocator
  end type ax_helm_allocator_entry

  !> Registry of allocators for user-defined matrix-vector product types.
  type(ax_helm_allocator_entry), allocatable :: ax_helm_registry(:)

  !> The size of `ax_helm_registry`.
  integer :: ax_helm_registry_size = 0

  public :: ax_helm_allocator, register_ax_helm, ax_helm_allocate

  !> Abstract interface for computing\f$ Ax \f$ inside a Krylov method
  !!
  !! @param w Vector of size @a (lx,ly,lz,nelv).
  !! @param u Vector of size @a (lx,ly,lz,nelv).
  !! @param coef Coefficients.
  !! @param msh Mesh.
  !! @param Xh Function space \f$ X_h \f$.
  abstract interface
     subroutine ax_compute(this, w, u, coef, msh, Xh)
       import space_t
       import mesh_t
       import coef_t
       import ax_t
       import rp
       implicit none
       class(ax_t), intent(in) :: this
       type(space_t), intent(in) :: Xh
       type(mesh_t), intent(in) :: msh
       type(coef_t), intent(in) :: coef
       real(kind=rp), intent(inout) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
       real(kind=rp), intent(in) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
     end subroutine ax_compute
  end interface

  !> Abstract interface for computing\f$ Ax \f$ inside a Krylov method,
  !! taking 3 components of a vector field in a coupled manner.
  !! @param au Result for the first component of the vector.
  !! @param av Result for the second component of the vector.
  !! @param aw Result for the third component of the vector.
  !! @param u The first component of the vector.
  !! @param v The second component of the vector.
  !! @param w The third component of the vector.
  !! @param coef Coefficients.
  !! @param msh Mesh.
  !! @param Xh Function space \f$ X_h \f$.
  abstract interface
     subroutine ax_compute_vector(this, au, av, aw, u, v, w, coef, msh, Xh)
       import space_t
       import mesh_t
       import coef_t
       import ax_t
       import rp
       implicit none
       class(ax_t), intent(in) :: this
       type(space_t), intent(in) :: Xh
       type(mesh_t), intent(in) :: msh
       type(coef_t), intent(in) :: coef
       real(kind=rp), intent(inout) :: au(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
       real(kind=rp), intent(inout) :: av(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
       real(kind=rp), intent(inout) :: aw(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
       real(kind=rp), intent(in) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
       real(kind=rp), intent(in) :: v(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
       real(kind=rp), intent(in) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
     end subroutine ax_compute_vector
  end interface

contains

  !> Release resources held by a matrix-vector product.
  !!
  !! This default implementation permits resource-free operators to require no
  !! special handling. Descendants with pointer or allocatable components
  !! should override it.
  subroutine ax_free(this)
    class(ax_t), intent(inout) :: this
  end subroutine ax_free

  !> Default implementation of the device compute method.
  !! This default implementation just throws an error.
  subroutine ax_compute_device(this, w_d, u_d, coef, msh, Xh)
    class(ax_t), intent(in) :: this
    type(space_t), intent(in) :: Xh
    type(mesh_t), intent(in) :: msh
    type(coef_t), intent(in) :: coef
    type(c_ptr), intent(inout) :: w_d
    type(c_ptr), intent(in) :: u_d
    call neko_error("compute_device not implemented for this ax_t type")
  end subroutine ax_compute_device

  !> Default implementation of the device compute method.
  !! This default implementation just throws an error.
  subroutine ax_compute_vector_device(this, au_d, av_d, aw_d, u_d, v_d, w_d, &
       coef, msh, Xh)
    class(ax_t), intent(in) :: this
    type(space_t), intent(in) :: Xh
    type(mesh_t), intent(in) :: msh
    type(coef_t), intent(in) :: coef
    type(c_ptr), intent(inout) :: au_d
    type(c_ptr), intent(inout) :: av_d
    type(c_ptr), intent(inout) :: aw_d
    type(c_ptr), intent(in) :: u_d
    type(c_ptr), intent(in) :: v_d
    type(c_ptr), intent(in) :: w_d
    call neko_error("compute_vector_device not implemented for this ax_t type")
  end subroutine ax_compute_vector_device
end module ax_product
