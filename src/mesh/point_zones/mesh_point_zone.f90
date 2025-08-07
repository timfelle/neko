! Copyright (c) 2019-2021, The Neko Authors
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
! Implements a mesh based point zone
module mesh_point_zone
  use point_zone, only : point_zone_t
  use num_types, only : rp
  use json_utils, only : json_get, json_get_or_default
  use json_module, only : json_file
  use math, only : abscmp
  use tri_mesh, only : tri_mesh_t
  use aabb_tree, only : aabb_tree_t
  use signed_distance_module, only : signed_distance
  use file, only : file_t
  implicit none
  private

  !> A point zone defined by a mesh
  !! @details As defined here, a sphere is described by its center of
  !! coordinates `x0,y0,z0` and its radius, specified in the json file
  !! as e.g. `"center": [<x0>, <y0>, <z0>]", "radius": <r>`.
  type, public, extends(point_zone_t) :: mesh_point_zone_t

     type(tri_mesh_t) :: mesh
     type(aabb_tree_t) :: search_tree
     real(kind=rp) :: threshold

   contains
     !> Constructor from json object file.
     procedure, pass(this) :: init => mesh_point_zone_init_from_json
     !> Destructor.
     procedure, pass(this) :: free => mesh_point_zone_free
     !> Defines the criterion of selection of a GLL point in the sphere point
     !! zone.
     procedure, pass(this) :: criterion => mesh_point_zone_criterion
  end type mesh_point_zone_t

contains

  !> Constructor from json object file.
  !! @param json Json object file.
  !! @param size Size with which to initialize the stack
  subroutine mesh_point_zone_init_from_json(this, json, size)
    class(mesh_point_zone_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    integer, intent(in) :: size

    character(len=:), allocatable :: name
    character(len=:), allocatable :: filename
    real(kind=rp) :: threshold
    logical :: invert

    call json_get(json, "name", name)
    call json_get(json, "filename", filename)
    call json_get_or_default(json, "invert", invert, .false.)
    call json_get_or_default(json, "threshold", threshold, 0.0_rp)


    call mesh_point_zone_init_common(this, size, trim(name), invert, &
         filename, threshold)

  end subroutine mesh_point_zone_init_from_json

  !> Initializes a sphere point zone from its center coordinates and radius.
  !! @param size Size of the scratch stack.
  !! @param name Name of the sphere point zone.
  !! @param x0 Sphere center's x-coordinate.
  !! @param y0 Sphere center's y-coordinate.
  !! @param z0 Sphere center's z-coordinate.
  !! @param radius Sphere radius.
  subroutine mesh_point_zone_init_common(this, size, name, invert, &
       filename, threshold)
    class(mesh_point_zone_t), intent(inout) :: this
    integer, intent(in), optional :: size
    character(len=*), intent(in) :: name
    logical, intent(in) :: invert
    character(len=*), intent(in) :: filename
    real(kind=rp), intent(in) :: threshold
    type(file_t) :: mesh_file

    call this%init_base(size, name, invert)

    call mesh_file%init(filename)
    call mesh_file%read(this%mesh)

    call this%search_tree%init(this%mesh%nelv)
    call this%search_tree%build(this%mesh%el)

  end subroutine mesh_point_zone_init_common

  !> Destructor.
  subroutine mesh_point_zone_free(this)
    class(mesh_point_zone_t), intent(inout) :: this

    call this%free_base()
    call this%mesh%free()
  end subroutine mesh_point_zone_free

  !> Defines the criterion of selection of a GLL point in the sphere point zone.
  !! A GLL point of coordinates \f$ \vec{X} = (x, y, z) \f$ is considered as
  !! being inside the zone if:
  !! \f{eqnarray*}{
  !!    |\vec{X} - \vec{X_0}|^2 \le r
  !! \f}
  !! Where \f$ r \f$ is the radius of the sphere and
  !! \f$ \vec{X_0} = (x_0, y_0, z_0) \f$
  !! the coordinates of its center.
  !! @param x x-coordinate of the GLL point.
  !! @param y y-coordinate of the GLL point.
  !! @param z z-coordinate of the GLL point.
  !! @param j 1st nonlinear index of the GLL point.
  !! @param k 2nd nonlinear index of the GLL point.
  !! @param l 3rd nonlinear index of the GLL point.
  !! @param e element index of the GLL point.
  pure function mesh_point_zone_criterion(this, x, y, z, j, k, l, e) &
       result(is_inside)
    class(mesh_point_zone_t), intent(in) :: this
    real(kind=rp), intent(in) :: x
    real(kind=rp), intent(in) :: y
    real(kind=rp), intent(in) :: z
    integer, intent(in) :: j
    integer, intent(in) :: k
    integer, intent(in) :: l
    integer, intent(in) :: e
    logical :: is_inside


    integer :: total_size
    integer :: id


    real(kind=dp), dimension(3) :: p
    real(kind=dp) :: distance


    p(1) = x
    p(2) = y
    p(3) = z

    distance = signed_distance(search_tree, mesh%el, p, max_distance)







    ! Inside if distance from center <= radius
    is_inside = (distance .le. this%threshold)

  end function mesh_point_zone_criterion

end module mesh_point_zone
