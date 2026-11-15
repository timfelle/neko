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
!> @brief Module containing Signed Distance Functions.
!! @details The sign of the distance to a triangle mesh follows the winding of
!! its triangles: with every triangle wound counter-clockwise when seen from
!! outside, the distance is negative inside the mesh and positive outside.
!! The sign is taken from the angle-weighted pseudo-normal of the mesh feature
!! (face, edge or vertex) that holds the closest point, which gives the right
!! sign for every point when the mesh is closed and consistently oriented
!! (J. A. Baerentzen and H. Aanaes, Signed distance computation using the
!! angle weighted pseudonormal, IEEE TVCG 11(3), 2005). Edges and vertices are
!! shared between triangles whose vertex coordinates are identical.
module signed_distance
  use num_types, only : dp, rp
  use field, only : field_t
  use tri, only : tri_t
  use tri_mesh, only : tri_mesh_t
  use point, only : point_t
  use tuple, only : tuple_i4_t
  use htable, only : htable_pt_t, htable_i4t2_t
  use aabb, only : aabb_t
  use aabb_tree, only : aabb_tree_t, aabb_node_t, AABB_NULL_NODE
  use device, only : device_memcpy, HOST_TO_DEVICE
  use neko_config, only : NEKO_BCKND_DEVICE
  use utils, only : neko_error, neko_warning
  use math, only : rzero
  use stack, only : stack_i4_t

  implicit none
  private
  public :: signed_distance_field

  !> Triangles with twice their area below this are degenerate and ignored.
  real(kind=dp), parameter :: degenerate_tol = 1e-10_dp

  !> Angle-weighted pseudo-normals of a triangle mesh.
  !! @details A feature of a triangle is encoded as 0 for its face, k = 1..3
  !! for its edge from vertex k to vertex mod(k, 3) + 1 and 3 + k for its
  !! vertex k. The pseudo-normal of an edge is the sum of the unit normals of
  !! the triangles sharing it, the one of a vertex the sum of the unit normals
  !! of the triangles sharing it, weighted by their angle at the vertex.
  type :: pseudo_normals_t
     !> Unit normal of each triangle, zero for a degenerate one
     real(kind=dp), allocatable :: face(:,:)
     !> Pseudo-normal of each distinct edge
     real(kind=dp), allocatable :: edge(:,:)
     !> Pseudo-normal of each distinct vertex
     real(kind=dp), allocatable :: vertex(:,:)
     !> Distinct edge of each local edge of each triangle
     integer, allocatable :: tri_edge(:,:)
     !> Distinct vertex of each local vertex of each triangle
     integer, allocatable :: tri_vertex(:,:)
   contains
     procedure, pass(this) :: init => pseudo_normals_init
     procedure, pass(this) :: free => pseudo_normals_free
     procedure, pass(this) :: get => pseudo_normals_get
  end type pseudo_normals_t

contains

  !> @brief Signed distance field
  !! @details This routine computes the signed distance field for a given
  !! object.
  !!
  !! Currently supported objects are:
  !! - Triangular mesh (tri_mesh_t)
  !!
  !! @param[in,out] field_data Field data
  !! @param[in] object Object
  !! @param[in,optional] max_distance Maximum distance outside the mesh
  subroutine signed_distance_field(field_data, object, max_distance)
    type(field_t), intent(inout) :: field_data
    class(*), intent(in) :: object
    real(kind=dp), intent(in), optional :: max_distance

    real(kind=dp) :: max_dist

    if (present(max_distance)) then
       max_dist = max_distance
    else
       max_dist = huge(0.0_dp)
    end if

    select type (object)
    type is (tri_mesh_t)
       call signed_distance_field_tri_mesh(field_data, object, max_dist)

    class default
       call neko_error("signed_distance_field: Object type not supported.")
    end select

  end subroutine signed_distance_field

  !> @brief Signed distance field for a triangular mesh
  !! @details This routine computes the signed distance field for a given
  !! triangular mesh. The algorithm utilizes an AABB tree to accelerate the
  !! earch for potential elements. The signed distance is computed using the
  !! brute force approach, where we compute the distance to each element
  !! found through the AABB tree, and return the minimum distance, signed by
  !! the pseudo-normal at the closest point.
  !!
  !! @param[in,out] field_data Field data
  !! @param[in] mesh Triangular mesh
  !! @param[in] max_distance Maximum distance outside the mesh
  subroutine signed_distance_field_tri_mesh(field_data, mesh, max_distance)
    type(field_t), intent(inout) :: field_data
    type(tri_mesh_t), intent(in) :: mesh
    real(kind=dp), intent(in) :: max_distance

    integer :: total_size
    integer :: id
    type(aabb_tree_t) :: search_tree
    type(pseudo_normals_t) :: normals
    real(kind=dp), dimension(3) :: p
    real(kind=dp) :: distance

    ! Zero the field
    call rzero(field_data%x, field_data%size())
    total_size = field_data%dof%size()

    call search_tree%init(mesh%nelv)
    call search_tree%build(mesh%el)

    if (search_tree%get_size() .ne. mesh%nelv) then
       call neko_error("signed_distance_field_tri_mesh: " // &
            "Error building the search tree.")
    end if

    call normals%init(mesh)

    !$omp parallel do private(id, p, distance)
    do id = 1, total_size
       p(1) = field_data%dof%x%x(id, 1, 1, 1)
       p(2) = field_data%dof%y%x(id, 1, 1, 1)
       p(3) = field_data%dof%z%x(id, 1, 1, 1)

       distance = tri_mesh_aabb_tree(search_tree, mesh%el, normals, p, &
            max_distance)

       field_data%x(id, 1, 1, 1) = real(distance, kind=rp)
    end do
    !$omp end parallel do

    call normals%free()

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_warning("signed_distance_field_tri_mesh:&
       & Device version not implemented.")
       call device_memcpy(field_data%x, field_data%x_d, field_data%size(), &
            HOST_TO_DEVICE, sync = .false.)
    end if

  end subroutine signed_distance_field_tri_mesh

  !> @brief Signed distance function
  !! @details This routine computes the signed distance function for the
  !! boundary mesh, to a given point (x, y, z). The algorithm is a
  !! brute force approach, where we compute the distance to each
  !! element in the mesh, and return the minimum distance.
  !!
  !! @note This routine is here for completeness, and is not optimized for
  !! performance. It is recommended to use the AABB tree version for large
  !! meshes.
  !!
  !! @param mesh Boundary mesh
  !! @param normals Pseudo-normals of the boundary mesh
  !! @param p Point
  !! @param max_distance Maximum distance outside the mesh
  !! @return Signed distance value
  function tri_mesh_brute_force(mesh, normals, p, max_distance) &
       result(distance)
    type(tri_mesh_t), intent(in) :: mesh
    type(pseudo_normals_t), intent(in) :: normals
    real(kind=dp), intent(in) :: p(3)
    real(kind=dp), intent(in) :: max_distance

    integer :: id, closest_id, closest_feature, cf
    real(kind=dp) :: distance
    real(kind=dp) :: cd, cp(3), closest_point(3)

    distance = huge(0.0_dp)
    closest_id = 0
    closest_feature = 0
    closest_point = p

    do id = 1, mesh%nelv
       call element_distance(mesh%el(id), p, cd, cp, cf)

       if (cd .lt. distance) then
          distance = cd
          closest_id = id
          closest_feature = cf
          closest_point = cp
       end if
    end do

    distance = min(distance, max_distance)
    if (closest_id .gt. 0) then
       distance = sign(distance, dot_product(p - closest_point, &
            normals%get(closest_id, closest_feature)))
    end if

  end function tri_mesh_brute_force

  !> @brief Signed distance function using an AABB tree
  !! @details This routine computes the signed distance function for the
  !! boundary mesh, to a given point (x, y, z). The algorithm utilizes an
  !! AABB tree to accelerate the search for potential elements. The signed
  !! distance is computed using the brute force approach, where we compute the
  !! distance to each element found through the AABB tree, and return the
  !! minimum distance, signed by the pseudo-normal at the closest point.
  !!
  !! @param tree AABB tree
  !! @param object_list List of objects
  !! @param normals Pseudo-normals of the objects
  !! @param p Point
  !! @param max_distance Maximum distance outside the mesh
  !! @return Signed distance value
  function tri_mesh_aabb_tree(tree, object_list, normals, p, max_distance) &
       result(distance)
    class(aabb_tree_t), intent(in) :: tree
    class(tri_t), contiguous, dimension(:), intent(in) :: object_list
    type(pseudo_normals_t), intent(in) :: normals
    real(kind=dp), dimension(3), intent(in) :: p
    real(kind=dp), intent(in) :: max_distance

    real(kind=dp) :: distance

    type(stack_i4_t) :: simple_stack
    integer :: current_index

    type(aabb_node_t) :: current_node
    type(aabb_t) :: current_aabb
    integer :: current_object_index
    real(kind=dp) :: current_distance
    real(kind=dp), dimension(3) :: current_point
    integer :: current_feature

    integer :: closest_object_index
    integer :: closest_feature
    real(kind=dp), dimension(3) :: closest_point

    type(aabb_node_t) :: left_node
    type(aabb_node_t) :: right_node

    type(aabb_t) :: root_box
    type(aabb_t) :: search_box

    integer :: root_index, left_index, right_index
    real(kind=dp) :: random_value

    ! Initialize the stack and the search box
    call simple_stack%init(size(object_list) * 2)
    call search_box%init(p - max_distance, p + max_distance)

    ! Check if the root node overlaps the search box, if it does, push it to
    ! the stack and update the search box to a randomly selected object.
    root_index = tree%get_root_index()
    root_box = tree%get_aabb(root_index)

    if (.not. root_box%overlaps(search_box)) then
       distance = max_distance
       return
    end if

    ! Grab a random object and compute the distance to it
    call random_number(random_value)
    current_object_index = floor(random_value * size(object_list) + 1)
    call element_distance(object_list(current_object_index), p, distance, &
         closest_point, closest_feature)
    closest_object_index = current_object_index
    distance = distance + object_list(current_object_index)%diameter()

    ! Update the search box to the new distance and push the root node
    call search_box%init(p - distance, p + distance)
    call simple_stack%push(root_index)

    ! Traverse the tree and compute the distance to the elements
    do while (.not. simple_stack%is_empty())
       current_index = simple_stack%pop()
       if (current_index .eq. AABB_NULL_NODE) cycle

       current_node = tree%get_node(current_index)
       current_aabb = current_node%get_aabb()

       if (current_node%is_leaf()) then
          if (distance .lt. current_node%min_distance(p)) then
             cycle
          end if

          current_object_index = current_node%get_object_index()
          call element_distance(object_list(current_object_index), p, &
               current_distance, current_point, current_feature)

          if (current_distance .lt. distance) then
             distance = current_distance
             closest_object_index = current_object_index
             closest_feature = current_feature
             closest_point = current_point
          end if

          ! Update the search box to the new distance
          if (distance .gt. current_aabb%get_diameter()) then
             call search_box%init(p - distance, p + distance)
          end if
       else

          left_index = tree%get_left_index(current_index)
          if (left_index .ne. AABB_NULL_NODE) then
             left_node = tree%get_left_node(current_index)
             if (left_node%aabb%overlaps(search_box)) then
                call simple_stack%push(left_index)
             end if
          end if

          right_index = tree%get_right_index(current_index)
          if (right_index .ne. AABB_NULL_NODE) then
             right_node = tree%get_right_node(current_index)
             if (right_node%aabb%overlaps(search_box)) then
                call simple_stack%push(right_index)
             end if
          end if
       end if
    end do

    if (distance .gt. max_distance) then
       distance = max_distance
    end if
    distance = sign(distance, dot_product(p - closest_point, &
         normals%get(closest_object_index, closest_feature)))

  end function tri_mesh_aabb_tree

  ! ========================================================================== !
  ! Pseudo-normals
  ! ========================================================================== !

  !> Compute the pseudo-normals of a triangle mesh
  !! @param this Pseudo-normals
  !! @param mesh Triangular mesh
  subroutine pseudo_normals_init(this, mesh)
    class(pseudo_normals_t), intent(inout) :: this
    type(tri_mesh_t), intent(in) :: mesh

    type(htable_pt_t) :: vertex_table
    type(htable_i4t2_t) :: edge_table
    type(point_t) :: vertex_key
    type(tuple_i4_t) :: edge_key
    real(kind=dp) :: v(3, 3)
    real(kind=dp), dimension(3) :: normal, e1, e2
    real(kind=dp) :: normal_length, angle
    integer :: i, k, id, a, b, n_vertices, n_edges

    call this%free()

    allocate(this%face(3, mesh%nelv))
    allocate(this%tri_edge(3, mesh%nelv))
    allocate(this%tri_vertex(3, mesh%nelv))

    ! Number the distinct vertices. Adding zero turns -0 into +0, which the
    ! hash of a point would otherwise tell apart.
    id = 0
    call vertex_table%init(3 * mesh%nelv, id)
    n_vertices = 0
    do i = 1, mesh%nelv
       do k = 1, 3
          call vertex_key%init(mesh%el(i)%pts(k)%p%x + 0.0_dp)
          if (vertex_table%get(vertex_key, id) .gt. 0) then
             n_vertices = n_vertices + 1
             id = n_vertices
             call vertex_table%set(vertex_key, id)
          end if
          this%tri_vertex(k, i) = id
       end do
    end do
    call vertex_table%free()

    ! Number the distinct edges by their pair of vertices
    call edge_table%init(3 * mesh%nelv, id)
    n_edges = 0
    do i = 1, mesh%nelv
       do k = 1, 3
          a = this%tri_vertex(k, i)
          b = this%tri_vertex(mod(k, 3) + 1, i)
          edge_key%x = [min(a, b), max(a, b)]
          if (edge_table%get(edge_key, id) .gt. 0) then
             n_edges = n_edges + 1
             id = n_edges
             call edge_table%set(edge_key, id)
          end if
          this%tri_edge(k, i) = id
       end do
    end do
    call edge_table%free()

    allocate(this%edge(3, n_edges))
    allocate(this%vertex(3, n_vertices))
    this%edge = 0.0_dp
    this%vertex = 0.0_dp

    do i = 1, mesh%nelv
       do k = 1, 3
          v(:, k) = mesh%el(i)%pts(k)%p%x
       end do

       normal = cross(v(:, 2) - v(:, 1), v(:, 3) - v(:, 1))
       normal_length = norm2(normal)
       if (normal_length .lt. degenerate_tol) then
          this%face(:, i) = 0.0_dp
          cycle
       end if
       normal = normal / normal_length
       this%face(:, i) = normal

       do k = 1, 3
          e1 = v(:, mod(k, 3) + 1) - v(:, k)
          e2 = v(:, mod(k + 1, 3) + 1) - v(:, k)
          angle = atan2(norm2(cross(e1, e2)), dot_product(e1, e2))

          this%edge(:, this%tri_edge(k, i)) = &
               this%edge(:, this%tri_edge(k, i)) + normal
          this%vertex(:, this%tri_vertex(k, i)) = &
               this%vertex(:, this%tri_vertex(k, i)) + angle * normal
       end do
    end do

  end subroutine pseudo_normals_init

  !> Deallocate pseudo-normals
  !! @param this Pseudo-normals
  subroutine pseudo_normals_free(this)
    class(pseudo_normals_t), intent(inout) :: this

    if (allocated(this%face)) deallocate(this%face)
    if (allocated(this%edge)) deallocate(this%edge)
    if (allocated(this%vertex)) deallocate(this%vertex)
    if (allocated(this%tri_edge)) deallocate(this%tri_edge)
    if (allocated(this%tri_vertex)) deallocate(this%tri_vertex)

  end subroutine pseudo_normals_free

  !> Pseudo-normal of a feature of a triangle
  !! @param this Pseudo-normals
  !! @param element Index of the triangle
  !! @param feature Feature of the triangle, see pseudo_normals_t
  !! @return Pseudo-normal, not normalised
  pure function pseudo_normals_get(this, element, feature) result(normal)
    class(pseudo_normals_t), intent(in) :: this
    integer, intent(in) :: element
    integer, intent(in) :: feature
    real(kind=dp), dimension(3) :: normal

    select case (feature)
    case (1:3)
       normal = this%edge(:, this%tri_edge(feature, element))
    case (4:6)
       normal = this%vertex(:, this%tri_vertex(feature - 3, element))
    case default
       normal = this%face(:, element)
    end select

  end function pseudo_normals_get

  ! ========================================================================== !
  ! Element distance functions
  ! ========================================================================== !

  !> @brief Main interface for the distance function for an element.
  !! @param element Element
  !! @param p Point
  !! @param distance Distance from the point to the element
  !! @param closest Closest point on the element
  !! @param feature Feature of the element that holds the closest point
  subroutine element_distance(element, p, distance, closest, feature)
    class(*), intent(in) :: element
    real(kind=dp), dimension(3), intent(in) :: p
    real(kind=dp), intent(out) :: distance
    real(kind=dp), dimension(3), intent(out) :: closest
    integer, intent(out) :: feature

    select type (element)
    type is (tri_t)
       call element_distance_triangle(element, p, distance, closest, feature)

    class default
       print *, "Error: Element type not supported."
       stop
    end select
  end subroutine element_distance

  ! -------------------------------------------------------------------------- !
  ! Specific element distance functions

  !> @brief Distance function for a triangle
  !! @details This routine computes the closest point of the triangle to the
  !! point and the feature of the triangle that holds it, by testing the
  !! Voronoi regions of the vertices, the edges and the face in turn
  !! (C. Ericson, Real-Time Collision Detection, 2005, section 5.1.5).
  !! A degenerate triangle is reported at the largest representable distance.
  !!
  !! @param triangle Triangle
  !! @param p Point
  !! @param distance Distance from the point to the triangle
  !! @param closest Closest point on the triangle
  !! @param feature Feature of the triangle that holds the closest point,
  !! encoded as in pseudo_normals_t
  subroutine element_distance_triangle(triangle, p, distance, closest, &
       feature)
    type(tri_t), intent(in) :: triangle
    real(kind=dp), dimension(3), intent(in) :: p
    real(kind=dp), intent(out) :: distance
    real(kind=dp), dimension(3), intent(out) :: closest
    integer, intent(out) :: feature

    real(kind=dp), dimension(3) :: a, b, c, ab, ac
    real(kind=dp) :: d1, d2, d3, d4, d5, d6, va, vb, vc

    a = triangle%pts(1)%p%x
    b = triangle%pts(2)%p%x
    c = triangle%pts(3)%p%x
    ab = b - a
    ac = c - a

    if (norm2(cross(ab, ac)) .lt. degenerate_tol) then
       distance = huge(1.0_dp)
       closest = a
       feature = 0
       return
    end if

    d1 = dot_product(ab, p - a)
    d2 = dot_product(ac, p - a)
    d3 = dot_product(ab, p - b)
    d4 = dot_product(ac, p - b)
    d5 = dot_product(ab, p - c)
    d6 = dot_product(ac, p - c)
    va = d3 * d6 - d5 * d4
    vb = d5 * d2 - d1 * d6
    vc = d1 * d4 - d3 * d2

    if (d1 .le. 0.0_dp .and. d2 .le. 0.0_dp) then
       closest = a
       feature = 4
    else if (d3 .ge. 0.0_dp .and. d4 .le. d3) then
       closest = b
       feature = 5
    else if (vc .le. 0.0_dp .and. d1 .ge. 0.0_dp .and. d3 .le. 0.0_dp) then
       closest = a + d1 / (d1 - d3) * ab
       feature = 1
    else if (d6 .ge. 0.0_dp .and. d5 .le. d6) then
       closest = c
       feature = 6
    else if (vb .le. 0.0_dp .and. d2 .ge. 0.0_dp .and. d6 .le. 0.0_dp) then
       closest = a + d2 / (d2 - d6) * ac
       feature = 3
    else if (va .le. 0.0_dp .and. d4 - d3 .ge. 0.0_dp .and. &
         d5 - d6 .ge. 0.0_dp) then
       closest = b + (d4 - d3) / ((d4 - d3) + (d5 - d6)) * (c - b)
       feature = 2
    else
       closest = a + (vb * ab + vc * ac) / (va + vb + vc)
       feature = 0
    end if

    distance = norm2(p - closest)

  end subroutine element_distance_triangle

  !> Compute cross product of two vectors
  !> @param[in] a First vector
  !> @param[in] b Second vector
  !> @return Cross product \f$ a \times b \f$
  pure function cross(a, b) result(c)
    real(kind=dp), dimension(3), intent(in) :: a
    real(kind=dp), dimension(3), intent(in) :: b
    real(kind=dp), dimension(3) :: c

    c(1) = a(2) * b(3) - a(3) * b(2)
    c(2) = a(3) * b(1) - a(1) * b(3)
    c(3) = a(1) * b(2) - a(2) * b(1)

  end function cross

end module signed_distance
