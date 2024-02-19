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
module filters
  use field, only: field_t
  use num_types, only: rp
  implicit none

  private
  public :: smooth_step_field, permeability_field, step_function_field

contains

  !> @brief Apply a smooth step function to a field.
  !! @details The smooth step function is defined as:
  !! \f[
  !! t = (x - edge0) / (edge1 - edge0)
  !!  f(t) = \begin{cases}
  !!            t^3 (t (6x - 15) + 10), & t \in [0, 1] \\
  !!              0, & t \leq 0 \\
  !!              1, & t \geq 1 \\
  !!          \end{cases}
  !! \f]
  !! @note The step can be inverted by swapping edge0 and edge1.
  !!
  !! @param[in,out] F Field to be modified.
  !! @param[in] edge0 Edge giving output 0.
  !! @param[in] edge1 Edge giving output 1.
  subroutine smooth_step_field(F, edge0, edge1)
    type(field_t), pointer, intent(inout) :: F
    real(kind=rp), intent(in) :: edge0, edge1

    F%x = smooth_step_cpu(F%x, edge0, edge1)
  end subroutine smooth_step_field

  !> @brief Apply a permeability function to a field.
  !! @details The permeability function is defined as:
  !! \f[ k(x) = k_0 + (k_1 - k_0) x \frac{penalty + 1}{penalty + x}} \f]
  !! @param[in,out] F Field to be modified.
  !! @param[in] perm_0 Permeability at x=0.
  !! @param[in] perm_1 Permeability at x=1.
  !! @param[in] penalty Penalty factor.
  subroutine permeability_field(F, perm_0, perm_1, penalty)
    type(field_t), pointer, intent(inout) :: F
    real(kind=rp), intent(in) :: perm_0, perm_1
    real(kind=rp), intent(in) :: penalty

    F%x = permeability_cpu(F%x, perm_0, perm_1, penalty)
  end subroutine permeability_field

  !> @brief Apply a step function to a field.
  !! @param[in,out] F Field to be modified.
  !! @param[in] x0 Position of the step.
  !! @param[in] value0 Value of the field before the step.
  !! @param[in] value1 Value of the field after the step.
  subroutine step_function_field(F, x0, value0, value1)
    type(field_t), pointer, intent(inout) :: F
    real(kind=rp), intent(in) :: x0, value0, value1

    F%x = step_function_cpu(F%x, x0, value0, value1)
  end subroutine step_function_field

  ! ========================================================================== !
  ! Internal functions and subroutines
  ! ========================================================================== !

  !> @brief Apply a smooth step function to a scalar.
  elemental function smooth_step_cpu(x, edge0, edge1) result(res)
    real(kind=rp), intent(in) :: x
    real(kind=rp), intent(in) :: edge0, edge1
    real(kind=rp) :: res, t

    t = clamp((x - edge0) / (edge1 - edge0), 0.0_rp, 1.0_rp)

    res = t**3 * (t * (6.0_rp * t - 15.0_rp) + 10.0_rp)

  end function smooth_step_cpu

  !> @brief Clamp a value between two limits.
  elemental function clamp(x, lowerlimit, upperlimit) result(res)
    real(kind=rp), intent(in) :: x
    real(kind=rp), intent(in) :: lowerlimit, upperlimit
    real(kind=rp) :: res

    res = max(lowerlimit, min(upperlimit, x))
  end function clamp

  !> @brief Apply a step function to a scalar.
  elemental function step_function_cpu(x, x0, value0, value1) result(res)
    real(kind=rp), intent(in) :: x, x0, value0, value1
    real(kind=rp) :: res

    res = merge(value0, value1, x > x0)

  end function step_function_cpu

  !> @brief Apply a permeability function to a scalar.
  elemental function permeability_cpu(x, perm_0, perm_1, penalty) result(perm)
    real(kind=rp), intent(in) :: x, perm_0, perm_1, penalty
    real(kind=rp) :: perm

    perm = x * (penalty + 1.0_rp) / (penalty + x) * (perm_1 - perm_0) - perm_1

  end function permeability_cpu


end module filters
