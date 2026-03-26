module test_mpm_shapefn_linear_suite
   !! Tests for mpm_shapefn_linear: cell location and 1D linear shape functions.
   use testdrive,             only: new_unittest, unittest_type, error_type, check
   use mpm_precision,         only: wp
   use mpm_shapefn_linear,    only: locate_cell_1d, eval_linear_1d
   implicit none
   private

   public :: collect_mpm_shapefn_linear_suite

contains

   subroutine collect_mpm_shapefn_linear_suite(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("locate_cell_interior",       test_locate_cell_interior),       &
         new_unittest("locate_cell_at_x_min",       test_locate_cell_at_x_min),       &
         new_unittest("locate_cell_at_x_max",       test_locate_cell_at_x_max),       &
         new_unittest("locate_cell_near_boundary",  test_locate_cell_near_boundary),  &
         new_unittest("eval_at_left_node",          test_eval_at_left_node),          &
         new_unittest("eval_at_right_node",         test_eval_at_right_node),         &
         new_unittest("eval_at_midpoint",           test_eval_at_midpoint),           &
         new_unittest("eval_partition_of_unity",    test_eval_partition_of_unity),    &
         new_unittest("eval_gradient_values",       test_eval_gradient_values),       &
         new_unittest("eval_gradient_sum_zero",     test_eval_gradient_sum_zero)      &
      ]
   end subroutine collect_mpm_shapefn_linear_suite

   ! ── locate_cell_1d ────────────────────────────────────────────────────────

   subroutine test_locate_cell_interior(error)
      !! Particle well inside the domain maps to the correct cell.
      type(error_type), allocatable, intent(out) :: error
      ! Grid: x in [0, 2], 4 cells of dx=0.5 → cells 1..4
      ! xp=0.6 lies in [0.5, 1.0] → cell 2
      call check(error, locate_cell_1d(0.6_wp, 0.0_wp, 0.5_wp, 4) == 2)
   end subroutine

   subroutine test_locate_cell_at_x_min(error)
      !! Particle at the left boundary maps to cell 1.
      type(error_type), allocatable, intent(out) :: error
      call check(error, locate_cell_1d(0.0_wp, 0.0_wp, 0.5_wp, 4) == 1)
   end subroutine

   subroutine test_locate_cell_at_x_max(error)
      !! Particle exactly at the right boundary is clamped to the last cell.
      type(error_type), allocatable, intent(out) :: error
      call check(error, locate_cell_1d(2.0_wp, 0.0_wp, 0.5_wp, 4) == 4)
   end subroutine

   subroutine test_locate_cell_near_boundary(error)
      !! Particle just inside the last cell maps correctly (not clamped).
      type(error_type), allocatable, intent(out) :: error
      call check(error, locate_cell_1d(1.99_wp, 0.0_wp, 0.5_wp, 4) == 4)
   end subroutine

   ! ── eval_linear_1d ────────────────────────────────────────────────────────

   subroutine test_eval_at_left_node(error)
      !! Particle sitting on the left node: N=[1,0].
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: N(2), dN_dx(2)
      call eval_linear_1d(0.0_wp, 0.0_wp, 0.5_wp, N, dN_dx)
      call check(error, abs(N(1) - 1.0_wp) < 1.0e-14_wp); if (allocated(error)) return
      call check(error, abs(N(2) - 0.0_wp) < 1.0e-14_wp)
   end subroutine

   subroutine test_eval_at_right_node(error)
      !! Particle sitting on the right node: N=[0,1].
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: N(2), dN_dx(2)
      call eval_linear_1d(0.5_wp, 0.0_wp, 0.5_wp, N, dN_dx)
      call check(error, abs(N(1) - 0.0_wp) < 1.0e-14_wp); if (allocated(error)) return
      call check(error, abs(N(2) - 1.0_wp) < 1.0e-14_wp)
   end subroutine

   subroutine test_eval_at_midpoint(error)
      !! Particle at cell midpoint: N=[0.5, 0.5].
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: N(2), dN_dx(2)
      call eval_linear_1d(0.25_wp, 0.0_wp, 0.5_wp, N, dN_dx)
      call check(error, abs(N(1) - 0.5_wp) < 1.0e-14_wp); if (allocated(error)) return
      call check(error, abs(N(2) - 0.5_wp) < 1.0e-14_wp)
   end subroutine

   subroutine test_eval_partition_of_unity(error)
      !! N(1) + N(2) == 1 at an arbitrary interior point.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: N(2), dN_dx(2)
      call eval_linear_1d(0.3_wp, 0.0_wp, 0.5_wp, N, dN_dx)
      call check(error, abs(N(1) + N(2) - 1.0_wp) < 1.0e-14_wp)
   end subroutine

   subroutine test_eval_gradient_values(error)
      !! dN_dx(1) = -1/dx,  dN_dx(2) = +1/dx.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: N(2), dN_dx(2)
      real(wp), parameter :: dx = 0.5_wp
      call eval_linear_1d(0.2_wp, 0.0_wp, dx, N, dN_dx)
      call check(error, abs(dN_dx(1) - (-1.0_wp/dx)) < 1.0e-14_wp); if (allocated(error)) return
      call check(error, abs(dN_dx(2) - ( 1.0_wp/dx)) < 1.0e-14_wp)
   end subroutine

   subroutine test_eval_gradient_sum_zero(error)
      !! Partition of unity → sum of gradients is zero.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: N(2), dN_dx(2)
      call eval_linear_1d(0.3_wp, 0.0_wp, 0.5_wp, N, dN_dx)
      call check(error, abs(dN_dx(1) + dN_dx(2)) < 1.0e-14_wp)
   end subroutine

end module test_mpm_shapefn_linear_suite
