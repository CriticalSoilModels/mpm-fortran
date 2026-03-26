module test_mpm_problem_config_suite
   !! Tests for mpm_problem_config: problem type constants, n_voigt, n_dims.
   use testdrive,          only: new_unittest, unittest_type, error_type, check
   use mpm_problem_config, only: n_voigt, n_dims, &
                                  PROB_1D, PROB_2D_PLANE_STRAIN, PROB_3D
   implicit none
   private

   public :: collect_mpm_problem_config_suite

contains

   subroutine collect_mpm_problem_config_suite(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("n_voigt_1d",             test_n_voigt_1d),              &
         new_unittest("n_voigt_2d_plane_strain", test_n_voigt_2d_plane_strain), &
         new_unittest("n_voigt_3d",             test_n_voigt_3d),              &
         new_unittest("n_dims_1d",              test_n_dims_1d),               &
         new_unittest("n_dims_2d_plane_strain", test_n_dims_2d_plane_strain),  &
         new_unittest("n_dims_3d",              test_n_dims_3d)                &
      ]
   end subroutine collect_mpm_problem_config_suite

   subroutine test_n_voigt_1d(error)
      type(error_type), allocatable, intent(out) :: error
      call check(error, n_voigt(PROB_1D) == 1)
   end subroutine

   subroutine test_n_voigt_2d_plane_strain(error)
      type(error_type), allocatable, intent(out) :: error
      call check(error, n_voigt(PROB_2D_PLANE_STRAIN) == 4)
   end subroutine

   subroutine test_n_voigt_3d(error)
      type(error_type), allocatable, intent(out) :: error
      call check(error, n_voigt(PROB_3D) == 6)
   end subroutine

   subroutine test_n_dims_1d(error)
      type(error_type), allocatable, intent(out) :: error
      call check(error, n_dims(PROB_1D) == 1)
   end subroutine

   subroutine test_n_dims_2d_plane_strain(error)
      type(error_type), allocatable, intent(out) :: error
      call check(error, n_dims(PROB_2D_PLANE_STRAIN) == 2)
   end subroutine

   subroutine test_n_dims_3d(error)
      type(error_type), allocatable, intent(out) :: error
      call check(error, n_dims(PROB_3D) == 3)
   end subroutine

end module test_mpm_problem_config_suite
