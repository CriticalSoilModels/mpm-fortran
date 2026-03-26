module test_mpm_particles_suite
   !! Tests for mpm_particles: allocation shapes, zero-init, deallocation.
   use testdrive,    only: new_unittest, unittest_type, error_type, check
   use mpm_precision,      only: wp
   use mpm_problem_config, only: PROB_1D, PROB_3D
   use mpm_particles, only: mpm_particle_set_t, &
                             allocate_particle_set, deallocate_particle_set
   implicit none
   private

   public :: collect_mpm_particles_suite

contains

   subroutine collect_mpm_particles_suite(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("alloc_1d_shape",        test_alloc_1d_shape),        &
         new_unittest("alloc_3d_shape",        test_alloc_3d_shape),        &
         new_unittest("alloc_zero_init",       test_alloc_zero_init),       &
         new_unittest("dealloc_resets_counts", test_dealloc_resets_counts)  &
      ]
   end subroutine collect_mpm_particles_suite

   subroutine test_alloc_1d_shape(error)
      type(error_type), allocatable, intent(out) :: error
      type(mpm_particle_set_t) :: set
      call allocate_particle_set(set, 10, PROB_1D, 0)
      call check(error, set%n_mp   == 10); if (allocated(error)) return
      call check(error, set%ndim   == 1);  if (allocated(error)) return
      call check(error, set%nvoigt == 1);  if (allocated(error)) return
      call check(error, size(set%xp,  1) == 1);  if (allocated(error)) return
      call check(error, size(set%xp,  2) == 10); if (allocated(error)) return
      call check(error, size(set%sig, 1) == 1);  if (allocated(error)) return
      call check(error, size(set%sig, 2) == 10)
      call deallocate_particle_set(set)
   end subroutine

   subroutine test_alloc_3d_shape(error)
      type(error_type), allocatable, intent(out) :: error
      type(mpm_particle_set_t) :: set
      call allocate_particle_set(set, 5, PROB_3D, 3)
      call check(error, set%ndim   == 3); if (allocated(error)) return
      call check(error, set%nvoigt == 6); if (allocated(error)) return
      call check(error, size(set%xp,     1) == 3); if (allocated(error)) return
      call check(error, size(set%sig,    1) == 6); if (allocated(error)) return
      call check(error, size(set%statev, 1) == 3)
      call deallocate_particle_set(set)
   end subroutine

   subroutine test_alloc_zero_init(error)
      type(error_type), allocatable, intent(out) :: error
      type(mpm_particle_set_t) :: set
      call allocate_particle_set(set, 5, PROB_1D, 0)
      call check(error, all(set%xp  == 0.0_wp)); if (allocated(error)) return
      call check(error, all(set%vp  == 0.0_wp)); if (allocated(error)) return
      call check(error, all(set%mp  == 0.0_wp)); if (allocated(error)) return
      call check(error, all(set%sig == 0.0_wp))
      call deallocate_particle_set(set)
   end subroutine

   subroutine test_dealloc_resets_counts(error)
      type(error_type), allocatable, intent(out) :: error
      type(mpm_particle_set_t) :: set
      call allocate_particle_set(set, 10, PROB_3D, 2)
      call deallocate_particle_set(set)
      call check(error, set%n_mp   == 0); if (allocated(error)) return
      call check(error, set%ndim   == 0); if (allocated(error)) return
      call check(error, set%nvoigt == 0)
   end subroutine

end module test_mpm_particles_suite
