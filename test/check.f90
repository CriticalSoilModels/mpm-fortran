program test_runner
   !! Test-drive runner — collects and executes all test suites.
   use iso_fortran_env,                only: error_unit
   use testdrive,                      only: run_testsuite, new_testsuite, testsuite_type
   use test_mpm_problem_config_suite,  only: collect_mpm_problem_config_suite
   use test_mpm_particles_suite,       only: collect_mpm_particles_suite
   use test_mpm_update_config_suite,   only: collect_mpm_update_config_suite
   implicit none

   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
      new_testsuite("mpm_problem_config", collect_mpm_problem_config_suite), &
      new_testsuite("mpm_particles",      collect_mpm_particles_suite),      &
      new_testsuite("mpm_update_config",  collect_mpm_update_config_suite)   &
   ]

   do is = 1, size(testsuites)
      write(error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if

end program test_runner
