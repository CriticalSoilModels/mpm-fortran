program mpm
   !! MPM solver entry point.
   use mpm_fortran, only: wp, PROB_1D
   implicit none

   write(*, '(a)') "mpm-fortran: Material Point Method solver"
   write(*, '(a)') "Run `fpm test` to execute the test suite."

end program mpm
