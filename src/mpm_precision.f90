module mpm_precision
   !! Working real kind for mpm-fortran.
   !!
   !! Change `wp = dp` to `wp = sp` here to switch the entire codebase to
   !! single precision. All other modules import `wp` from this module only.
   use stdlib_kinds, only: sp, dp
   implicit none
   private

   integer, parameter, public :: wp = dp

end module mpm_precision
