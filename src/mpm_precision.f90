module mpm_precision
   !! Precision control and problem-type constants for mpm-fortran.
   !!
   !! All modules use `wp` as the working real kind. To switch the entire
   !! codebase to single precision, change `wp = dp` to `wp = sp` here.
   !!
   !! Supported problem types: PROB_1D, PROB_2D_PLANE_STRAIN, PROB_3D.
   !! Axisymmetric variants are deferred — see notes/future_plan.md.
   use stdlib_kinds, only: sp, dp
   implicit none
   private

   ! --- working precision ---
   integer, parameter, public :: wp = dp  !! working real kind (default: double)

   ! --- problem type constants ---
   integer, parameter, public :: PROB_1D              = 1
   integer, parameter, public :: PROB_2D_PLANE_STRAIN = 2
   integer, parameter, public :: PROB_3D              = 3

   ! --- Voigt index labels (canonical ordering [11,22,33,12,13,23]) ---
   integer, parameter, public :: V_11 = 1
   integer, parameter, public :: V_22 = 2
   integer, parameter, public :: V_33 = 3
   integer, parameter, public :: V_12 = 4
   integer, parameter, public :: V_13 = 5
   integer, parameter, public :: V_23 = 6

   public :: n_voigt
   public :: n_dims

contains

   pure function n_voigt(problem_type) result(nv)
      !! Return the number of independent stress/strain components.
      !! 1D→1, 2D plane strain→4, 3D→6.
      integer, intent(in) :: problem_type
      integer :: nv
      select case (problem_type)
         case (PROB_1D)
            nv = 1
         case (PROB_2D_PLANE_STRAIN)
            nv = 4
         case (PROB_3D)
            nv = 6
         case default
            nv = -1
      end select
   end function n_voigt

   pure function n_dims(problem_type) result(nd)
      !! Return the number of spatial dimensions.
      !! 1D→1, 2D plane strain→2, 3D→3.
      integer, intent(in) :: problem_type
      integer :: nd
      select case (problem_type)
         case (PROB_1D)
            nd = 1
         case (PROB_2D_PLANE_STRAIN)
            nd = 2
         case (PROB_3D)
            nd = 3
         case default
            nd = -1
      end select
   end function n_dims

end module mpm_precision
