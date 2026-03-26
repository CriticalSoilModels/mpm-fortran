module mpm_shapefn_linear
   !! Linear (2-node) shape functions for 1D MPM.
   !!
   !! All procedures are pure and marked !$acc routine seq so they can be
   !! called directly from inside !$acc parallel regions without modification.
   !! No allocations, no I/O — safe for GPU kernels.
   use mpm_precision, only: wp
   implicit none
   private

   public :: locate_cell_1d
   public :: eval_linear_1d

contains

   pure function locate_cell_1d(xp, x_min, dx, n_cells) result(c)
      !! Return the 1-based index of the cell containing particle position xp.
      !!
      !! For a uniform grid with cell size dx starting at x_min, the host cell
      !! is floor((xp - x_min) / dx) + 1, clamped to [1, n_cells] so that a
      !! particle exactly on the right boundary maps to the last cell rather
      !! than an out-of-bounds index.
      !$acc routine seq
      real(wp), intent(in) :: xp       !! particle position [m]
      real(wp), intent(in) :: x_min    !! domain left boundary [m]
      real(wp), intent(in) :: dx       !! uniform cell size [m]
      integer,  intent(in) :: n_cells  !! total number of cells
      integer :: c

      c = int((xp - x_min) / dx) + 1
      if (c > n_cells) c = n_cells

   end function locate_cell_1d

   pure subroutine eval_linear_1d(xp, x_a, dx, N, dN_dx)
      !! Evaluate 1D linear shape functions and their gradients at xp.
      !!
      !! The local coordinate is xi = (xp - x_a) / dx in [0, 1].
      !!
      !!   N(1) = 1 - xi        (left  node weight)
      !!   N(2) = xi            (right node weight)
      !!   dN_dx(1) = -1 / dx
      !!   dN_dx(2) = +1 / dx
      !!
      !! Partition of unity holds: N(1) + N(2) = 1, dN_dx(1) + dN_dx(2) = 0.
      !$acc routine seq
      real(wp), intent(in)  :: xp     !! particle position [m]
      real(wp), intent(in)  :: x_a    !! left node position [m]
      real(wp), intent(in)  :: dx     !! cell size [m]
      real(wp), intent(out) :: N(2)     !! shape function values    [-]
      real(wp), intent(out) :: dN_dx(2) !! shape function gradients [1/m]

      real(wp) :: xi

      xi = (xp - x_a) / dx

      N(1) = 1.0_wp - xi
      N(2) = xi

      dN_dx(1) = -1.0_wp / dx
      dN_dx(2) =  1.0_wp / dx

   end subroutine eval_linear_1d

end module mpm_shapefn_linear
