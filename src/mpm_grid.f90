module mpm_grid
   !! Background grid storage and operations.
   !!
   !! Arrays use the same SOA convention as mpm_particles: first index is the
   !! component, second is the node number. This keeps all per-node data for
   !! a given component contiguous in memory (Fortran column-major).
   use mpm_precision,      only: wp
   use mpm_problem_config, only: n_dims, PROB_1D, PROB_2D_PLANE_STRAIN, PROB_3D, MASS_TOL
   implicit none
   private

   public :: mpm_grid_t
   public :: nodes_per_cell
   public :: allocate_grid, deallocate_grid
   public :: build_regular_grid_1d
   public :: reset_grid, apply_bcs, update_grid

   type :: mpm_grid_t
      !! Background mesh and nodal field storage.

      ! --- configuration (set at allocation, read-only thereafter) ---
      integer  :: n_nodes = 0   !! total number of nodes
      integer  :: n_cells = 0   !! total number of cells
      integer  :: ndim    = 0   !! spatial dimension
      integer  :: npc     = 0   !! nodes per cell (2 / 4 / 8)

      ! --- geometry (set by build_regular_grid_*, fixed thereafter) ---
      real(wp) :: dx    = 0.0_wp   !! uniform cell size [m]
      real(wp) :: x_min = 0.0_wp   !! domain lower bound [m]
      real(wp) :: x_max = 0.0_wp   !! domain upper bound [m]

      real(wp), allocatable :: xI(:,:)          !! nodal positions [m]       (ndim, n_nodes)
      integer,  allocatable :: cell_nodes(:,:)  !! connectivity              (npc,  n_cells)

      ! --- nodal fields (zeroed each step by reset_grid) ---
      real(wp), allocatable :: mass_I(:)        !! nodal mass     [kg]       (n_nodes)
      real(wp), allocatable :: momentum_I(:,:)  !! nodal momentum [kg·m/s]   (ndim, n_nodes)
      real(wp), allocatable :: force_I(:,:)     !! nodal force    [N]        (ndim, n_nodes)
      real(wp), allocatable :: vel_I(:,:)       !! nodal velocity [m/s]      (ndim, n_nodes)
      real(wp), allocatable :: dv_I(:,:)        !! velocity increment [m/s]  (ndim, n_nodes)
                                                !! = vel_I^{n+1} - vel_I^n; used by FLIP G2P

      ! --- boundary conditions ---
      logical, allocatable  :: fixed_dof(:,:)   !! fixed degree of freedom   (ndim, n_nodes)
                                                !! .true. → zero force and momentum here
   end type mpm_grid_t

contains

   pure function nodes_per_cell(problem_type) result(npc)
      !! Number of nodes per cell: 1D → 2,  2D plane strain → 4,  3D → 8.
      integer, intent(in) :: problem_type
      integer :: npc
      select case (problem_type)
         case (PROB_1D)             ; npc = 2
         case (PROB_2D_PLANE_STRAIN); npc = 4
         case (PROB_3D)             ; npc = 8
         case default               ; npc = -1
      end select
   end function nodes_per_cell

   subroutine allocate_grid(grid, n_cells, problem_type)
      !! Allocate all arrays in a grid and initialise to zero / .false.
      !!
      !! Only the 1D formula n_nodes = n_cells + 1 is used here.
      !! For 2D/3D call build_regular_grid_2d/3d after allocation.
      type(mpm_grid_t), intent(out) :: grid
      integer,          intent(in)  :: n_cells
      integer,          intent(in)  :: problem_type

      integer :: nd, npc_val, nn

      nd      = n_dims(problem_type)
      npc_val = nodes_per_cell(problem_type)
      nn      = n_cells + 1   ! 1D; 2D/3D builders will override n_nodes

      grid%n_cells = n_cells
      grid%n_nodes = nn
      grid%ndim    = nd
      grid%npc     = npc_val

      allocate(grid%xI(nd, nn),              source=0.0_wp)
      allocate(grid%cell_nodes(npc_val, n_cells), source=0)
      allocate(grid%mass_I(nn),              source=0.0_wp)
      allocate(grid%momentum_I(nd, nn),      source=0.0_wp)
      allocate(grid%force_I(nd, nn),         source=0.0_wp)
      allocate(grid%vel_I(nd, nn),           source=0.0_wp)
      allocate(grid%dv_I(nd, nn),            source=0.0_wp)
      allocate(grid%fixed_dof(nd, nn),       source=.false.)

   end subroutine allocate_grid

   subroutine deallocate_grid(grid)
      !! Deallocate all arrays and reset scalar fields to zero.
      type(mpm_grid_t), intent(inout) :: grid

      if (allocated(grid%xI))         deallocate(grid%xI)
      if (allocated(grid%cell_nodes)) deallocate(grid%cell_nodes)
      if (allocated(grid%mass_I))     deallocate(grid%mass_I)
      if (allocated(grid%momentum_I)) deallocate(grid%momentum_I)
      if (allocated(grid%force_I))    deallocate(grid%force_I)
      if (allocated(grid%vel_I))      deallocate(grid%vel_I)
      if (allocated(grid%dv_I))       deallocate(grid%dv_I)
      if (allocated(grid%fixed_dof))  deallocate(grid%fixed_dof)

      grid%n_nodes = 0
      grid%n_cells = 0
      grid%ndim    = 0
      grid%npc     = 0
      grid%dx      = 0.0_wp
      grid%x_min   = 0.0_wp
      grid%x_max   = 0.0_wp

   end subroutine deallocate_grid

   subroutine build_regular_grid_1d(grid, x_min, x_max, n_cells)
      !! Fill xI and cell_nodes for a uniform 1D mesh.
      !!
      !! Node i sits at x_min + (i-1)*dx.
      !! Cell c connects nodes [c, c+1].
      type(mpm_grid_t), intent(inout) :: grid
      real(wp),         intent(in)    :: x_min
      real(wp),         intent(in)    :: x_max
      integer,          intent(in)    :: n_cells

      integer  :: i, c
      real(wp) :: dx

      dx = (x_max - x_min) / real(n_cells, wp)

      grid%dx    = dx
      grid%x_min = x_min
      grid%x_max = x_max

      do i = 1, grid%n_nodes
         grid%xI(1, i) = x_min + real(i - 1, wp) * dx
      end do

      do c = 1, n_cells
         grid%cell_nodes(1, c) = c
         grid%cell_nodes(2, c) = c + 1
      end do

   end subroutine build_regular_grid_1d

   subroutine reset_grid(grid)
      !! Zero all nodal fields. Called at the start of every time step.
      !!
      !! Does NOT touch xI, cell_nodes, fixed_dof, or geometry scalars.
      type(mpm_grid_t), intent(inout) :: grid

      grid%mass_I     = 0.0_wp
      grid%momentum_I = 0.0_wp
      grid%force_I    = 0.0_wp
      grid%vel_I      = 0.0_wp
      grid%dv_I       = 0.0_wp

   end subroutine reset_grid

   subroutine apply_bcs(grid)
      !! Zero force and momentum on fixed degrees of freedom.
      !!
      !! Called twice per USL step: once before and once after the momentum update.
      type(mpm_grid_t), intent(inout) :: grid

      integer :: i, n

      do n = 1, grid%n_nodes
         do i = 1, grid%ndim
            if (grid%fixed_dof(i, n)) then
               grid%force_I(i, n)    = 0.0_wp
               grid%momentum_I(i, n) = 0.0_wp
            end if
         end do
      end do

   end subroutine apply_bcs

   subroutine update_grid(grid, dt)
      !! Advance nodal momentum; compute new velocity and velocity increment.
      !!
      !! For each non-empty node:
      !!   vel_old      = momentum_I / mass_I
      !!   momentum_I  += force_I * dt
      !!   vel_I        = momentum_I / mass_I
      !!   dv_I         = vel_I - vel_old      (used by FLIP in G2P)
      !!
      !! Nodes with mass_I < MASS_TOL are skipped to avoid division by zero.
      type(mpm_grid_t), intent(inout) :: grid
      real(wp),         intent(in)    :: dt

      integer  :: i, n
      real(wp) :: vel_old

      do n = 1, grid%n_nodes
         if (grid%mass_I(n) < MASS_TOL) cycle
         do i = 1, grid%ndim
            vel_old              = grid%momentum_I(i, n) / grid%mass_I(n)
            grid%momentum_I(i,n) = grid%momentum_I(i,n) + grid%force_I(i,n) * dt
            grid%vel_I(i,n)      = grid%momentum_I(i,n) / grid%mass_I(n)
            grid%dv_I(i,n)       = grid%vel_I(i,n) - vel_old
         end do
      end do

   end subroutine update_grid

end module mpm_grid
