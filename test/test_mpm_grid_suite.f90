module test_mpm_grid_suite
   !! Tests for mpm_grid: allocation shapes, geometry, reset, BCs, momentum update.
   use testdrive,          only: new_unittest, unittest_type, error_type, check
   use mpm_precision,      only: wp
   use mpm_problem_config, only: PROB_1D, PROB_2D_PLANE_STRAIN, PROB_3D, MASS_TOL
   use mpm_grid,           only: mpm_grid_t, nodes_per_cell, &
                                  allocate_grid, deallocate_grid, &
                                  build_regular_grid_1d, &
                                  reset_grid, apply_bcs, update_grid
   implicit none
   private

   public :: collect_mpm_grid_suite

contains

   subroutine collect_mpm_grid_suite(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("nodes_per_cell_1d",           test_nodes_per_cell_1d),           &
         new_unittest("nodes_per_cell_2d",           test_nodes_per_cell_2d),           &
         new_unittest("nodes_per_cell_3d",           test_nodes_per_cell_3d),           &
         new_unittest("alloc_1d_counts",             test_alloc_1d_counts),             &
         new_unittest("alloc_1d_shapes",             test_alloc_1d_shapes),             &
         new_unittest("alloc_zero_init",             test_alloc_zero_init),             &
         new_unittest("alloc_fixed_dof_false",       test_alloc_fixed_dof_false),       &
         new_unittest("build_1d_node_positions",     test_build_1d_node_positions),     &
         new_unittest("build_1d_connectivity",       test_build_1d_connectivity),       &
         new_unittest("build_1d_geometry_scalars",   test_build_1d_geometry_scalars),   &
         new_unittest("reset_zeroes_nodal_fields",   test_reset_zeroes_nodal_fields),   &
         new_unittest("reset_preserves_geometry",    test_reset_preserves_geometry),    &
         new_unittest("apply_bcs_zeroes_fixed_dofs", test_apply_bcs_zeroes_fixed_dofs), &
         new_unittest("apply_bcs_leaves_free_dofs",  test_apply_bcs_leaves_free_dofs),  &
         new_unittest("update_grid_momentum",        test_update_grid_momentum),        &
         new_unittest("update_grid_vel_and_dv",      test_update_grid_vel_and_dv),      &
         new_unittest("update_grid_skips_empty",     test_update_grid_skips_empty)      &
      ]
   end subroutine collect_mpm_grid_suite

   ! ── nodes_per_cell ────────────────────────────────────────────────────────

   subroutine test_nodes_per_cell_1d(error)
      type(error_type), allocatable, intent(out) :: error
      call check(error, nodes_per_cell(PROB_1D) == 2)
   end subroutine

   subroutine test_nodes_per_cell_2d(error)
      type(error_type), allocatable, intent(out) :: error
      call check(error, nodes_per_cell(PROB_2D_PLANE_STRAIN) == 4)
   end subroutine

   subroutine test_nodes_per_cell_3d(error)
      type(error_type), allocatable, intent(out) :: error
      call check(error, nodes_per_cell(PROB_3D) == 8)
   end subroutine

   ! ── allocate_grid ─────────────────────────────────────────────────────────

   subroutine test_alloc_1d_counts(error)
      type(error_type), allocatable, intent(out) :: error
      type(mpm_grid_t) :: grid
      call allocate_grid(grid, n_cells=10, problem_type=PROB_1D)
      call check(error, grid%n_cells == 10); if (allocated(error)) return
      call check(error, grid%n_nodes == 11); if (allocated(error)) return
      call check(error, grid%ndim    == 1);  if (allocated(error)) return
      call check(error, grid%npc     == 2)
      call deallocate_grid(grid)
   end subroutine

   subroutine test_alloc_1d_shapes(error)
      type(error_type), allocatable, intent(out) :: error
      type(mpm_grid_t) :: grid
      call allocate_grid(grid, n_cells=10, problem_type=PROB_1D)
      call check(error, size(grid%xI,         1) == 1);  if (allocated(error)) return
      call check(error, size(grid%xI,         2) == 11); if (allocated(error)) return
      call check(error, size(grid%cell_nodes,  1) == 2);  if (allocated(error)) return
      call check(error, size(grid%cell_nodes,  2) == 10); if (allocated(error)) return
      call check(error, size(grid%mass_I)         == 11); if (allocated(error)) return
      call check(error, size(grid%momentum_I, 1) == 1);  if (allocated(error)) return
      call check(error, size(grid%momentum_I, 2) == 11); if (allocated(error)) return
      call check(error, size(grid%force_I,    1) == 1);  if (allocated(error)) return
      call check(error, size(grid%force_I,    2) == 11); if (allocated(error)) return
      call check(error, size(grid%vel_I,      1) == 1);  if (allocated(error)) return
      call check(error, size(grid%vel_I,      2) == 11); if (allocated(error)) return
      call check(error, size(grid%dv_I,       1) == 1);  if (allocated(error)) return
      call check(error, size(grid%dv_I,       2) == 11); if (allocated(error)) return
      call check(error, size(grid%fixed_dof,  1) == 1);  if (allocated(error)) return
      call check(error, size(grid%fixed_dof,  2) == 11)
      call deallocate_grid(grid)
   end subroutine

   subroutine test_alloc_zero_init(error)
      type(error_type), allocatable, intent(out) :: error
      type(mpm_grid_t) :: grid
      call allocate_grid(grid, n_cells=5, problem_type=PROB_1D)
      call check(error, all(grid%xI         == 0.0_wp)); if (allocated(error)) return
      call check(error, all(grid%mass_I     == 0.0_wp)); if (allocated(error)) return
      call check(error, all(grid%momentum_I == 0.0_wp)); if (allocated(error)) return
      call check(error, all(grid%force_I    == 0.0_wp)); if (allocated(error)) return
      call check(error, all(grid%vel_I      == 0.0_wp)); if (allocated(error)) return
      call check(error, all(grid%dv_I       == 0.0_wp))
      call deallocate_grid(grid)
   end subroutine

   subroutine test_alloc_fixed_dof_false(error)
      type(error_type), allocatable, intent(out) :: error
      type(mpm_grid_t) :: grid
      call allocate_grid(grid, n_cells=5, problem_type=PROB_1D)
      call check(error, all(.not. grid%fixed_dof))
      call deallocate_grid(grid)
   end subroutine

   ! ── build_regular_grid_1d ─────────────────────────────────────────────────

   subroutine test_build_1d_node_positions(error)
      type(error_type), allocatable, intent(out) :: error
      type(mpm_grid_t) :: grid
      real(wp) :: dx
      call allocate_grid(grid, n_cells=4, problem_type=PROB_1D)
      call build_regular_grid_1d(grid, x_min=0.0_wp, x_max=2.0_wp, n_cells=4)
      dx = 0.5_wp
      call check(error, abs(grid%xI(1,1) - 0.0_wp) < 1.0e-14_wp); if (allocated(error)) return
      call check(error, abs(grid%xI(1,2) - dx     ) < 1.0e-14_wp); if (allocated(error)) return
      call check(error, abs(grid%xI(1,5) - 2.0_wp ) < 1.0e-14_wp)
      call deallocate_grid(grid)
   end subroutine

   subroutine test_build_1d_connectivity(error)
      type(error_type), allocatable, intent(out) :: error
      type(mpm_grid_t) :: grid
      call allocate_grid(grid, n_cells=4, problem_type=PROB_1D)
      call build_regular_grid_1d(grid, x_min=0.0_wp, x_max=1.0_wp, n_cells=4)
      ! first cell connects nodes 1 and 2
      call check(error, grid%cell_nodes(1,1) == 1); if (allocated(error)) return
      call check(error, grid%cell_nodes(2,1) == 2); if (allocated(error)) return
      ! last cell connects nodes 4 and 5
      call check(error, grid%cell_nodes(1,4) == 4); if (allocated(error)) return
      call check(error, grid%cell_nodes(2,4) == 5)
      call deallocate_grid(grid)
   end subroutine

   subroutine test_build_1d_geometry_scalars(error)
      type(error_type), allocatable, intent(out) :: error
      type(mpm_grid_t) :: grid
      call allocate_grid(grid, n_cells=4, problem_type=PROB_1D)
      call build_regular_grid_1d(grid, x_min=1.0_wp, x_max=3.0_wp, n_cells=4)
      call check(error, abs(grid%dx    - 0.5_wp) < 1.0e-14_wp); if (allocated(error)) return
      call check(error, abs(grid%x_min - 1.0_wp) < 1.0e-14_wp); if (allocated(error)) return
      call check(error, abs(grid%x_max - 3.0_wp) < 1.0e-14_wp)
      call deallocate_grid(grid)
   end subroutine

   ! ── reset_grid ────────────────────────────────────────────────────────────

   subroutine test_reset_zeroes_nodal_fields(error)
      type(error_type), allocatable, intent(out) :: error
      type(mpm_grid_t) :: grid
      call allocate_grid(grid, n_cells=4, problem_type=PROB_1D)
      call build_regular_grid_1d(grid, x_min=0.0_wp, x_max=1.0_wp, n_cells=4)
      ! pollute nodal fields
      grid%mass_I     = 5.0_wp
      grid%momentum_I = 3.0_wp
      grid%force_I    = 2.0_wp
      grid%vel_I      = 1.0_wp
      grid%dv_I       = 0.5_wp
      call reset_grid(grid)
      call check(error, all(grid%mass_I     == 0.0_wp)); if (allocated(error)) return
      call check(error, all(grid%momentum_I == 0.0_wp)); if (allocated(error)) return
      call check(error, all(grid%force_I    == 0.0_wp)); if (allocated(error)) return
      call check(error, all(grid%vel_I      == 0.0_wp)); if (allocated(error)) return
      call check(error, all(grid%dv_I       == 0.0_wp))
      call deallocate_grid(grid)
   end subroutine

   subroutine test_reset_preserves_geometry(error)
      type(error_type), allocatable, intent(out) :: error
      type(mpm_grid_t) :: grid
      call allocate_grid(grid, n_cells=4, problem_type=PROB_1D)
      call build_regular_grid_1d(grid, x_min=0.0_wp, x_max=1.0_wp, n_cells=4)
      call reset_grid(grid)
      ! node positions must survive a reset
      call check(error, abs(grid%xI(1,1) - 0.00_wp) < 1.0e-14_wp); if (allocated(error)) return
      call check(error, abs(grid%xI(1,5) - 1.00_wp) < 1.0e-14_wp)
      call deallocate_grid(grid)
   end subroutine

   ! ── apply_bcs ─────────────────────────────────────────────────────────────

   subroutine test_apply_bcs_zeroes_fixed_dofs(error)
      type(error_type), allocatable, intent(out) :: error
      type(mpm_grid_t) :: grid
      call allocate_grid(grid, n_cells=4, problem_type=PROB_1D)
      grid%fixed_dof(1,1)   = .true.
      grid%force_I(1,1)     = 9.0_wp
      grid%momentum_I(1,1)  = 7.0_wp
      call apply_bcs(grid)
      call check(error, grid%force_I(1,1)    == 0.0_wp); if (allocated(error)) return
      call check(error, grid%momentum_I(1,1) == 0.0_wp)
      call deallocate_grid(grid)
   end subroutine

   subroutine test_apply_bcs_leaves_free_dofs(error)
      type(error_type), allocatable, intent(out) :: error
      type(mpm_grid_t) :: grid
      call allocate_grid(grid, n_cells=4, problem_type=PROB_1D)
      grid%force_I(1,3)    = 4.0_wp
      grid%momentum_I(1,3) = 6.0_wp
      call apply_bcs(grid)
      call check(error, grid%force_I(1,3)    == 4.0_wp); if (allocated(error)) return
      call check(error, grid%momentum_I(1,3) == 6.0_wp)
      call deallocate_grid(grid)
   end subroutine

   ! ── update_grid ───────────────────────────────────────────────────────────

   subroutine test_update_grid_momentum(error)
      !! momentum_I^{n+1} = momentum_I^n + force_I * dt
      type(error_type), allocatable, intent(out) :: error
      type(mpm_grid_t) :: grid
      call allocate_grid(grid, n_cells=4, problem_type=PROB_1D)
      grid%mass_I(3)       = 2.0_wp
      grid%momentum_I(1,3) = 4.0_wp
      grid%force_I(1,3)    = 1.0_wp
      call update_grid(grid, dt=0.5_wp)
      call check(error, abs(grid%momentum_I(1,3) - 4.5_wp) < 1.0e-14_wp)
      call deallocate_grid(grid)
   end subroutine

   subroutine test_update_grid_vel_and_dv(error)
      !! vel_I = momentum_I^{n+1} / mass_I;  dv_I = vel_I^{n+1} - vel_I^n
      type(error_type), allocatable, intent(out) :: error
      type(mpm_grid_t) :: grid
      call allocate_grid(grid, n_cells=4, problem_type=PROB_1D)
      grid%mass_I(3)       = 2.0_wp
      grid%momentum_I(1,3) = 4.0_wp   ! vel_old = 2.0
      grid%force_I(1,3)    = 1.0_wp
      call update_grid(grid, dt=0.5_wp)
      ! momentum_new = 4.5; vel_new = 2.25; dv = 0.25
      call check(error, abs(grid%vel_I(1,3) - 2.25_wp) < 1.0e-14_wp); if (allocated(error)) return
      call check(error, abs(grid%dv_I(1,3)  - 0.25_wp) < 1.0e-14_wp)
      call deallocate_grid(grid)
   end subroutine

   subroutine test_update_grid_skips_empty(error)
      !! Nodes with mass < MASS_TOL must be skipped (no division by zero).
      type(error_type), allocatable, intent(out) :: error
      type(mpm_grid_t) :: grid
      call allocate_grid(grid, n_cells=4, problem_type=PROB_1D)
      grid%mass_I(2)    = 0.0_wp   ! empty node
      grid%force_I(1,2) = 1.0_wp
      call update_grid(grid, dt=1.0_wp)
      call check(error, grid%vel_I(1,2) == 0.0_wp); if (allocated(error)) return
      call check(error, grid%dv_I(1,2)  == 0.0_wp)
      call deallocate_grid(grid)
   end subroutine

end module test_mpm_grid_suite
