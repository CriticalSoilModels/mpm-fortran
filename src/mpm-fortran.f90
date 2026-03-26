module mpm_fortran
   !! Top-level re-export module for the mpm-fortran library.
   !!
   !! Consumers can either `use mpm_fortran` to get everything, or import
   !! individual modules directly for finer control.
   use mpm_precision,      only: wp
   use mpm_problem_config, only: PROB_1D, PROB_2D_PLANE_STRAIN, PROB_3D, &
                                  V_11, V_22, V_33, V_12, V_13, V_23,    &
                                  n_voigt, n_dims, MASS_TOL
   use mpm_particles,      only: mpm_particle_set_t, &
                                  allocate_particle_set, deallocate_particle_set
   use mpm_update_config,  only: mpm_update_config_t, &
                                  SCHEME_FLIP_PIC, SCHEME_APIC, SCHEME_XPIC, &
                                  is_valid_config
   use mpm_shapefn_linear, only: locate_cell_1d, eval_linear_1d
   use mpm_grid,           only: mpm_grid_t, nodes_per_cell, &
                                  allocate_grid, deallocate_grid, &
                                  build_regular_grid_1d, &
                                  reset_grid, apply_bcs, update_grid
   implicit none
   private

   public :: wp
   public :: PROB_1D, PROB_2D_PLANE_STRAIN, PROB_3D
   public :: V_11, V_22, V_33, V_12, V_13, V_23
   public :: n_voigt, n_dims, MASS_TOL
   public :: mpm_particle_set_t
   public :: allocate_particle_set, deallocate_particle_set
   public :: mpm_update_config_t
   public :: SCHEME_FLIP_PIC, SCHEME_APIC, SCHEME_XPIC
   public :: is_valid_config
   public :: mpm_grid_t, nodes_per_cell
   public :: allocate_grid, deallocate_grid
   public :: build_regular_grid_1d
   public :: reset_grid, apply_bcs, update_grid
   public :: locate_cell_1d, eval_linear_1d

end module mpm_fortran
