# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

A Modern Fortran Material Point Method (MPM) solver targeting CPU and GPU execution. The primary GPU target is AMD (via OpenACC/HIP), but GPU portability is a design goal — no vendor-specific intrinsics in core kernels. Constitutive models are provided by the `critical-soil-models` library (sibling repo at `../critical-soil-models`).

### Sibling Repositories

| Repo | Path | Role |
|------|------|------|
| `critical-soil-models` | `../critical-soil-models` | Constitutive model library; provides UMAT subroutines and `csm_model_t` OOP interface |
| `Incremental_Driver` | `../Incremental_Driver` | Single-element test harness for constitutive models |
| `fumat` | `../fumat` | Calibration, analysis, and plotting of element test results |

The MPM solver calls into `critical-soil-models` for all stress updates. New constitutive models should be developed and tested there before being used here.

## Build and Test Commands

```bash
# Build the library
fpm build

# Run all tests
fpm test

# Run a specific test by name (glob supported: ? = any char, * = any string)
fpm test --target test_mpm_grid

# Pass arguments to the test executable
fpm test -- -x 10 -y 20

# List available test targets without running them
fpm test --list

# Run with release optimisations
fpm test --profile release

# Run the main application
fpm run

# Run a specific executable
fpm run --target mpm_solver

# Build with OpenACC for AMD GPU (ROCm/HIP backend)
fpm build --flag "-fopenacc -foffload=amdgcn-amdhsa"

# Build with a different compiler
fpm build --compiler amdflang

# Clean build artefacts (keeps dependencies)
fpm clean

# Clean everything including dependencies
fpm clean --all

# Update dependencies
fpm update

# Generate documentation (requires ford)
ford fpm.toml
```

### Key fpm flags

| Flag | Description |
|------|-------------|
| `--profile debug\|release` | Compilation profile; default is `debug` |
| `--flag FFLAGS` | Extra Fortran compiler flags (appended to profile defaults) |
| `--compiler NAME` | Override Fortran compiler (default: `gfortran`; env: `FPM_FC`) |
| `--no-prune` | Disable tree-shaking of unused module dependencies |
| `--verbose` | Display additional build/run information |
| `--build-dir DIR` | Override build output directory (env: `FPM_BUILD_DIR`) |

### Environment variables

| Variable | Purpose |
|----------|---------|
| `FPM_FC` | Default Fortran compiler |
| `FPM_FFLAGS` | Default Fortran compiler flags |
| `FPM_BUILD_DIR` | Default build directory |

### Environment Setup
```bash
conda env create -f environment.yml
conda activate mpm
```
Requires: gfortran (or nvfortran/amdflang for GPU), fpm, fortls, ford.

## Test-Driven Development Cycle

This project follows a strict TDD cycle. Every module gets its tests written **before** the implementation.

### Cycle

1. **Write a failing test** — add a test case in `test/test_<module>_suite.f90` that expresses the desired behaviour
2. **Confirm it fails** — `fpm test -- test_<suite_name>`; if it compiles and passes without implementation, the test is wrong
3. **Write minimum code** — implement only enough to make the test pass; no extra features
4. **Confirm it passes** — `fpm test -- test_<suite_name>`
5. **Refactor** — clean up, then confirm tests still pass
6. **Commit** — one commit per completed red→green→refactor cycle (see Git Workflow below)

### Test file pattern (test-drive framework)

```fortran
module test_mpm_grid_suite
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mpm_precision, only: wp, PROB_1D
   use mpm_grid,      only: mpm_grid_t, allocate_grid
   implicit none
   private

   public :: collect_mpm_grid_suite

contains

   subroutine collect_mpm_grid_suite(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("grid_allocates_correct_size", test_grid_allocates_correct_size), &
         new_unittest("grid_reset_zeroes_nodal_mass", test_grid_reset_zeroes_nodal_mass)  &
      ]
   end subroutine

   subroutine test_grid_allocates_correct_size(error)
      type(error_type), allocatable, intent(out) :: error
      ! ...
      call check(error, size(grid%xI, 2) == n_nodes); if (allocated(error)) return
   end subroutine

end module
```

The main test runner (`test/main.f90`) collects all suites:

```fortran
program test_runner
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mpm_precision_suite, only: collect_mpm_precision_suite
   ! ... other suites
   implicit none
   ! register and run all suites
end program
```

### What to test

- **Allocation:** arrays have correct shapes after `allocate_*`
- **Initialisation:** zero-initialised fields are zero; non-zero defaults (e.g. `Sr = 1.0`) are correct
- **Pure functions:** `n_voigt`, `n_dims`, shape function values/gradients against analytical results
- **Conservation:** energy, mass, and momentum conserved to machine precision for elastic problems
- **Regression:** known solutions (wave speed, consolidation) match analytical within tolerance

### What not to test

- Fortran language behaviour (allocation, assignment)
- Internal implementation details that may change during refactor

---

## Git Workflow

Commit after each completed TDD cycle and after each major module is validated. Never commit broken code — all tests must pass before committing.

```bash
# Stage specific files (never git add -A)
git add src/mpm_grid.f90 test/test_mpm_grid_suite.f90

# Commit
git commit -m "add mpm_grid_t allocation and reset"

# Check before pushing
fpm test
```

### Commit message conventions

- `add <module>` — new module or type
- `add <feature> to <module>` — new capability in existing module
- `fix <description>` — bug fix
- `test <module>` — test suite only (no implementation change)
- `refactor <description>` — behaviour-preserving cleanup

One logical change per commit. If a commit message needs "and", it should probably be two commits.

## Architecture

See `ARCHITECTURE.md` for full design notes. Summary of the intended layered structure:

```
┌──────────────────────────────────────────────────────────┐
│  Top-level solver loop                                   │
│  - time stepping, output, restart                        │
├──────────────────────────────────────────────────────────┤
│  MPM algorithm layer   (mpm_algorithm_t)                 │
│  - USL / USF / MUSL variants                             │
│  - particle-to-grid (P2G) mapping                        │
│  - grid momentum update + BCs                            │
│  - grid-to-particle (G2P) mapping                        │
│  - stress update (calls constitutive layer below)        │
├──────────────────────────────────────────────────────────┤
│  Constitutive layer                                      │
│  - calls csm_model_t (OOP, CPU research path)            │
│  - OR calls umat_* directly (flat, GPU path)             │
│  - pure functions: !$acc routine seq                     │
├──────────────────────────────────────────────────────────┤
│  Grid  (mpm_grid_t)          Particles (mpm_particles_t) │
│  - nodal mass, momentum      - position, velocity        │
│  - nodal forces              - stress, strain            │
│  - boundary conditions       - internal state vars       │
├──────────────────────────────────────────────────────────┤
│  Shape functions  (mpm_shapefn_t)                        │
│  - GIMP, B-spline, linear                                │
│  - values + gradients at material points                 │
├──────────────────────────────────────────────────────────┤
│  Shared infrastructure                                   │
│  - precision (mpm_precision)                             │
│  - Voigt convention (mirrors critical-soil-models)       │
│  - I/O, restart, VTK output                              │
└──────────────────────────────────────────────────────────┘
```

### Key Modules (planned)

- `mpm_precision` — working precision `wp = real64`; single source of truth
- `mpm_particles` — `mpm_particles_t` SOA layout for GPU coalescing
- `mpm_grid` — `mpm_grid_t` background mesh, nodal arrays
- `mpm_shapefn` — shape function evaluation and gradients
- `mpm_algorithm` — abstract `mpm_algorithm_t` with concrete USL/USF variants
- `mpm_bc` — boundary conditions
- `mpm_io` — VTK/HDF5 output, restart files
- `mpm_solver` — top-level driver

### GPU Design Rules

All hot-path kernels (P2G, G2P, stress update loops) must be callable from OpenACC device regions:

- **No polymorphism in GPU kernels** — use flat pure functions, not `class()` dispatch
- **SOA (Structure of Arrays)** for particle storage, not AOS — enables coalesced access in Fortran column-major order
- **`!$acc routine seq`** on any pure function called inside a parallel region
- **Voigt convention** internal ordering `[11,22,33,12,13,23]` — same as `critical-soil-models`; translate only at I/O boundaries
- **No allocatables inside GPU regions** — pre-allocate everything before entering `!$acc data` region
- Copy data to GPU once per time step (or keep resident across steps with `!$acc enter data`); copy back only for output

### Constitutive Model Integration

The stress update loop calls into `critical-soil-models`:

- **CPU/research path:** Use `csm_model_t` polymorphic interface with `euler_substep` integrator
- **GPU path:** Call pure functions from `*_functions.f90` directly, marked `!$acc routine seq`; avoid dynamic dispatch inside `!$acc parallel` regions

---

# Modern Fortran Style Guide

A practical style guide for modern Fortran (2008+), with emphasis on GPU-accelerated scientific computing using OpenACC and CUDA Fortran.

## Naming Conventions

### Files

- Module files: `module_name.f90` (lowercase with underscores)
- Preprocessed files: `module_name.F90` (uppercase extension for files needing preprocessing)
- Test files: `test_module_name.f90`
- One module per file (submodules may be separate)

### Modules

Use the `mpm_` prefix to namespace all modules in this repo:

```fortran
! Good — prefixed with project abbreviation
module mpm_particles
module mpm_grid
module mpm_shapefn

! Bad
module ShallowWaterSolver   ! CamelCase
module particles            ! No prefix, collision risk
```

### Derived Types

All derived types use the `_t` suffix with snake_case:

```fortran
type :: mpm_particles_t
type :: mpm_grid_t
type :: mpm_shapefn_t
type :: mpm_algorithm_t
```

### Variable Naming for Mathematical Quantities

Key rules (consistent with `critical-soil-models` conventions):
- Increments: `d"var"` — `dsig`, `deps`, `dv`, `dt`
- Partial derivatives: `d"num"_by_d"denom"` — `dN_by_dx`, `dN_by_dy`
- Stress vector: `sig(6)`, strain increment: `deps(6)`, elastic stiffness: `stiff_e(6,6)`
- Shape function values: `N(n_nodes)`, gradients: `dN_dx(n_nodes,3)`
- Particle arrays: `xp(3,n_mp)` — first index spatial, second index particle (column-major)

### Variables and Procedures

Snake_case with descriptive names:

```fortran
integer  :: n_mp          ! number of material points
integer  :: n_nodes       ! number of grid nodes
real(wp) :: dt            ! time step [s]
real(wp) :: t_end         ! end time [s]
subroutine map_p2g(particles, grid, shapefn)
subroutine map_g2p(grid, particles, shapefn)
```

**Exception:** Single-letter variables are acceptable for loop indices (`i`, `j`, `k`, `p` for particle index) and mathematical formulas matching published notation.

### Function and Subroutine Naming

Verb + noun patterns:

```fortran
subroutine map_p2g(...)             ! particle-to-grid
subroutine map_g2p(...)             ! grid-to-particle
subroutine update_stress(...)
subroutine apply_boundary_conditions(...)
function compute_stable_dt(...) result(dt)

! Logical returns use is_/has_/can_ prefixes
function is_active(particle_id) result(active)
function has_converged(residual, tol) result(converged)
```

### Units in Comments

Document physical units in variable declarations:

```fortran
type :: mpm_particles_t
   real(wp), allocatable :: xp(:,:)   !! Position [m], shape (3, n_mp)
   real(wp), allocatable :: vp(:,:)   !! Velocity [m/s], shape (3, n_mp)
   real(wp), allocatable :: mp(:)     !! Mass [kg], shape (n_mp)
   real(wp), allocatable :: sig(:,:)  !! Cauchy stress [kPa], shape (6, n_mp) Voigt [11,22,33,12,13,23]
   real(wp), allocatable :: eps_p(:,:)!! Plastic strain [-], shape (6, n_mp)
end type

real(wp), parameter :: GRAVITY = 9.81_wp   ! [m/s^2]
```

### Constants

UPPERCASE with underscores:

```fortran
real(wp), parameter :: GRAVITY = 9.81_wp
real(wp), parameter :: DT_COURANT_FACTOR = 0.5_wp
integer,  parameter :: MAX_ITER = 1000
integer,  parameter :: VOIGT_SIZE = 6
```

---

## Required Practices

### Use Statements with `only` Clause

```fortran
! Good
use iso_fortran_env, only: real64, int32
use mpm_precision,   only: wp
use mpm_particles,   only: mpm_particles_t

! Bad — pollutes namespace, hides dependencies
use iso_fortran_env
```

### Implicit None

Always in modules and programs:

```fortran
module mpm_particles
   implicit none
   private
end module
```

### Intent Declarations

Always declare intent for all procedure arguments:

```fortran
subroutine map_p2g(particles, grid, shapefn)
   type(mpm_particles_t), intent(in)    :: particles
   type(mpm_grid_t),      intent(inout) :: grid
   type(mpm_shapefn_t),   intent(in)    :: shapefn
```

### Functions Should Have No Side Effects

```fortran
! Good — pure computation
pure function kinetic_energy(mp, vp) result(ke)
   real(wp), intent(in) :: mp(:), vp(:,:)
   real(wp) :: ke

! Bad — function mutates state (use a subroutine instead)
function update_and_return_energy(particles) result(energy)
   type(mpm_particles_t), intent(inout) :: particles    ! Side effect hidden in function
```

### Private by Default

```fortran
module mpm_grid
   implicit none
   private

   public :: mpm_grid_t
   public :: initialize_grid
   public :: reset_grid
```

### Limit Procedure Arguments

Public procedures should have **6 or fewer arguments**. Group related arguments into derived types:

```fortran
! Bad
subroutine run_mpm(nx, ny, nz, dx, dy, dz, dt, t_end, n_mp, ...)

! Good
subroutine run_mpm(grid, particles, config, output)
```

**Performance note:** For performance-critical GPU kernels (P2G/G2P loops), explicit scalar/array arguments are acceptable over derived type indirection when needed for `!$acc routine seq` compatibility.

---

## Forbidden Practices

- **No `goto`** — use structured control flow
- **No arithmetic IF** — use `if-then-else` or `select case`
- **No COMMON blocks** — use module variables or derived types
- **No EQUIVALENCE**
- **No fixed-form source** — all files must be `.f90` / `.F90`
- **No assumed-size arrays** (`arr(*)`) — use assumed-shape (`arr(:)`)
- **No `external` statements** — use modules for explicit interfaces
- **No implicit save** — avoid module-level variables with initialisation; if needed use explicit `save`

---

## Recommended Practices

### Working Precision

```fortran
module mpm_precision
   use iso_fortran_env, only: real64
   implicit none
   integer, parameter :: wp = real64
end module

! Never use literal kind numbers
real(8) :: x           ! Bad — non-portable
real(wp) :: x          ! Good
```

### Prefer Allocatable Over Pointer

```fortran
! Good
real(wp), allocatable :: xp(:,:)

! Bad — manual cleanup required
real(wp), pointer :: xp(:,:) => null()
```

Use pointers only when aliasing, linked structures, or polymorphic returns are needed.

### SOA Layout for GPU Particles

Store particle data as Structure of Arrays (SOA), not Array of Structures (AOS). Fortran is column-major so the particle index should be the **last** (slowest-varying) index, with the component index first:

```fortran
! Good — SOA; inner loop over particles hits contiguous memory
real(wp), allocatable :: xp(:,:)   ! (3, n_mp)
real(wp), allocatable :: sig(:,:)  ! (6, n_mp)

!$acc parallel loop
do p = 1, n_mp
   sig(1,p) = ...   ! [11] component
end do

! Bad — AOS; non-contiguous per component in GPU thread
type :: particle_t
   real(wp) :: x(3)
   real(wp) :: sig(6)
end type
type(particle_t), allocatable :: particles(:)   ! (n_mp)
```

### Pure and Elemental Procedures

```fortran
pure function shape_fn_value(xi, eta) result(N)
   real(wp), intent(in) :: xi, eta
   real(wp) :: N(4)
   ! ...
end function

elemental function degrees_to_radians(deg) result(rad)
   real(wp), intent(in) :: deg
   real(wp) :: rad
   rad = deg * PI / 180.0_wp
end function
```

### Associate for Readability

```fortran
associate(sig => particles%sig(:,p), &
          eps => particles%eps(:,p),  &
          m   => particles%mp(p))
   ke = 0.5_wp * m * dot_product(vp, vp)
end associate
```

**Compiler caveat:** `associate` support varies. Flang tends to be most robust. Test in performance-critical code and fall back to explicit temporaries if issues arise with gfortran/nvfortran/amdflang.

### No Magic Numbers

```fortran
real(wp), parameter :: DT_COURANT_FACTOR = 0.5_wp

! Bad
dt = 0.5_wp * dx / c_wave

! Good
dt = DT_COURANT_FACTOR * dx / c_wave
```

### Avoid Deep Nesting

Maximum 3-4 levels. Use early `cycle` and `return`:

```fortran
! Good — flat structure
do p = 1, n_mp
   if (.not. active(p)) cycle
   if (mp(p) < MP_MASS_TOL) cycle
   ! work here
end do
```

### Labeled Loops

When using `cycle`/`exit` with nested loops:

```fortran
outer: do p = 1, n_mp
   inner: do n = 1, n_nodes
      if (found(p,n)) exit outer
      if (skip(n)) cycle inner
   end do inner
end do outer
```

### Documentation (FORD compatible)

```fortran
type :: mpm_particles_t
   !! Material point (particle) storage in SOA layout.
   !!
   !! All arrays have shape (component, n_mp). Stress uses Voigt
   !! ordering [11,22,33,12,13,23] consistent with critical-soil-models.

   real(wp), allocatable :: xp(:,:)  !! Position [m]
   real(wp), allocatable :: vp(:,:)  !! Velocity [m/s]
   real(wp), allocatable :: mp(:)    !! Mass [kg]
   real(wp), allocatable :: sig(:,:) !! Cauchy stress [kPa] (6, n_mp)
end type
```

---

## GPU Programming

### OpenACC Basics

```fortran
!$acc parallel loop collapse(2) default(present)
do p = 1, n_mp
   do i = 1, 3
      xp(i,p) = xp(i,p) + dt * vp(i,p)
   end do
end do
```

### Data Management Strategy

Copy to GPU once per run (or keep resident), copy back only for output:

```fortran
!$acc enter data copyin(particles%xp, particles%vp, particles%sig, &
!$acc&                  grid%mass, grid%momentum)

do while (t < t_end)
   call map_p2g(particles, grid, shapefn)
   call update_grid(grid, dt)
   call map_g2p(grid, particles, shapefn, dt)
   call update_stress(particles)
end do

!$acc exit data copyout(particles%xp, particles%vp, particles%sig)
```

### Memory Coalescing

In Fortran (column-major), the inner loop index must be the **first** array index. For SOA particle arrays with shape `(component, n_mp)`, loop over particles in the outer loop and components in the inner loop:

```fortran
! Good — contiguous access per thread warp
!$acc parallel loop collapse(2)
do p = 1, n_mp
   do i = 1, 6
      sig(i,p) = sig(i,p) + dsig(i,p)
   end do
end do
```

### Routines for GPU-Callable Functions

Any pure function called inside an `!$acc parallel` region must be declared with `!$acc routine seq`:

```fortran
!$acc routine seq
pure function voigt_trace(sig) result(tr)
   real(wp), intent(in) :: sig(6)
   real(wp) :: tr
   tr = sig(1) + sig(2) + sig(3)
end function
```

This is identical to the pattern used in `critical-soil-models/src/models/*/le_functions.f90`.

### Do Concurrent (Use with Caution)

```fortran
do concurrent (p = 1:n_mp)
   mp(p) = rho0 * Vp(p)
end do
```

No `exit`, `cycle`, `return`, or `goto` inside. Compiler support varies — prefer regular `do` when in doubt. `do concurrent` with `reduce` is Fortran 2018 and not universally supported.

### AMD GPU Notes

- Primary compiler path: `gfortran -fopenacc` (GCC offload to AMDGPU) or `amdflang` (LLVM-based)
- OpenACC with GCC offload uses `libgomp`; link with `-fopenacc -foffload=amdgcn-amdhsa`
- Avoid `!$acc atomic` in inner loops — use reduction clauses where possible
- `nvfortran`/`pgfortran` CUDA Fortran extensions (`attributes(device)`, `<<<...>>>`) are NVIDIA-only — do not use

---

## Voigt Convention

Internal ordering is `[11, 22, 33, 12, 13, 23]` — identical to `critical-soil-models`. This must be consistent across all modules. Translation happens only at I/O boundaries (VTK output, restart files, UMAT wrappers).

```fortran
! Voigt indices: 1=11, 2=22, 3=33, 4=12, 5=13, 6=23
integer, parameter :: V_11 = 1
integer, parameter :: V_22 = 2
integer, parameter :: V_33 = 3
integer, parameter :: V_12 = 4
integer, parameter :: V_13 = 5
integer, parameter :: V_23 = 6
```

---

## Common AI/LLM Mistakes in Fortran

### Pi is Not a Built-in Constant

```fortran
real(wp), parameter :: PI = 4.0_wp * atan(1.0_wp)
```

### `random_number` is a Subroutine

```fortran
call random_number(x)   ! correct
x = random_number()     ! wrong
```

### No I/O in `pure` Procedures

```fortran
pure function compute(x) result(y)
   real(wp), intent(in) :: x
   real(wp) :: y
   ! print *, x   ! compiler error — I/O is a side effect
   y = x**2
end function
```

### Declarations Before Executable Code

```fortran
subroutine foo()
   real(wp) :: x   ! declarations first
   x = 1.0_wp      ! then executable statements
end subroutine
```

### Array Constructors

```fortran
integer :: arr(3) = [1, 2, 3]     ! modern (Fortran 2003+)
integer :: arr(3) = (/1, 2, 3/)   ! old style, still valid
integer :: arr(3) = (1, 2, 3)     ! wrong — compiler error
```

### Allocatable Components in GPU Regions

```fortran
! Allocatable components of derived types cannot be used inside !$acc parallel
! Pass raw arrays instead of types with allocatable components to GPU kernels
!$acc parallel loop
do p = 1, n_mp
   ! particles%xp(1,p)  ! BAD if xp is allocatable inside an acc region with a derived type arg
   xp(1,p) = ...         ! GOOD — pass raw array to kernel
end do
```

---

## Summary Table

| Category | Do | Don't |
|---|---|---|
| Module prefix | `mpm_particles`, `mpm_grid` | `particles`, `GridModule` |
| Types | `mpm_particles_t`, `mpm_grid_t` | `Particles`, `TGrid` |
| Variables | `n_mp`, `dt`, `xp` | `nMP`, `DeltaT`, `X_Particles` |
| Constants | `GRAVITY`, `DT_COURANT_FACTOR` | `g`, `courant` |
| Imports | `use mpm_particles, only: mpm_particles_t` | `use mpm_particles` |
| Arrays | `arr(:)` | `arr(*)` |
| Memory | `allocatable` | `pointer` (unless needed) |
| Particle layout | SOA: `sig(6, n_mp)` | AOS: `particles(n_mp)%sig(6)` |
| GPU functions | `!$acc routine seq` on pure fn | `class()` dispatch in GPU region |
| Control | `if/select/do` | `goto` |
| Voigt order | `[11,22,33,12,13,23]` | Any other ordering in core code |
