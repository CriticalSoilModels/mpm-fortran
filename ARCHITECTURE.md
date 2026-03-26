# ARCHITECTURE.md

Design notes and rationale for the `mpm-fortran` solver.

---

## Goals and Scope

A Modern Fortran MPM solver targeting home and research-level workstations. Maximum scope: a few CPUs and a few GPUs (primary target: AMD 7900 XTX and similar). Not designed for supercomputers or distributed memory. Both single-phase and multiphase (solid + water + air) problems. General-purpose — not locked to geotechnical applications.

Constitutive models are provided by the `critical-soil-models` sibling library (`../critical-soil-models`). That dependency is deferred until the core solver is validated.

---

## Dimensionality Roadmap

Problems are identified by a `problem_type` constant. Everything — Voigt size, shape function gradient dimension, grid topology — derives from this single value.

| Constant | Description | Voigt size | Spatial dim |
|---|---|---|---|
| `PROB_1D` | 1D bar | 1 | 1 |
| `PROB_2D_PLANE_STRAIN` | 2D plane strain | 4 | 2 |
| `PROB_2D_AXISYM` | 2D axisymmetric | 4 | 2 |
| `PROB_3D` | Full 3D | 6 | 3 |

The Voigt ordering for 3D and 2D plane strain/axisymmetric is `[11, 22, 33, 12, 13, 23]` — identical to `critical-soil-models`. For 1D only the `11` component is active.

---

## Layer Structure

```
┌──────────────────────────────────────────────────────────┐
│  Top-level driver                                        │
│  - time loop, output hooks, restart                      │
├──────────────────────────────────────────────────────────┤
│  mpm_algorithm_t  (abstract)                             │
│  - deferred: step(particle_sets, grid, shapefn, dt)      │
│  - concrete: mpm_usl_t, mpm_usf_t, mpm_musl_t           │
├──────────────────────────────────────────────────────────┤
│  Physics modules  (optional, checked via phase_config)   │
│  - solid mechanics  (always active)                      │
│  - pore water / Biot coupling  (has_liquid)              │
│  - pore air  (has_gas)                                   │
│  - contact mechanics  (has_contact)  [future]            │
├──────────────────────────────────────────────────────────┤
│  mpm_particle_set_t  ×N  +  mpm_grid_t  (shared)        │
├──────────────────────────────────────────────────────────┤
│  mpm_shapefn_t  (abstract)                               │
│  - concrete: linear, B-spline, GIMP                      │
│  - pure function layer: !$acc routine seq  (GPU path)    │
├──────────────────────────────────────────────────────────┤
│  mpm_precision  —  mpm_problem_config                    │
└──────────────────────────────────────────────────────────┘
```

---

## Algorithm Variants

USL, USF, and MUSL share the same underlying subroutines. The only difference is call ordering within the time step. Because the algorithm dispatch happens **once per time step from the CPU**, polymorphism here is essentially free — no GPU impact.

### USL (Update Stress Last) — starting point

```
reset_grid
P2G          (mass, momentum, internal forces)
apply_bc     (nodal boundary conditions)
update_grid  (solve momentum equation)
G2P          (update particle velocity, position)
update_stress
```

### USF (Update Stress First)

```
update_stress   ← moved before P2G
reset_grid
P2G
apply_bc
update_grid
G2P
```

### MUSL (Modified USL)

```
reset_grid
P2G
apply_bc
update_grid
G2P_velocity    ← velocity update only
update_stress
G2P_position    ← position update only
```

Abstract interface:

```fortran
type, abstract :: mpm_algorithm_t
contains
   procedure(step_interface), deferred :: step
end type

abstract interface
   subroutine step_interface(self, sets, grid, shapefn, config, dt)
      import :: mpm_algorithm_t, mpm_particle_set_t, mpm_grid_t
      import :: mpm_shapefn_t, mpm_phase_config_t, wp
      class(mpm_algorithm_t),   intent(inout) :: self
      type(mpm_particle_set_t), intent(inout) :: sets(:)
      type(mpm_grid_t),         intent(inout) :: grid
      class(mpm_shapefn_t),     intent(in)    :: shapefn
      type(mpm_phase_config_t), intent(in)    :: config
      real(wp),                 intent(in)    :: dt
   end subroutine
end interface
```

---

## Particle Sets and Multiphase

The fundamental data unit is `mpm_particle_set_t` — a named SOA array set belonging to one phase. A simulation owns one or more sets. The grid is always shared.

### Phase configuration

```fortran
type :: mpm_phase_config_t
   logical :: has_liquid = .false.   !! Enable pore water phase
   logical :: has_gas    = .false.   !! Enable pore air phase
   logical :: has_contact = .false.  !! Enable contact detection [future]
   integer :: n_sets                 !! Derived: number of particle sets
end type
```

The algorithm checks these flags before calling any phase-specific physics. Unset flags mean those arrays are never allocated and those routines are never called. Zero runtime overhead for single-phase problems.

### Single-point vs two-point

**Single-point mixture:** one particle set. Fluid quantities (pore pressure `uw`, saturation `Sr`) are extra fields within that set. Phases are coupled implicitly through the constitutive model.

**Two-point MPM:** two particle sets — solid and liquid. Each does its own P2G/G2P. Coupling is explicit drag forces resolved at the grid between the two momentum equations. This is the approach used in Anura3D.

Both are supported by the same `mpm_particle_set_t` + `mpm_phase_config_t` infrastructure. The algorithm selects the appropriate loop structure based on `n_sets` and `has_liquid`.

### Particle set type

```fortran
type :: mpm_particle_set_t
   ! --- identity ---
   integer  :: phase        !! PHASE_SOLID, PHASE_LIQUID, PHASE_GAS
   integer  :: body_id      !! for contact detection [future]
   integer  :: n_mp         !! number of material points

   ! --- kinematics (always present) ---
   real(wp), allocatable :: xp(:,:)    !! position  [m]      (ndim, n_mp)
   real(wp), allocatable :: vp(:,:)    !! velocity  [m/s]    (ndim, n_mp)
   real(wp), allocatable :: ap(:,:)    !! accel     [m/s^2]  (ndim, n_mp)
   real(wp), allocatable :: mp(:)      !! mass      [kg]     (n_mp)
   real(wp), allocatable :: vol(:)     !! volume    [m^3]    (n_mp)
   real(wp), allocatable :: rho(:)     !! density   [kg/m^3] (n_mp)

   ! --- solid mechanics (always present) ---
   real(wp), allocatable :: sig(:,:)   !! Cauchy stress [kPa]  (nvoigt, n_mp)
   real(wp), allocatable :: eps(:,:)   !! total strain  [-]    (nvoigt, n_mp)
   real(wp), allocatable :: eps_p(:,:) !! plastic strain [-]   (nvoigt, n_mp)
   real(wp), allocatable :: statev(:,:)!! constitutive state   (nstatev, n_mp)

   ! --- fluid (allocated only when has_liquid or has_gas) ---
   real(wp), allocatable :: uw(:)      !! pore water pressure [kPa] (n_mp)
   real(wp), allocatable :: ua(:)      !! pore air pressure   [kPa] (n_mp)
   real(wp), allocatable :: Sr(:)      !! degree of saturation [-]  (n_mp)
end type
```

SOA layout is mandatory — first index is component, last is particle number. This ensures contiguous memory access per component in GPU kernels (Fortran column-major).

---

## Grid

Regular Cartesian grid to start. Connectivity is trivially generated but stored explicitly in a connectivity array to avoid hardcoding the generation formula — this keeps the door open for unstructured grids later without retrofitting.

```fortran
type :: mpm_grid_t
   ! --- geometry ---
   integer  :: n_nodes              !! total nodes
   integer  :: n_cells              !! total cells
   real(wp) :: dx                   !! cell size [m] (uniform for now)
   real(wp), allocatable :: xI(:,:) !! nodal positions (ndim, n_nodes)

   ! --- connectivity ---
   integer, allocatable :: cell_nodes(:,:) !! (nodes_per_cell, n_cells)

   ! --- nodal fields (reset each step) ---
   real(wp), allocatable :: mass_I(:)       !! nodal mass     [kg]
   real(wp), allocatable :: momentum_I(:,:) !! nodal momentum [kg·m/s] (ndim, n_nodes)
   real(wp), allocatable :: force_I(:,:)    !! nodal force    [N]      (ndim, n_nodes)
   real(wp), allocatable :: vel_I(:,:)      !! nodal velocity [m/s]    (ndim, n_nodes)
   real(wp), allocatable :: acc_I(:,:)      !! nodal accel    [m/s^2]  (ndim, n_nodes)
end type
```

---

## Shape Functions

Shape functions are the performance-critical abstraction. They are called inside the P2G/G2P loops for every particle×node pair. Two paths — same as `critical-soil-models` constitutive models:

- **OOP path** (`class(mpm_shapefn_t)`): CPU, research, debugging
- **Pure function path** (`!$acc routine seq`): GPU kernels, no dynamic dispatch

```fortran
type, abstract :: mpm_shapefn_t
contains
   procedure(eval_interface), deferred :: evaluate
end type

abstract interface
   subroutine eval_interface(self, xp, cell_nodes, xI, N, dN_dx)
      !! Evaluate shape function values and gradients at particle position xp.
      !! N(nodes_per_cell), dN_dx(ndim, nodes_per_cell)
   end subroutine
end interface
```

### Roadmap

| Stage | Type | Notes |
|---|---|---|
| 1 | Linear hat functions | Standard MPM; cell-crossing instability |
| 2 | B-spline (quadratic/cubic) | Smoother; fixes cell-crossing; larger support |
| 3 | GIMP | Tracks particle domain; handles large deformation better |
| ~~4~~ | ~~IGA (Alsardi et al.)~~ | Deferred — revisit when stages 1–3 are stable |

GIMP requires a particle domain field (`lp`, the particle half-length/domain size) in `mpm_particle_set_t`. This field is present but unused by linear and B-spline implementations.

---

## Stress Update and Large Deformation

Even for the first 1D problem, use an objective stress update to avoid retrofitting. The velocity gradient `L` (or equivalently the rate of deformation `D` and spin `W`) must be passed to the stress update, not just the strain increment `deps`.

Jaumann stress rate is the standard starting point:

```
sig_dot_J = sig_dot + sig·W - W·sig
```

The stress update call signature:

```fortran
subroutine update_stress(sig, statev, D, W, dt, params)
   !! D = symmetric part of L (rate of deformation)
   !! W = antisymmetric part of L (spin tensor)
```

This is consistent with the UMAT interface in `critical-soil-models` when that dependency is eventually added.

---

## GPU Design Rules (deferred implementation, active design constraint)

GPU implementation is deferred. These rules constrain design decisions made now so the GPU path is not blocked later.

1. **No polymorphism inside hot loops.** P2G, G2P, and stress update loops must call pure functions, not `class()` procedures. OOP wrappers sit outside the loop.
2. **SOA everywhere.** `(component, n_mp)` array layout. No structs-of-arrays inside GPU kernels.
3. **No allocatables inside GPU regions.** All arrays pre-allocated before `!$acc data` regions.
4. **`!$acc routine seq`** on any pure function called inside a parallel region.
5. **One `!$acc enter data` per run** (or per time step for output). Keep data resident on device between steps.
6. **No NVIDIA-only extensions.** No `attributes(device)`, no `<<<...>>>`. OpenACC only.
7. **Compiler target:** `gfortran -fopenacc -foffload=amdgcn-amdhsa` or `amdflang`. Both must work.

---

## Contact Mechanics (future)

Not implemented. Design impact: `mpm_particle_set_t` carries `body_id` from day one. Contact detection operates on bodies (sets of particle sets), not individual particles or phases. The algorithm checks `config%has_contact` before running any contact detection pass.

---

## Build Roadmap

### Stage 1 — 1D elastic wave (current)

Minimum viable modules:

| Module | Type | Contents |
|---|---|---|
| `mpm_precision` | foundation | `wp`, problem type constants |
| `mpm_particles` | data | `mpm_particle_set_t`, allocate/deallocate |
| `mpm_grid` | data | `mpm_grid_t`, allocate/deallocate, reset |
| `mpm_shapefn_linear` | shape fn | 1D linear hat functions |
| `mpm_algorithm_usl` | algorithm | USL time step for single solid phase |
| `mpm_stress_elastic` | physics | Linear elastic stress update (1D) |
| test suite | tests | Wave speed, energy conservation |

### Stage 2 — 2D plane strain

- Extend `mpm_particle_set_t` to 2D
- `mpm_shapefn_linear_2d` — bilinear quadrilateral shape functions
- `mpm_grid` extended to 2D structured quad mesh
- Validation: 2D wave propagation, footing problem

### Stage 3 — B-spline shape functions

- `mpm_shapefn_bspline` — quadratic B-spline
- Requires extended support radius — affects particle-to-cell search

### Stage 4 — Multiphase (single-point)

- Extend `mpm_particle_set_t` with `uw`, `ua`, `Sr` fields
- `mpm_physics_biot` — Biot coupling module
- Validation: 1D consolidation (Terzaghi)

### Stage 5 — Two-point MPM

- Second particle set for fluid phase
- Coupled grid equations with drag force
- Validation: 1D consolidation comparison with single-point

### Stage 6 — GIMP

- `mpm_shapefn_gimp`
- `lp` (particle domain) field in `mpm_particle_set_t`

### Stage 7 — GPU (OpenACC)

- Pure function extraction from all hot loops
- `!$acc routine seq` annotations
- `!$acc data` regions in time loop
- Validation: performance comparison CPU vs AMD GPU

---

## Validation Problems

| Problem | Tests | Stage |
|---|---|---|
| 1D elastic wave | wave speed = `sqrt(E/rho)`, energy conservation | 1 |
| 2D wave propagation | isotropy, dispersion | 2 |
| Elastic footing | load-displacement, comparison to FEM | 2 |
| 1D Terzaghi consolidation | pore pressure dissipation vs analytical | 4 |
| Slope stability | failure mechanism, run-out distance | post stage 4 |
