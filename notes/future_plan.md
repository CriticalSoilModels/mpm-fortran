# Future Plans

Deferred features and design notes. These are not on the immediate roadmap but
the architecture should not actively prevent them.

---

## Axisymmetric Problem Types

`PROB_2D_AXISYM` and `PROB_3D_AXISYM` are deferred. When added:

- Restore `PROB_2D_AXISYM = 4` (and `PROB_3D_AXISYM = 5` if needed) in `mpm_precision`
- Add `case (PROB_2D_AXISYM)` back to `n_voigt` (returns 4) and `n_dims` (returns 2)
- Shape functions require the Jacobian in cylindrical coordinates
- The weak form gains a `1/r` weighting term — this affects P2G and grid force assembly
- Stress components map as: V_11=σ_rr, V_22=σ_zz, V_33=σ_θθ, V_12=σ_rz

---

## Multiphase (Solid + Water + Air)

### Motivation

Geotechnical problems often involve partially or fully saturated soils where
pore water and pore air pressures play a critical role (consolidation, slope
failure, liquefaction). A three-phase formulation (solid skeleton + pore water
+ pore air) is required for these cases.

### Design sketch

#### Phase configuration

```fortran
! In a future mpm_phase_config module:
type :: mpm_phase_config_t
   logical :: has_liquid  = .false.  ! enable pore water phase
   logical :: has_gas     = .false.  ! enable pore air phase
   logical :: has_contact = .false.  ! enable contact detection
end type
```

#### Extended particle set

Add fluid fields to `mpm_particle_set_t`, allocated only when the phase config
requests them (zero overhead for single-phase problems):

```fortran
! Additional fields in mpm_particle_set_t:
integer :: phase = PHASE_SOLID        ! PHASE_SOLID / PHASE_LIQUID / PHASE_GAS

! Fluid fields — allocate only when config%has_liquid / config%has_gas:
real(wp), allocatable :: uw(:)  ! pore water pressure [kPa] (n_mp)
real(wp), allocatable :: ua(:)  ! pore air pressure   [kPa] (n_mp)
real(wp), allocatable :: Sr(:)  ! degree of saturation [-]  (n_mp); default 1.0
```

#### Phase constants

```fortran
integer, parameter :: PHASE_SOLID  = 1
integer, parameter :: PHASE_LIQUID = 2
integer, parameter :: PHASE_GAS    = 3
```

#### Extended allocate_particle_set

```fortran
subroutine allocate_particle_set(set, n_mp, problem_type, nstatev, config)
   type(mpm_particle_set_t),          intent(out) :: set
   integer,                           intent(in)  :: n_mp
   integer,                           intent(in)  :: problem_type
   integer,                           intent(in)  :: nstatev
   type(mpm_phase_config_t), optional, intent(in) :: config

   ! ... allocate solid fields as now ...

   if (present(config)) then
      if (config%has_liquid) then
         allocate(set%uw(n_mp), source=0.0_wp)
         allocate(set%Sr(n_mp), source=1.0_wp)  ! fully saturated default
      end if
      if (config%has_gas) then
         allocate(set%ua(n_mp), source=0.0_wp)
      end if
   end if
end subroutine
```

### Single-point vs two-point MPM

**Single-point mixture:** one particle set carries averaged properties. Fluid
quantities (`uw`, `ua`, `Sr`) are extra fields on the solid particle set.
Coupling is implicit in the constitutive model.

**Two-point MPM:** separate particle sets for solid and fluid phases. Each
does its own P2G/G2P pass. Coupling is explicit drag forces at the grid.
This is how Anura3D implements saturated soil.

Both approaches share the same `mpm_particle_set_t` infrastructure. The
algorithm selects the appropriate loop structure based on `n_sets` and
`config%has_liquid`.

### Physics modules required

- `mpm_physics_biot` — Biot coupling module; drag forces between solid and
  fluid momentum equations; sits between grid update and G2P
- Geotechnical-specific: van Genuchten retention curve, hydraulic
  conductivity, Biot coefficient, air compressibility — all in this module,
  not in the core solver

### Validation problems

| Problem | Tests |
|---|---|
| 1D Terzaghi consolidation | pore pressure dissipation vs analytical |
| 2D consolidation | settlement and pore pressure distribution |
| Slope failure (unsaturated) | failure mechanism and run-out |

---

## Contact Mechanics

`mpm_particle_set_t` already carries `body_id` for this purpose. Contact
detection operates on bodies (groups of particle sets), not individual
particles.

```fortran
! Sketch of contact algorithm hook in the time loop:
if (config%has_contact) then
   call detect_contact(sets, grid, contact_pairs)
   call apply_contact_forces(grid, contact_pairs)
end if
```

---

## APIC Update Scheme

APIC (Affine Particle-In-Cell, Jiang et al. 2015) requires an affine velocity
matrix `C_p` per particle. Unlike FLIP/PIC, it also modifies the P2G step.

```fortran
! Additional field in mpm_particle_set_t (allocated when SCHEME_APIC active):
real(wp), allocatable :: Cp(:,:)  ! affine velocity matrix, flattened SOA
                                  ! shape: (ndim*ndim, n_mp)

! Modified P2G with APIC:
! p_i = Σ_p N_ip * m_p * (v_p + C_p * (x_i - x_p))

! G2P with APIC:
! v_p^{n+1} = Σ_i N_ip * v_i^{n+1}
! C_p^{n+1} = Σ_i N_ip * v_i^{n+1} ⊗ (x_i - x_p)
```

---

## IGA Shape Functions (Alsardi et al., Virginia Tech)

Isogeometric analysis shape functions for MPM. Deferred until B-spline and
GIMP implementations are stable and the original paper is well understood.
Revisit after stage 3 of the build roadmap.

---

## GPU (OpenACC)

See `ARCHITECTURE.md` — GPU design rules are active constraints now even
though implementation is deferred. The pure function / `!$acc routine seq`
pattern must be preserved throughout development so the GPU path is not
blocked when the time comes.
