# Next Steps

## Where we are

44 tests passing across 5 modules:

| Module | File | What it does |
|--------|------|-------------|
| `mpm_precision` | `src/mpm_precision.f90` | Working precision `wp = real64` |
| `mpm_problem_config` | `src/mpm_problem_config.f90` | Problem type constants, `n_voigt`, `n_dims`, `MASS_TOL` |
| `mpm_particles` | `src/mpm_particles.f90` | SOA particle storage, allocate/deallocate |
| `mpm_update_config` | `src/mpm_update_config.f90` | Update scheme constants (FLIP/PIC/APIC/XPIC), `is_valid_config` |
| `mpm_grid` | `src/mpm_grid.f90` | Grid storage, `build_regular_grid_1d`, `reset_grid`, `apply_bcs`, `update_grid` |
| `mpm_shapefn_linear` | `src/mpm_shapefn_linear.f90` | `locate_cell_1d`, `eval_linear_1d` — pure, `!$acc routine seq` |

---

## Immediate next step: `mpm_algorithm_usl`

This is where all the modules above come together for the first time into a
full USL (Update Stress Last) time step.

### The USL step sequence

```
reset_grid
P2G:       scatter mass, momentum, internal forces from particles to grid
apply_bcs  (zero fixed DOFs before momentum solve)
update_grid (momentum_I += force_I * dt; compute vel_I, dv_I)
apply_bcs  (zero fixed DOFs after momentum solve)
G2P:       gather velocity/velocity-increment back to particles
           update particle positions, volumes, densities
Stress update: compute strain increment, Jaumann correction, elastic update
```

### Design decisions already made

- **Flat subroutines**, not a type with deferred procedures — the type-based
  OOP layer (`mpm_algorithm_t`) can come later; a flat `step_usl` subroutine
  is enough to validate the physics
- **No allocations inside the step** — all arrays pre-allocated; the step just
  reads and writes them
- **GPU-ready inner loops** — P2G and G2P loops should be written with
  `!$acc parallel loop` annotations from the start (they compile and run
  correctly without OpenACC, and are already correct for the GPU path)
- **Linear elastic stress update** for the first implementation — no
  constitutive model library dependency yet; just `sig += E * deps` for 1D

### What to test

- **P2G:** mass conservation (sum of nodal masses == sum of particle masses)
- **P2G:** momentum conservation (sum of nodal momenta == sum of particle momenta)
- **G2P:** position update with uniform grid velocity (all particles advance by `v*dt`)
- **G2P:** FLIP velocity update (dv correctly interpolated)
- **Stress update:** elastic wave — stress increment proportional to strain increment
- **Full step:** energy conservation for a single elastic step with no body force

### Validation problem (after the step is working)

1D elastic wave in a bar:
- Fixed left end, free right end
- Initial velocity impulse on leftmost particles
- Wave speed `c = sqrt(E/rho)` should match analytical
- This is the first end-to-end test

---

## After `mpm_algorithm_usl`

In rough priority order:

1. **`mpm_io`** — VTK or CSV output so results can be visualised
2. **Example program** (`app/main.f90` or `examples/elastic_wave_1d.f90`) —
   first demonstration of the Fortran-as-API design
3. **`mpm_shapefn` abstract layer** — thin OOP wrapper over the flat shape
   function modules; needed when the algorithm layer needs to be
   shape-function-agnostic
4. **2D extension** — `build_regular_grid_2d`, `eval_linear_2d`, update
   `mpm_algorithm_usl` for 2D plane strain
5. **Constitutive model integration** — link to `critical-soil-models` for
   the stress update; replace the linear elastic stub
6. **OpenACC GPU kernels** — annotate P2G and G2P loops; test on AMD 7900 XTX

---

## Getting started on a new machine

```bash
# Clone the repo
git clone https://github.com/CriticalSoilModels/mpm-fortran.git
cd mpm-fortran

# Set up the conda environment
conda env create -f environment.yml
conda activate mpm

# Build and test
fpm test
```

All 44 tests should pass on a clean clone.
