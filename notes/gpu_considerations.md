# GPU Considerations

Design rules and gotchas for OpenACC GPU execution. Implementation is deferred
but these constraints are active — code written now must not block the GPU path
later. Primary target: AMD (ROCm/HIP) via `gfortran -fopenacc` or `amdflang`.

---

## Memory Layout

### SOA is correct — 2D and flat 1D are identical

`sig(nv, n_mp)` in Fortran column-major is a flat 1D array in memory:

```
sig(1,1), sig(2,1), ..., sig(6,1), sig(1,2), sig(2,2), ..., sig(6,2), ...
```

A flat `sig(nv*n_mp)` indexed as `sig(i + (p-1)*nv)` is byte-for-byte
identical. The 2D syntax is preferred — it is not a performance liability.

### Coalescing with collapse(2)

A pure particle loop has stride-`nv` access when warp threads each read the
same Voigt component of adjacent particles:

```fortran
! stride-nv — partially coalesced (nv=6 in 3D → ~17% efficiency)
!$acc parallel loop
do p = 1, n_mp
   sig(1, p) = sig(1, p) + dsig(1, p)
end do
```

Use `collapse(2)` to map threads to `(i, p)` pairs. Fortran column-major then
gives stride-1 access across the warp:

```fortran
! stride-1 — fully coalesced
!$acc parallel loop collapse(2)
do p = 1, n_mp
   do i = 1, nv
      sig(i, p) = sig(i, p) + dsig(i, p)
   end do
end do
```

**Rule:** All hot loops over particle arrays (P2G, G2P, stress update) should
use `collapse(2)` with the Voigt/spatial index as the inner loop.

---

## Derived Type Allocatables in ACC Regions

**This is the most common correctness bug in Fortran OpenACC code.**

Allocatable components of a derived type cannot be used directly inside
`!$acc parallel` regions. The device cannot follow the host-side allocation
pointer at runtime:

```fortran
! BAD — allocatable component of derived type inside acc region
!$acc parallel loop
do p = 1, n_mp
   set%sig(1, p) = 0.0_wp   ! undefined behaviour on device
end do
```

### Fix: associate before entering the region

Extract raw array references via `associate` on the host before the parallel
region. The compiler sees a plain array, not an allocatable component:

```fortran
associate(sig => set%sig, xp => set%xp, vp => set%vp)
   !$acc parallel loop collapse(2) default(present)
   do p = 1, n_mp
      do i = 1, nv
         sig(i, p) = 0.0_wp   ! OK — raw array reference
      end do
   end do
end associate
```

### Fix: pass raw arrays to GPU subroutines

GPU-callable subroutines should accept raw arrays, not derived types with
allocatable components. This also satisfies the `!$acc routine seq` requirement:

```fortran
! Good — pure function, GPU-callable
!$acc routine seq
pure function voigt_trace(sig) result(tr)
   real(wp), intent(in) :: sig(:)
   real(wp) :: tr
   tr = sig(V_11) + sig(V_22) + sig(V_33)
end function

! Bad — derived type argument in a routine seq context
!$acc routine seq
pure function voigt_trace(set, p) result(tr)
   type(mpm_particle_set_t), intent(in) :: set   ! allocatable inside → unsafe
   integer,                  intent(in) :: p
end function
```

**Rule:** GPU kernel subroutines and `!$acc routine seq` functions take raw
arrays (`real(wp), intent(in) :: sig(:)` or explicit shape) — never
`type(mpm_particle_set_t)`.

---

## Data Regions

### Keep data resident — copy once per run

Copying data to/from the GPU each time step dominates runtime for large
particle counts. Use `!$acc enter data` / `!$acc exit data` to keep arrays
resident across the entire time loop. Only copy back for output:

```fortran
! Extract raw arrays for acc data region
associate(xp => set%xp, vp => set%vp, sig => set%sig, &
          mass_I => grid%mass_I, momentum_I => grid%momentum_I)

   !$acc enter data copyin(xp, vp, sig, mass_I, momentum_I)

   do while (t < t_end)
      call map_p2g(...)
      call update_grid(...)
      call map_g2p(...)
      call update_stress(...)
   end do

   !$acc exit data copyout(xp, vp, sig)

end associate
```

### No allocatables inside data regions

All arrays must be allocated on the host **before** entering an
`!$acc enter data` region. Never allocate or deallocate inside a data region.

### default(present)

Use `default(present)` on parallel regions when arrays are already in a data
region. This prevents accidental silent copies:

```fortran
!$acc parallel loop default(present)
do p = 1, n_mp
   ...
end do
```

If an array is not present on the device and `default(present)` is set, the
runtime raises an error rather than silently copying — easier to debug.

---

## Routine Directives

Any pure function called inside an `!$acc parallel` region must be declared
with `!$acc routine seq`. Without it, the compiler will not generate device
code for the function and the call will fail or silently fall back to host:

```fortran
!$acc routine seq
pure function linear_shape_fn(xi) result(N)
   real(wp), intent(in) :: xi
   real(wp) :: N(2)
   N(1) = 0.5_wp * (1.0_wp - xi)
   N(2) = 0.5_wp * (1.0_wp + xi)
end function
```

This is identical to the pattern used in
`critical-soil-models/src/models/*/le_functions.f90`.

**Rule:** Any function that will be called from a hot loop gets
`!$acc routine seq` from day one, even before GPU implementation. It costs
nothing on CPU and documents intent.

---

## Polymorphism

Dynamic dispatch (`class()`, deferred procedures) cannot be used inside
`!$acc parallel` regions. The GPU cannot resolve vtable pointers at runtime.

The two-path design addresses this:

| Path | When | How |
|---|---|---|
| OOP (`class(mpm_shapefn_t)`) | CPU, research, debugging | Normal polymorphic dispatch |
| Pure functions (`!$acc routine seq`) | GPU hot loops | Static dispatch, no `class()` |

The algorithm-level dispatch (which variant: USL/USF/MUSL) happens once per
time step on the CPU — polymorphism there is fine and essentially free.

---

## AMD-Specific Notes

- Compiler: `gfortran -fopenacc -foffload=amdgcn-amdhsa` or `amdflang`
- OpenACC via GCC offload uses `libgomp` — link with `-fopenacc`
- Avoid `!$acc atomic` in inner loops — prefer `reduction` clauses
- NVIDIA-only extensions (`attributes(device)`, `<<<...>>>`, `threadIdx`) are
  forbidden — they will not compile with AMD toolchains
- Test correctness on CPU first (`fpm test`), then validate on GPU with
  `-fopenacc` — the CPU path is the ground truth

---

## Checklist for New Hot-Loop Code

Before writing a new kernel (P2G, G2P, stress update loop):

- [ ] Arrays passed as raw arguments — no `type(mpm_particle_set_t)` args
- [ ] `collapse(2)` on any loop over `(voigt/spatial component, particle)`
- [ ] `default(present)` on parallel regions
- [ ] Helper functions marked `!$acc routine seq`
- [ ] No `class()` dispatch inside the parallel region
- [ ] No `allocate` / `deallocate` inside the region
- [ ] No I/O inside the region (also violates `pure`)
