# MPM Particle Update Schemes

Design notes and comparison of particle velocity update schemes. These govern
the G2P (grid-to-particle) transfer and in some cases also the P2G transfer.
The choice affects dissipation, stability, angular momentum conservation, and
particle data requirements.

---

## The Problem

After solving the nodal momentum equation we have `v_i^{n+1}`. We need to
transfer that back to particles. How we do that is the update scheme.

---

## Schemes

### PIC (Particle-In-Cell)

Set particle velocity directly from the interpolated nodal velocity. The
particle discards its own velocity history entirely.

```
v_p^{n+1} = Σ_i N_ip * v_i^{n+1}
```

**Pros:** Maximally stable. No particle noise.
**Cons:** Maximally dissipative. High-frequency modes (waves) are killed in a
handful of steps. Unusable for elastic dynamics.
**Extra particle storage:** None.

---

### FLIP (Fluid-Implicit Particle)

Update particle velocity by the nodal velocity *increment* rather than the
absolute velocity. Particles carry their own velocity history forward.

```
v_p^{n+1} = v_p^n + Σ_i N_ip * (v_i^{n+1} - v_i^n)
```

Requires `v_i^n` (old nodal velocity) to be saved before the grid update.
In `mpm_grid_t` this is the `vel_I` array — store it before `update_grid`,
use it during G2P.

**Pros:** Near-zero numerical dissipation. Good for elastic wave problems.
**Cons:** Particles can accumulate uncorrelated "particle-mode" velocities over
time, leading to noisy instability in long runs.
**Extra particle storage:** None. Extra grid storage: `v_i^n` (already planned
as `vel_I` in `mpm_grid_t`).

---

### α-FLIP (FLIP/PIC Blend) — **default for this solver**

Blend PIC and FLIP with a scalar parameter α ∈ [0, 1]:

```
v_p^{n+1} = α * [v_p^n + Σ_i N_ip * Δv_i]    ← FLIP part
           + (1-α) * [Σ_i N_ip * v_i^{n+1}]   ← PIC  part
```

where `Δv_i = v_i^{n+1} - v_i^n`.

- α = 1.0 → pure FLIP
- α = 0.0 → pure PIC
- α = 0.99 → standard geotechnical default (low dissipation, slow noise bleed)
- α = 1.0 → use for elastic wave validation (no artificial damping)

α lives in `mpm_update_config_t` and defaults to 0.99. Set to 1.0 for the
1D elastic wave test.

**Pros:** Tunable tradeoff between stability and dissipation.
**Cons:** α is problem-dependent; wrong choice degrades results.
**Extra particle storage:** None.

---

### APIC (Affine Particle-In-Cell)

Jiang et al. (2015). Each particle carries a velocity *and* an affine velocity
matrix `C_p` (the local velocity gradient). Both G2P **and** P2G differ from
the standard formulation — this is not a drop-in swap for G2P alone.

**G2P:**
```
v_p^{n+1} = Σ_i N_ip * v_i^{n+1}
C_p^{n+1} = Σ_i N_ip * v_i^{n+1} ⊗ (x_i - x_p)
```

**P2G (modified — affine term contributes to nodal momentum):**
```
p_i = Σ_p N_ip * m_p * (v_p + C_p * (x_i - x_p))
```

**Pros:**
- Zero numerical dissipation (like FLIP)
- No particle noise (like PIC)
- Exactly conserves angular momentum (FLIP does not)
- `C_p` also provides the velocity gradient needed for the strain rate, which
  can simplify the stress update

**Cons:**
- Extra particle storage: `C_p` is `ndim × ndim` per particle
- P2G loop is more expensive (extra dot product per particle-node pair)
- Requires care with GPU: `C_p` must be stored in SOA layout

**Extra particle storage:**

| Problem type | `C_p` components | Array shape |
|---|---|---|
| 1D | 1 (scalar) | `(1, n_mp)` |
| 2D | 4 | `(4, n_mp)` or `(2, 2, n_mp)` |
| 3D | 9 | `(9, n_mp)` or `(3, 3, n_mp)` |

Store as `Cp(ndim*ndim, n_mp)` (flattened, SOA) for GPU coalescing.

---

### XPIC (eXtended PIC)

Hammerquist & Nairn (2017). Higher-order correction terms added to PIC that
progressively reduce numerical dissipation.

```
XPIC(m): v_p^{n+1} = Σ_i N_ip * [v_i^{n+1} - Σ_{k=1}^{m} (-1)^k * correction_k]
```

XPIC(1) = standard PIC. Higher orders approach FLIP-like behaviour.

**Pros:** No extra particle storage. More accurate than PIC.
**Cons:** More complex implementation. Less common than α-FLIP or APIC.
**Status:** Deferred. Not planned for initial implementation.

---

## Comparison Table

| Scheme | Dissipation | Noise risk | Angular momentum | Extra particle data | P2G change? |
|---|---|---|---|---|---|
| PIC | High | None | No | None | No |
| FLIP | ~Zero | High | No | None | No |
| α-FLIP | Tunable | Low (α<1) | No | None | No |
| APIC | ~Zero | None | Yes | `C_p (ndim²)` | **Yes** |
| XPIC | Low–Medium | Low | No | None | No |

---

## Interaction with Algorithm Variants

Update scheme and algorithm variant (USL/USF/MUSL) are orthogonal — any
combination is valid. Common pairings:

| Algorithm | Scheme | Use case |
|---|---|---|
| USL | α-FLIP (α=0.99) | Standard geotechnical |
| USL | α-FLIP (α=1.0) | Elastic wave validation |
| USL | APIC | High-accuracy, angular momentum critical |
| MUSL | α-FLIP | Reduced velocity oscillation in geotechnics |

---

## Implementation Plan

### Stage 1 — α-FLIP (current)

Add `mpm_update_config_t` to hold α:

```fortran
type :: mpm_update_config_t
   real(wp) :: alpha = 0.99_wp  !! FLIP/PIC blend; 1.0=pure FLIP, 0.0=pure PIC
end type
```

G2P for α-FLIP:

```fortran
! save v_i^n before grid update (store in grid%vel_I before update_grid call)
! then in G2P:
do p = 1, n_mp
   v_flip = vp(:,p) + matmul(N_ip, dv_I)    ! FLIP increment
   v_pic  = matmul(N_ip, vel_I_new)          ! PIC interpolation
   vp(:,p) = alpha * v_flip + (1-alpha) * v_pic
end do
```

### Stage 2 — APIC (future)

Add `Cp(ndim*ndim, n_mp)` to `mpm_particle_set_t` — allocated only when
`update_config%scheme == SCHEME_APIC`. Requires modifying both G2P and P2G.

### Recommended parameter for each validation problem

| Problem | α | Reason |
|---|---|---|
| 1D elastic wave | 1.0 | No artificial damping of wave |
| 2D footing (quasi-static) | 0.99 | Suppress particle noise in long run |
| Slope failure (large deform.) | 0.99 | Stability over large deformations |
| Consolidation | 0.99 | Low dissipation needed for pore pressure |

---

## References

- Brackbill & Ruppel (1986) — original FLIP for fluids
- Sulsky et al. (1994, 1995) — MPM, original PIC transfer
- Jiang et al. (2015) — APIC: *The Affine Particle-in-Cell Method*, SIGGRAPH
- Hammerquist & Nairn (2017) — XPIC: *A new method for material point method
  particle updates that reduces noise and enhances stability*
- de Vaucorbeil et al. (2020) — MPM review covering FLIP/PIC/APIC variants
