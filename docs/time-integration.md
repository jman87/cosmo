---
title: Time Integration
nav_order: 8
description: "SSP Runge-Kutta schemes, lumped mass matrix, and CFL stability condition."
---

# Time Integration
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Method of lines

The spatial discretization converts the PDEs into the semi-discrete system

$$
\mathbf{M}\,\dot{\mathbf{U}} = \mathbf{R}(\mathbf{U}, t),
\tag{MOL}
$$

where:

* $$\mathbf{U}(t)$$ is the global vector of all DOF values for $$(\rho, \rho u_r, \rho u_z, \rho E)$$,
* $$\mathbf{M}$$ is the consistent mass matrix (or its row-sum lumped approximation), and
* $$\mathbf{R}(\mathbf{U})$$ is the assembled spatial residual from the weak form.

Explicit time integration of equation (MOL) avoids solving a nonlinear system each step
and is natural for hyperbolic problems.

---

## Lumped mass matrix

### Why lumping?

The **consistent mass matrix** $$M_{ij} = \int_\Omega \phi_i \phi_j\, r\,dr\,dz$$ is
banded but not diagonal; inverting it at each time step is expensive.  **Row-sum lumping**
replaces $$M_{ij}$$ with the diagonal matrix

$$
M_i^L = \sum_j M_{ij} = \int_\Omega \phi_i \sum_j \phi_j\, r\,dr\,dz = \int_\Omega \phi_i\, r\,dr\,dz,
$$

where the last equality uses the partition of unity property
$$\sum_j \phi_j(\mathbf{x}) = 1$$.  The explicit update then reduces to

$$
U_i^{n+1} = U_i^n + \Delta t\, \frac{R_i}{M_i^L} \quad \forall i,
$$

which is a per-DOF scalar division — no global linear solve required.

### Symmetry-axis entries

At nodes on the symmetry axis $$r = 0$$, the cylindrical measure $$r\,dr\,dz = 0$$ makes
the row-sum exactly zero.  These entries are floored to a small fraction of the smallest
off-axis row-sum before the division, preventing blow-up.  The axisymmetric residual
$$R_i$$ is also very small at axis nodes (the same $$r$$ factor appears in the flux integral),
so the resulting update $$\Delta U_i \approx 0$$ is physically correct.

### Assembly

The lumped mass is assembled once at the start of the run by constructing the consistent
mass matrix $$\mathbf{M}$$ and multiplying by the all-ones vector:

$$
\mathbf{m}^L = \mathbf{M}\,\mathbf{1}.
$$

Separate lumped masses are built for the scalar (Q2) and vector (Q2v) spaces.

---

## Explicit time-stepping schemes

### Forward Euler (1st order)

$$
\mathbf{U}^{n+1} = \mathbf{U}^n + \Delta t\,(\mathbf{M}^L)^{-1}\mathbf{R}(\mathbf{U}^n).
$$

First-order accurate in time.  Provided for debugging and verification only; it is
inefficient in practice because the CFL constraint requires the same small $$\Delta t$$
as higher-order schemes.

---

### SSP-RK2 (2-stage, 2nd order) — default

The **Strong-Stability-Preserving Runge-Kutta** scheme of Shu & Osher (1988), also
called the **Heun predictor-corrector**, is the default time integrator:

$$
\mathbf{U}^{(1)} = \mathbf{U}^n + \Delta t\,(\mathbf{M}^L)^{-1}\mathbf{R}(\mathbf{U}^n),
$$

$$
\mathbf{U}^{n+1} = \tfrac{1}{2}\,\mathbf{U}^n
+ \tfrac{1}{2}\!\left[\mathbf{U}^{(1)} + \Delta t\,(\mathbf{M}^L)^{-1}\mathbf{R}(\mathbf{U}^{(1)})\right].
$$

This is equivalent to the standard 2nd-order Runge-Kutta method but expressed in
Shu-Osher (convex-combination) form to make the SSP property transparent:
$$\mathbf{U}^{n+1}$$ is a convex combination of $$\mathbf{U}^n$$ and a forward-Euler update
of $$\mathbf{U}^{(1)}$$.

**SSP property**: If the forward-Euler update is strongly stable (does not increase a
chosen norm) under a time-step restriction $$\Delta t \leq \Delta t_\text{FE}$$, then
SSP-RK2 is strongly stable under the same restriction.  For TVD (total-variation-
diminishing) norms this guarantees no spurious oscillation growth.

**Cost**: 2 residual evaluations per step.

---

### SSP-RK3 (3-stage, 3rd order)

The 3-stage, 3rd-order Shu-Osher scheme:

$$
\mathbf{U}^{(1)} = \mathbf{U}^n + \Delta t\,(\mathbf{M}^L)^{-1}\mathbf{R}(\mathbf{U}^n),
$$

$$
\mathbf{U}^{(2)} = \tfrac{3}{4}\,\mathbf{U}^n
+ \tfrac{1}{4}\!\left[\mathbf{U}^{(1)} + \Delta t\,(\mathbf{M}^L)^{-1}\mathbf{R}(\mathbf{U}^{(1)})\right],
$$

$$
\mathbf{U}^{n+1} = \tfrac{1}{3}\,\mathbf{U}^n
+ \tfrac{2}{3}\!\left[\mathbf{U}^{(2)} + \Delta t\,(\mathbf{M}^L)^{-1}\mathbf{R}(\mathbf{U}^{(2)})\right].
$$

**SSP coefficient** $$\mathcal{C} = 1$$ (same CFL restriction as forward Euler).
**Cost**: 3 residual evaluations per step.

SSP-RK3 provides one order of accuracy higher than SSP-RK2 at a 50% increase in
computational cost per step.  The improvement is most visible in smooth regions of the
flow; across shocks the order reduction (due to the discontinuity) makes the practical
difference small.

---

## CFL stability condition

Explicit time integrators are conditionally stable.  The maximum stable time step is
bounded by the **Courant–Friedrichs–Lewy (CFL) condition**:

$$
\Delta t = \text{CFL} \cdot \frac{h_\text{min}}{\max_{i}(|\mathbf{u}_i| + c_i)},
$$

where:

| Symbol | Meaning |
|--------|---------|
| $$\text{CFL}$$ | Safety factor $$\in (0, 1]$$; default 0.4 |
| $$h_\text{min}$$ | Global minimum cell diameter (MPI-reduced) |
| $$\|\mathbf{u}_i\|$$ | Velocity magnitude at DOF $$i$$ |
| $$c_i$$ | Local sound speed at DOF $$i$$ |

The maximum wave speed $$\max(|\mathbf{u}| + c)$$ is computed over all scalar DOFs and
MPI-reduced with `MPI.MAX`.

### Detonation velocity bound

The air-$$\gamma$$ sound speed is a poor estimate inside JWL-product cells.  Whenever any
material indicator $$\lambda > 0$$ is present on any rank, the solver additionally bounds
the wave speed by the user-supplied detonation velocity $$D$$:

$$
v_\text{max} = \max\!\left(\max_i(|\mathbf{u}_i| + c_i),\; D\right).
$$

This prevents the CFL step from becoming larger than $$h_\text{min}/D$$ while explosive
products are still expanding from the charge.

### Hard upper bound

An optional `dt_max` parameter in the input JSON imposes a hard ceiling on $$\Delta t$$,
independent of the CFL condition.  This is useful when early-time dynamics (e.g., the
initial detonation pulse) require a coarser time step than the CFL condition would impose
for stability but a finer step than physics requires.

---

## Burn indicator update

Before each outer time step the material indicator is updated from the analytic
programmed-burn formula (see [Programmed Burn Model](burn-model)):

```
_update_burn(state, burn, t)
```

This call re-evaluates $$\lambda(r, z, t)$$ at all DG0 DOF coordinates.  The `mat` field
is **not** updated within the RK sub-stages; it is kept fixed at the value corresponding
to the start of the time step.  This avoids inconsistencies between the burn front
position used in different stage evaluations.

---

## Symmetry boundary condition enforcement

After each stage update, the radial momentum $$\rho u_r$$ is zeroed at all DOFs lying on
the $$r = 0$$ axis, and (when $$z_\text{min} = 0$$) the axial momentum $$\rho u_z$$ is zeroed
at all DOFs on the $$z = 0$$ plane:

$$
(\rho u_r)\big|_{r=0} = 0, \qquad (\rho u_z)\big|_{z=0} = 0.
$$

These constraints enforce the symmetry conditions exactly, preventing numerical drift
of the velocity components that should be zero by symmetry.
