---
title: Time Integration
nav_order: 8
description: "Explicit SSP Runge–Kutta time stepping and the CFL stability condition for the cell-centred finite-volume residual."
---

# Time Integration
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Method of lines

The finite-volume spatial discretisation
([Finite-Volume Discretisation](finite-volume)) reduces the PDE system to
a system of ordinary differential equations in time:

$$
\frac{d \mathbf{U}_{ij}}{dt} \;=\; \mathbf{L}(\mathbf{U})_{ij},
\tag{MOL}
$$

where $$\mathbf{U}_{ij}$$ is the cell-averaged conservative state and
$$\mathbf{L}$$ is the **spatial residual**

$$
\mathbf{L}(\mathbf{U})_{ij}
\;=\;
- \frac{1}{V_{ij}}
  \Bigl[
    \bigl(\mathbf{F}^*_{i+1/2,j}\,r_{i+1/2} \;-\; \mathbf{F}^*_{i-1/2,j}\,r_{i-1/2}\bigr)\,\Delta z
    \;+\;
    \bigl(\mathbf{G}^*_{i,j+1/2} \;-\; \mathbf{G}^*_{i,j-1/2}\bigr)\,r_i\,\Delta r
  \Bigr]
\;+\; \mathbf{S}_{p,\,ij}.
$$

Because the FV cell volume is the only "mass" weight that appears, the
update in (MOL) is **already pre-multiplied** by $$V_{ij}^{-1}$$. There
is no global mass matrix to invert and no row-sum lumping is required —
the per-cell update is a scalar arithmetic operation. The same residual
structure applies to the passive scalar $$\rho Y$$ and is computed in
the same `compute_residual!` call.

Explicit time integration of (MOL) avoids any nonlinear solve per step
and is the natural choice for hyperbolic problems whose maximum stable
$$\Delta t$$ is set by the CFL bound rather than by stiffness.

---

## Explicit time-stepping schemes

Two explicit schemes are implemented in
[`src/timeint.jl`](https://github.com/jman87/cosmo/blob/main/src/timeint.jl);
both apply the same per-cell update kernel to the conservative state
$$\mathbf{U}$$ and the passive scalar $$\rho Y$$ in lock-step.

### Forward Euler (1st order)

$$
\mathbf{U}^{n+1}_{ij}
\;=\;
\mathbf{U}^n_{ij} \;+\; \Delta t\,\mathbf{L}(\mathbf{U}^n)_{ij}.
$$

First-order accurate in time. Provided for debugging and verification
only — its formal accuracy is mismatched with the 2nd-order MUSCL
spatial reconstruction, and it has the same CFL constraint as the
higher-order schemes, so it is strictly worse than SSP-RK2 in production.

The implementation is `step_forward_euler!`. The `forward_euler` /
`fe` / `euler` aliases in the JSON input all select this scheme.

---

### SSP-RK2 / Heun (2-stage, 2nd order) — default

The **strong-stability-preserving** Runge–Kutta scheme of Shu & Osher
(1988), also known as the Heun predictor-corrector, is the default
integrator:

$$
\mathbf{U}^{(1)}_{ij}
\;=\; \mathbf{U}^n_{ij} \;+\; \Delta t\,\mathbf{L}(\mathbf{U}^n)_{ij},
$$

$$
\mathbf{U}^{n+1}_{ij}
\;=\; \tfrac{1}{2}\,\mathbf{U}^n_{ij}
\;+\; \tfrac{1}{2}\,\bigl[\,\mathbf{U}^{(1)}_{ij}
\;+\; \Delta t\,\mathbf{L}(\mathbf{U}^{(1)})_{ij}\,\bigr].
$$

Equivalent to the standard 2nd-order RK method but expressed in
**Shu–Osher (convex-combination) form** so that the SSP property is
manifest: the final update is a convex combination of $$\mathbf{U}^n$$
and a forward-Euler step of $$\mathbf{U}^{(1)}$$.

**SSP property.** If the forward-Euler update is strongly stable
(does not increase a chosen norm) under a step bound
$$\Delta t \leq \Delta t_\text{FE}$$, then SSP-RK2 is strongly stable
under the same bound. Combined with a TVD spatial reconstruction
(MUSCL-minmod, see
[Riemann Solver and Reconstruction](riemann-reconstruction)) this
guarantees no spurious oscillation growth across shocks: the discrete
total variation of the solution does not increase from one time step to
the next.

**Cost.** Two residual evaluations per step. Each residual evaluation
performs one halo exchange on $$\mathbf{U}$$, one on $$\rho Y$$, one BC
pass, one full MUSCL+HLLC face-flux assembly, and one residual
assembly.

The implementation is `step_ssp_rk2!`. The aliases `ssp_rk2` /
`rk2` / `heun` all select this scheme.

---

### Why no SSP-RK3

The previous Python sibling solver also offered an SSP-RK3 option
(Shu–Osher 3-stage, 3rd order). The Julia implementation deliberately
omits SSP-RK3:

* The spatial scheme is at most 2nd-order accurate (MUSCL+minmod), so a
  3rd-order time stepper does not improve the order of the combined
  scheme — at shocks both schemes drop to 1st order and at smooth
  extrema the limiter caps spatial accuracy at 1st order.
* SSP-RK3 costs 3 residual evaluations per step versus 2 for SSP-RK2,
  a 50 % runtime increase for no measurable gain on practical blast
  problems.

Adding SSP-RK3 would be straightforward — the convex-combination
structure is identical — but the cost/benefit favours SSP-RK2.

---

## CFL stability condition

The maximum stable explicit step is bounded by the
**Courant–Friedrichs–Lewy (CFL)** condition. For a 2-D Euler-type
hyperbolic system on a structured grid the standard bound is

$$
\boxed{\;\;
\Delta t \;\leq\;
\sigma\,
\frac{h_\text{min}}{\displaystyle \max_{i,\,j}\,\bigl(\,\sqrt{u_r^2 + u_z^2} \;+\; c\,\bigr)_{ij}}
\;\;}
\tag{CFL}
$$

with

* $$\sigma$$ the CFL safety factor (`time.cfl` in JSON; default 0.4),
* $$h_\text{min} = \min(\Delta r,\, \Delta z)$$,
* the maximum taken over **interior cells** of the rank,
  then **MPI-reduced** with `MPI.Allreduce(..., MAX)` to a global maximum.

The face-normal-only form of the bound (the more relaxed
$$\sigma\,\min[\Delta r/\max(|u_r|+c),\,\Delta z/\max(|u_z|+c)]$$) would
give a slightly larger admissible step on highly anisotropic flow
fields, but the simpler $$h_\text{min}/\max(|\mathbf{u}|+c)$$ form is
robust in 2-D and is what the production code uses.

The local sound speed $$c$$ is the EOS-blended estimate
([Equations of State](eos))

$$
c \;=\; (1 - \lambda)\,c_\text{air} \;+\; \lambda\,c_\text{JWL},
$$

evaluated from $$(\rho, p)$$ at the cell centre. The implementation is
`compute_dt` in `src/timeint.jl`.

### Detonation-velocity safeguard

The gamma-equivalent JWL sound speed used by the CFL estimator can
under-estimate the true detonation-front speed by a wide margin: the
ratio $$D / c_\text{JWL}$$ for unburned TNT at $$e_\text{CJ}$$ exceeds
unity by roughly an order of magnitude. Using only the cell-centre
$$|\mathbf{u}| + c$$ during the burn would let the explicit step grow
larger than $$h_\text{min}/D$$, which would fail to advance the
burn front correctly and crash within a few steps.

The solver therefore **floors** the maximum global wave speed by the
user-supplied detonation velocity $$D$$ whenever any cell on any rank
still holds explosive material:

$$
v_\text{max}
\;=\;
\max\!\Bigl(\,
  \max_{i,j}\bigl(|\mathbf{u}_{ij}| + c_{ij}\bigr),\;\; D
\,\Bigr)
\quad
\text{if}\;\;
\max_{i,j}\,\lambda_{ij} > 0
\;\;\text{(global max via Allreduce).}
$$

Once the burn finishes everywhere ($$\max\,\lambda_{ij} = 0$$ on all
ranks — the burn front has swept all charge cells and lambda has been
reset by the IC, which only holds in Brode mode where lambda is
identically 1; in programmed-burn mode lambda remains 1 inside the
charge for the rest of the run, so this safeguard remains active for
the duration of any run that contains an explosive charge), the
$$|\mathbf{u}| + c$$ part is used unmodified.

### Hard upper bound: `dt_max`

An optional `time.dt_max` JSON parameter imposes a hard ceiling on
$$\Delta t$$ independent of the CFL bound. This is occasionally useful
when the early-time dynamics produce a much larger admissible step than
the user wants (e.g. for stable temporal sampling of an output field).
By default `dt_max` is `Inf` and is inactive.

### Last-step shrink

The actual step taken each iteration is

$$
\Delta t_\text{step} \;=\; \min\bigl(\Delta t_\text{CFL},\;\; \Delta t_\text{max},\;\; t_\text{end} - t\bigr).
$$

The third term ensures the run terminates **exactly** at $$t_\text{end}$$
without a final overshoot.

---

## Update of the burn fraction

Before each outer time step the reaction-progress field is refreshed
analytically:

```
update_burn!(state.lambda, grid, charge, t)
exchange_halo_scalar!(state.lambda, grid, topo)
```

This call evaluates the kinematic burn rule (see
[Programmed Burn Model](burn-model)) at every cell centre and overwrites
$$\lambda_{ij}$$ where the front has just arrived. The burn fraction
is then **frozen for the duration of the outer step**, so both stages
of the SSP-RK2 update see the same $$\lambda_{ij}$$. Holding $$\lambda$$
fixed within the RK is consistent because the programmed-burn front is
purely kinematic and does not depend on the flow field; updating it
mid-step would couple the two RK stages through a non-smooth function
and destroy the formal 2nd-order accuracy.

After the analytic update the lambda field is halo-exchanged so that
each rank's ghost layers contain the correct neighbour values for the
HLLC EOS evaluation on inter-rank faces.

---

## Boundary-condition enforcement within the RK

The conservative-state ghost layers are filled at the start of every
RK stage (not just every outer step) so that each residual evaluation
sees a consistent state on every face:

```
exchange_halos!(U, grid, topo)            # fill internal-seam ghosts from neighbours
exchange_halo_scalar!(rhoY, grid, topo)
apply_bcs!(U, grid, topo;                 # fill physical-edge ghosts (mirror or zero-grad)
           reflect_r_min=..., reflect_z_min=...)
apply_scalar_bcs!(rhoY, grid, topo)
compute_fluxes!(state, grid)              # MUSCL + HLLC at every interior face
compute_residual!(state, grid)            # divergence + geometric source
```

The enforced BCs (mirror for the symmetry axis and the rigid ground,
zero-gradient outflow elsewhere) are detailed in
[Boundary Conditions](boundary-conditions). Because the BCs operate on
ghost cells rather than on interior DOFs, no separate constraint
projection on the interior state is required — the BCs propagate into
the interior through the MUSCL stencil and the HLLC flux on the next
face evaluation.
