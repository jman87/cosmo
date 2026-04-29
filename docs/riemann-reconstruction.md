---
title: Riemann Solver and Reconstruction
nav_order: 7
description: "HLLC approximate Riemann solver, MUSCL spatial reconstruction, and the minmod slope limiter."
---

# Riemann Solver and Reconstruction
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Why a Riemann solver?

The compressible Euler equations admit discontinuous (shock and contact)
solutions even from smooth initial data. A naive central-difference
discretisation of these equations produces spurious post-shock
oscillations (a Gibbs-like phenomenon at strong gradients) and is
unconditionally unstable in the absence of artificial dissipation.

The Godunov framework — of which finite-volume + Riemann is the
canonical realisation — replaces the local Taylor expansion with the
exact (or approximate) solution of a 1-D **Riemann problem** at each
cell face:

$$
\mathbf{U}(r, t = 0)
\;=\;
\begin{cases}
  \mathbf{U}^L, & r < 0, \\
  \mathbf{U}^R, & r > 0.
\end{cases}
$$

Sampling that self-similar solution along $$r/t = 0$$ gives the
numerical flux $$\mathbf{F}^*$$ used in the cell-update equation
([Finite-Volume Discretisation](finite-volume)). Upwinding is therefore
imposed *exactly* by the wave structure of the conservation law itself,
which removes the need for any artificial-viscosity stabilisation. This
is the central reason the Julia solver does not use a von
Neumann–Richtmyer artificial bulk viscosity (in contrast to a
continuous-Galerkin FEM approach).

Solving the exact Riemann problem at every face every step is too
expensive; **approximate Riemann solvers** retain the upwind structure
while replacing the full nonlinear wave system with a small number of
simple waves. COSMO uses the **HLLC** solver of Toro, Spruce & Speares
(1994) — a three-wave (left, contact, right) approximation that
preserves the contact discontinuity exactly and is robust at very strong
shocks, both of which are essential for blast applications.

---

## HLLC approximate Riemann solver

### Wave structure

HLLC assumes the exact solution to the Riemann problem on the face
consists of **three waves** separating four constant states:

```
   t                  S_L   S*    S_R
   │                   ╲    │    ╱
   │   U_L              ╲   │   ╱     U_R
   │                     ╲  │  ╱
   │              U*_L    ╲ │ ╱   U*_R
   │                       ╲│╱
   ─┴─────────────────────────────────► r
                          r = 0
```

* $$S_L,\,S_R$$ are the left- and right-going acoustic / shock waves;
* $$S^*$$ is the contact wave (slip line);
* $$\mathbf{U}^*_L,\,\mathbf{U}^*_R$$ are the two intermediate "star"
  states adjacent to the contact.

The contact has the property $$u^*_n = S^*$$ and $$p^*_L = p^*_R$$. HLLC
constructs the two star states algebraically from the Rankine–Hugoniot
conditions across each acoustic wave, given the wave speeds.

### Wave-speed estimate (Davis)

The solver uses the **Davis** wave-speed estimate (Toro 2009 Eq. 10.48),
which is robust at very strong contrasts:

$$
S_L \;=\; \min\!\bigl(u_n^L - c^L,\;\; u_n^R - c^R\bigr),
$$

$$
S_R \;=\; \max\!\bigl(u_n^L + c^L,\;\; u_n^R + c^R\bigr),
$$

where $$u_n$$ is the face-normal velocity and $$c$$ is the local sound
speed evaluated from the **EOS-blended** sound-speed formula
([Equations of State](eos)) using the local burn fraction $$\lambda$$:

$$
c \;=\; (1 - \lambda)\,c_\text{air} \;+\; \lambda\,c_\text{JWL}.
$$

The contact-wave speed (Toro 2009 Eq. 10.37) is

$$
S^* \;=\;
\frac{p^R - p^L
      \;+\; \rho^L\,u_n^L\,(S_L - u_n^L)
      \;-\; \rho^R\,u_n^R\,(S_R - u_n^R)}
     {\rho^L\,(S_L - u_n^L) \;-\; \rho^R\,(S_R - u_n^R)}.
$$

### Star-state assembly and flux

For axis-aligned faces with normal $$\hat{\mathbf{n}}$$ pointing in
direction `dim` (1 = $$r$$, 2 = $$z$$), denote by $$u_n$$ the normal
component and $$u_t$$ the tangential component of velocity. The
intermediate density and energy on the $$L$$-side of the contact are

$$
\rho^*_L \;=\; \rho^L\,\frac{S_L - u_n^L}{S_L - S^*},
$$

$$
E^*_L \;=\; E^L \;+\; (S^* - u_n^L)\!\left[\,S^* \;+\; \frac{p^L}{\rho^L\,(S_L - u_n^L)}\,\right],
$$

with the analogous expressions on the $$R$$-side. The intermediate
**conservative state** has the normal-momentum component aligned with
$$S^*$$ and the tangential component preserved from the upwind side; the
density and total energy come from the relations above. The HLLC
**numerical flux** is then the closed-form selection rule

$$
\mathbf{F}^*_\text{HLLC} \;=\;
\begin{cases}
  \mathbf{F}^L,                                                & 0 \;\leq\; S_L, \\
  \mathbf{F}^L \;+\; S_L\,(\mathbf{U}^*_L - \mathbf{U}^L),     & S_L \;<\; 0 \;\leq\; S^*, \\
  \mathbf{F}^R \;+\; S_R\,(\mathbf{U}^*_R - \mathbf{U}^R),     & S^* \;<\; 0 \;\leq\; S_R, \\
  \mathbf{F}^R,                                                & S_R \;<\; 0,
\end{cases}
$$

where $$\mathbf{F}^{L,R}$$ are the *physical* fluxes evaluated at
$$\mathbf{U}^{L,R}$$ in the face-normal direction. The implementation
lives in `hllc_flux` in
[`src/riemann.jl`](https://github.com/jman87/cosmo/blob/main/src/riemann.jl).

### Mixed-EOS handling at the burn front

Within a few cells of the detonation front the air-side and
products-side EOS differ. HLLC is invoked with separate burn fractions
$$\lambda^L,\lambda^R$$ on each side of the face; the primitives and
sound speeds are evaluated with the local blend on each side. Treating
the burn front as a contact-like discontinuity in this way is the
simplest physically-consistent way to apply HLLC across an EOS jump and
preserves the contact-wave structure that the solver was chosen to
resolve.

The burn fraction itself is **piecewise-constant per cell**
(cell-centred); the face value is taken from the **owning cell** on each
side of the face. A more sophisticated approach would interpolate
$$\lambda$$ to the face, but the burn front is binary
([Programmed Burn Model](burn-model)) so any sub-cell interpolation
would introduce its own approximation. Within a single time step at most
a few cells experience a $$\lambda$$ transition, so the
piecewise-constant choice has negligible global effect.

---

## MUSCL spatial reconstruction

A 1st-order Godunov scheme would use the cell averages directly as
left/right states at the face: $$\mathbf{U}^L = \mathbf{U}_{i,j}$$,
$$\mathbf{U}^R = \mathbf{U}_{i+1,j}$$. This is unconditionally robust
but produces an excessively diffusive shock (a 1st-order scheme has
a numerical viscosity proportional to $$\Delta x$$).

**MUSCL** (Monotone Upstream-centred Schemes for Conservation Laws,
van Leer 1979) raises the spatial order to 2 by reconstructing each
cell as a **piecewise-linear** profile and extrapolating to the face
along the linear interpolant. To preserve monotonicity at shocks and
prevent overshoots, the linear slope is **limited** by a TVD slope
limiter that drops the slope to zero near extrema.

### Centred slope and minmod limiter

For each component of the conservative state and each cell $$(i, j)$$,
form the *backward* and *forward* differences along the reconstruction
direction:

$$
\Delta^- \;=\; \mathbf{U}_{i,j} \;-\; \mathbf{U}_{i-1,j},
\qquad
\Delta^+ \;=\; \mathbf{U}_{i+1,j} \;-\; \mathbf{U}_{i,j}.
$$

The **minmod-limited centred slope** is

$$
\sigma_{i,j}
\;=\; \mathrm{minmod}(\Delta^-,\,\Delta^+)
\;\equiv\;
\begin{cases}
  \Delta^-, & \Delta^-\,\Delta^+ > 0 \;\;\text{and}\;\; |\Delta^-| \;\leq\; |\Delta^+|, \\
  \Delta^+, & \Delta^-\,\Delta^+ > 0 \;\;\text{and}\;\; |\Delta^+| \;<\; |\Delta^-|, \\
  0,         & \Delta^-\,\Delta^+ \;\leq\; 0.
\end{cases}
$$

The limiter has three purposes: (i) when both differences have the same
sign and finite magnitude, take the smaller (in absolute value) — this
prevents overshoots; (ii) when they have opposite signs the cell sits at
a local extremum and the slope is killed (TVD); (iii) when one is zero
the slope is zero (preserves constants).

Among the standard TVD limiters (van Leer, monotonised central, superbee,
minmod), **minmod is the most diffusive** but also the most robust:
it never produces overshoots, never amplifies negative-density excursions
near very strong shocks, and is the limiter of choice when correctness
matters more than crispness. It is the natural default for explosive-blast
applications where the shock strengths span tens of orders of magnitude
in pressure.

### Face extrapolation

With slope $$\sigma$$ defined per cell, the face states on the radial
face between cells $$(i, j)$$ and $$(i+1, j)$$ are

$$
\mathbf{U}^L_{i+1/2,\,j} \;=\; \mathbf{U}_{i,j}   \;+\; \tfrac{1}{2}\,\sigma_{i,j},
$$

$$
\mathbf{U}^R_{i+1/2,\,j} \;=\; \mathbf{U}_{i+1,j} \;-\; \tfrac{1}{2}\,\sigma_{i+1,j}.
$$

The full **4-point stencil** for a single face is therefore

$$
(\mathbf{U}_{i-1,j},\;\mathbf{U}_{i,j},\;\mathbf{U}_{i+1,j},\;\mathbf{U}_{i+2,j}),
$$

which is why the solver requires $$n_g = 2$$ ghost layers
([Mesh and Grid](mesh-discretization)). The implementation lives in
`reconstruct_face` in
[`src/reconstruct.jl`](https://github.com/jman87/cosmo/blob/main/src/reconstruct.jl).

The same extrapolation is used on $$z$$-faces with the obvious change
of indices, and the four conservative-state components are
**reconstructed independently**. Reconstructing on conservative variables
is the simplest robust choice for a fixed-grid Eulerian solver; it is
mass- and momentum-monotone by construction, but can in principle let the
local pressure dip non-physically low through cancellation of the
momentum and total-energy reconstructions. Reconstructing on primitives
($$\rho,\,u_r,\,u_z,\,p$$) or on characteristic variables avoids that
mode but requires an extra cell-by-cell EOS evaluation. For TNT-class
blast problems the conservative-variable choice has been thoroughly
exercised and has not produced spurious states.

---

## TVD property and shock width

A scheme is **total-variation-diminishing (TVD)** if

$$
\sum_i \bigl|\,\mathbf{U}^{n+1}_{i+1} - \mathbf{U}^{n+1}_{i}\,\bigr|
\;\leq\;
\sum_i \bigl|\,\mathbf{U}^{n}_{i+1}  - \mathbf{U}^{n}_{i}\,\bigr|
$$

for every initial datum. The Harten/Sweby theorem states that minmod
combined with an SSP time integrator and a CFL-bounded step is TVD on a
1-D scalar conservation law, and the result extends in practice to the
multi-dimensional Euler system through the dimension-by-dimension
update. TVD implies that **no new local extrema are created** by the
update — i.e. no spurious post-shock oscillations.

The MUSCL+HLLC+minmod combination resolves an isolated normal shock over
roughly **3–4 cells**, with no post-shock ringing for any reasonable
choice of CFL and no need for additional artificial dissipation. The
shock width is independent of the shock strength (it is set by the
limiter, not by the wave amplitude).

---

## Why no artificial viscosity

Some compressible-flow codes — particularly continuous-Galerkin FEM and
classical Lagrangian hydrocodes — rely on a von Neumann–Richtmyer (VNR)
artificial bulk viscosity to smear shocks and damp oscillations. COSMO
does **not** use an artificial-viscosity term:

* The HLLC numerical flux supplies the upwind physics that VNR
  approximates by a phenomenological compressive viscosity.
* MUSCL+minmod supplies the high-order interior accuracy without
  introducing the dispersive Gibbs-like ringing that VNR is added to
  damp.
* The combined scheme is TVD, which VNR cannot guarantee.

Removing artificial viscosity also removes its three main costs:
(i) the unphysical heating from $$q\,\nabla\!\cdot\!\mathbf{u}$$ in
compressive zones, (ii) two additional tuning constants
($$c_q,\,c_l$$), and (iii) an additional CFL-like time-step bound from
the parabolic-like viscous term. The Riemann/MUSCL framework is
self-contained.
