---
title: Governing Equations
nav_order: 2
description: "Axisymmetric compressible Euler equations in strong conservative form."
---

# Governing Equations
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Problem geometry

COSMO models a detonation in a fluid domain that is rotationally symmetric
about the $$z$$-axis. It is therefore sufficient to solve on the 2-D
meridional half-plane $$\Omega \subset \{(r, z) \mid r \geq 0\}$$, where

* $$r$$ is the radial distance from the symmetry axis (horizontal in the mesh),
* $$z$$ is the axial coordinate (vertical in the mesh).

The full 3-D solution is recovered by rotating the 2-D solution through
$$2\pi$$ radians. Because all fields are independent of the azimuthal angle
$$\theta$$, no $$\partial / \partial \theta$$ terms appear. A 3-D fluid
element occupies the wedge $$dV = r\,dr\,d\theta\,dz$$, so a finite-volume
cell of footprint $$\Delta r \times \Delta z$$ centred at radius $$r_c$$ has
azimuthally-integrated volume

$$
V_{ij} = 2\pi\, r_c\, \Delta r\, \Delta z,
$$

and the $$2\pi$$ cancels uniformly across the conservation laws and the
divergence theorem. The solver therefore works with the **annular volume**

$$
V_{ij} \;=\; r_c\,\Delta r\,\Delta z
$$

and **annular face areas**

$$
A^{(r)}_{i\pm 1/2} = r_{i\pm 1/2}\,\Delta z,
\qquad
A^{(z)}_{j\pm 1/2} = r_c\,\Delta r,
$$

where $$r_{i\pm 1/2}$$ is the radius of the radial face between cells $$i$$ and
$$i\pm 1$$ and $$r_c$$ is the radius of the cell centre.

---

## Conserved variables

The solver advances four conserved scalars per cell as the components of the
state vector

$$
\mathbf{U} =
\begin{pmatrix}
  \rho \\
  \rho u_r \\
  \rho u_z \\
  \rho E
\end{pmatrix},
$$

stored in code as `state.U[k, i, j]` with $$k \in \{1,2,3,4\}$$, $$i$$ the
radial cell index, and $$j$$ the axial cell index. The interior cell range
spans $$i \in [n_g+1,\, n_g+n_r]$$ with $$n_g = 2$$ ghost cells on each side.

| Symbol | Meaning | Units |
|--------|---------|-------|
| $$\rho$$ | Mass density | $$\text{lbf}\cdot\text{s}^2/\text{in}^4$$ |
| $$\rho u_r$$ | Radial momentum density | $$\text{lbf}\cdot\text{s}/\text{in}^3$$ |
| $$\rho u_z$$ | Axial momentum density | $$\text{lbf}\cdot\text{s}/\text{in}^3$$ |
| $$\rho E$$ | Total energy density | $$\text{lbf}/\text{in}^2$$ |

where $$E = e + \tfrac{1}{2}(u_r^2 + u_z^2)$$ is the specific total energy,
$$e$$ is the specific internal energy ($$\text{in}^2/\text{s}^2$$), and
$$u_r$$, $$u_z$$ are the velocity components ($$\text{in}/\text{s}$$).

In addition to $$\mathbf{U}$$, the solver carries

* a **reaction-progress variable** $$\lambda(r, z, t) \in [0, 1]$$ used to
  blend the JWL and ideal-gas pressures, set by the programmed-burn rule;
* a **passive scalar** $$\rho Y$$, where $$Y \in [0, 1]$$ is the mass fraction
  of explosive-origin material. $$\rho Y$$ is conserved and advected
  separately so the post-burn material tag travels with the products as they
  expand into the air.

---

## Strong conservative form

The inviscid compressible Euler equations in 2-D axisymmetric coordinates are

$$
\frac{\partial \rho}{\partial t}
  + \frac{1}{r}\frac{\partial (r\,\rho u_r)}{\partial r}
  + \frac{\partial (\rho u_z)}{\partial z}
  = 0
\tag{mass}
$$

$$
\frac{\partial (\rho u_r)}{\partial t}
  + \frac{1}{r}\frac{\partial\!\left(r\,(\rho u_r^2 + p)\right)}{\partial r}
  + \frac{\partial (\rho u_r u_z)}{\partial z}
  - \frac{p}{r}
  = 0
\tag{$$r$$-momentum}
$$

$$
\frac{\partial (\rho u_z)}{\partial t}
  + \frac{1}{r}\frac{\partial (r\,\rho u_r u_z)}{\partial r}
  + \frac{\partial (\rho u_z^2 + p)}{\partial z}
  = 0
\tag{$$z$$-momentum}
$$

$$
\frac{\partial (\rho E)}{\partial t}
  + \frac{1}{r}\frac{\partial\!\left(r\,(\rho E + p)\,u_r\right)}{\partial r}
  + \frac{\partial\!\left((\rho E + p)\,u_z\right)}{\partial z}
  = 0
\tag{energy}
$$

The passive-scalar advection equation for the explosive mass fraction is the
analogue of the mass equation:

$$
\frac{\partial (\rho Y)}{\partial t}
  + \frac{1}{r}\frac{\partial(r\,\rho u_r Y)}{\partial r}
  + \frac{\partial(\rho u_z Y)}{\partial z}
  = 0.
\tag{species}
$$

The reaction-progress variable $$\lambda$$ is **not** evolved by a
PDE — it is set point-wise by the programmed-burn rule
(see [Programmed Burn Model](burn-model)) and provides closure for the EOS.

---

## Compact flux–source form

Define the radial and axial physical flux vectors

$$
\mathbf{F}(\mathbf{U}) =
\begin{pmatrix}
  \rho u_r \\
  \rho u_r^2 + p \\
  \rho u_r u_z \\
  (\rho E + p)\,u_r
\end{pmatrix},
\qquad
\mathbf{G}(\mathbf{U}) =
\begin{pmatrix}
  \rho u_z \\
  \rho u_r u_z \\
  \rho u_z^2 + p \\
  (\rho E + p)\,u_z
\end{pmatrix}
$$

and the geometric pressure source

$$
\mathbf{S}_p(\mathbf{U}) =
\begin{pmatrix} 0 \\ p/r \\ 0 \\ 0 \end{pmatrix}.
$$

The system is then

$$
\frac{\partial \mathbf{U}}{\partial t}
  + \frac{1}{r}\frac{\partial (r\,\mathbf{F})}{\partial r}
  + \frac{\partial \mathbf{G}}{\partial z}
  = \mathbf{S}_p.
$$

The geometric source $$p/r$$ in the radial-momentum equation is the
**hoop-stress** (or centripetal pressure) term. It arises because the
divergence of the radial flux in cylindrical coordinates contains a $$1/r$$
metric factor:

$$
\frac{1}{r}\frac{\partial}{\partial r}\!\left[r\,(\rho u_r^2 + p)\right]
\;=\;
\frac{1}{r}\frac{\partial(r\,\rho u_r^2)}{\partial r}
\;+\;
\frac{\partial p}{\partial r}
\;+\;
\frac{p}{r}.
$$

The first two terms are absorbed into a divergence on the left-hand side,
leaving $$-p/r$$ as an apparent source on the LHS — equivalently $$+p/r$$
on the RHS. Physically it represents the net outward pressure force per unit
volume on a wedge-shaped fluid element of finite angular extent. The
discretisation of this term on each cell is documented in
[Finite-Volume Discretisation](finite-volume).

---

## Primitive variable recovery

The Riemann solver, the EOS, and the CFL estimator all require the
primitives $$(\rho, u_r, u_z, p)$$ recovered from $$\mathbf{U}$$:

$$
u_r = \frac{(\rho u_r)}{\rho}, \qquad
u_z = \frac{(\rho u_z)}{\rho},
$$

$$
e = \frac{(\rho E)}{\rho} - \tfrac{1}{2}(u_r^2 + u_z^2),
$$

$$
p = p(\rho,\,e;\,\lambda) \quad \text{from the EOS blend (see Equations of State).}
$$

In code these calculations live in `primitives(U, lambda)` in
[`src/eos.jl`](https://github.com/jman87/cosmo/blob/main/src/eos.jl). To
prevent NaN propagation in extremely cold or near-vacuum cells the solver
applies floors before evaluation:

| Floor | Value | Purpose |
|-------|-------|---------|
| $$\rho_\text{floor}$$ | $$10^{-12}\ \text{lbf}\cdot\text{s}^2/\text{in}^4$$ | Prevent division by zero in $$u = (\rho u)/\rho$$ |
| $$e_\text{floor}$$ | $$10^{-6}\ \text{in}^2/\text{s}^2$$ | Keep $$e > 0$$ for the EOS evaluation |
| $$p_\text{floor}$$ | $$10^{-6}\ \text{psi}$$ | Returned from `mix_pressure` if the blend yields a negative value |

The floors only activate in pathological, clearly non-physical cells; for any
properly-initialised case they are inactive everywhere on the interior mesh.

---

## Reaction-progress-weighted pressure

A single-fluid Eulerian description requires a single pressure per cell. The
solver uses a linear blend of the air and JWL pressures weighted by the
local reaction-progress variable:

$$
p \;=\; (1 - \lambda)\, p_\text{air}(\rho, e) \;+\; \lambda\, p_\text{JWL}(\rho, e),
\qquad \lambda \in [0, 1].
$$

A blend of two analytic EOS is the simplest physically-consistent way to
treat the brief region where the burn front is cutting through a cell —
it avoids the expense and the numerical fragility of solving a stiff
reactive flow system while still giving the products their correct stiff
(JWL) response once the burn has finished. The full discussion is in
[Equations of State](eos).

---

## Sound speed used by the CFL bound

The maximum stable explicit time step (see [Time Integration](time-integration))
is set by the largest cell-centre wave speed $$|u| + c$$. The squared
isentropic sound speed for a general EOS is

$$
c^2 = \left.\frac{\partial p}{\partial \rho}\right|_s
    = \left.\frac{\partial p}{\partial \rho}\right|_e
    + \frac{p}{\rho^2}\left.\frac{\partial p}{\partial e}\right|_\rho,
$$

which simplifies for the **ideal gas** ($$p = (\gamma-1)\rho e$$) to

$$
c_\text{air}^2 = \gamma\,(\gamma - 1)\,e \;=\; \frac{\gamma\,p}{\rho}.
$$

For **JWL** the exact $$c^2$$ is the sum of the exponentials' partial
derivatives plus the Grüneisen contribution and is expensive to evaluate.
Because the CFL bound only needs an *over-estimate* of the wave speed, the
solver uses the gamma-equivalent

$$
c_\text{JWL}^2 \;\approx\; (\omega + 1)\,\frac{p}{\rho},
$$

which is the standard production-blast-code approximation and is bounded
above the true value across the JWL parameter ranges encountered for TNT.
The blended sound speed is the same convex combination as the pressure:

$$
c \;=\; (1 - \lambda)\,c_\text{air} \;+\; \lambda\,c_\text{JWL}.
$$

A separate **detonation-velocity safeguard** floors the maximum wave speed by
the user-supplied $$D$$ whenever any cell on any rank still holds explosive
material ($$\lambda > 0$$ on a non-empty subset). This protects against the
gamma-equivalent $$c_\text{JWL}$$ under-estimating the true detonation-front
speed during the burn. See [Time Integration](time-integration).

---

## Assumptions and limitations

| Assumption | Implication |
|------------|-------------|
| Inviscid flow | No viscous stress; no thermal conduction |
| No body forces | Gravity and electromagnetic forces neglected |
| Axisymmetry | No azimuthal velocity $$u_\theta$$; $$\partial / \partial \theta = 0$$ |
| Single-pressure cells | Air and product pressures linearly blended by $$\lambda$$; no sharp interface tracking |
| Eulerian fixed mesh | Material passes through cell faces; no mesh motion |
| Symmetry axis at $$r = 0$$ | Domain must include $$r = 0$$ for the reflective axis BC; off-axis problems use zero-gradient outflow there |
| Programmed burn | Detonation front is a kinematic sphere from the initiation point at speed $$D$$; reactive-flow effects (ignition delay, growth) are ignored |
| Constant $$\gamma = 1.4$$ for air | No temperature-dependent specific heats; valid up to a few thousand K |
