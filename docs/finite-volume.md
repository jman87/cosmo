---
title: Finite-Volume Discretisation
nav_order: 6
description: "Cell-volume averages, the discrete divergence theorem in axisymmetric coordinates, residual assembly, and the geometric source term."
---

# Finite-Volume Discretisation
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Cell-volume averaged state

The cell-centred finite-volume method stores the **average** of the
conservative state over each cell rather than a point value:

$$
\overline{\mathbf{U}}_{ij}(t) \;=\;
\frac{1}{V_{ij}} \int_{C_{ij}} \mathbf{U}(r, z, t)\,r\,dr\,dz,
\qquad
V_{ij} \;=\; r_i\,\Delta r\,\Delta z.
$$

In the absence of ambiguity the overbar is dropped and
$$\mathbf{U}_{ij}$$ is understood to denote the cell average.

The $$r$$-weight in the integrand reflects the cylindrical volume
element. Because $$2\pi$$ cancels uniformly across all balance equations
once the azimuthal angle is integrated out, the solver works with the
*annular* cell volume $$V_{ij} = r_i\,\Delta r\,\Delta z$$ (no $$2\pi$$).

---

## Integral conservation law

Integrating the strong-conservative form
([Governing Equations](governing-equations))

$$
\frac{\partial \mathbf{U}}{\partial t}
\;+\; \frac{1}{r}\frac{\partial(r\,\mathbf{F})}{\partial r}
\;+\; \frac{\partial \mathbf{G}}{\partial z}
\;=\; \mathbf{S}_p
$$

over a single cell $$C_{ij}$$ and applying the divergence theorem in
cylindrical coordinates yields the **integral conservation law**:

$$
\frac{d}{dt}\!\int_{C_{ij}}\!\mathbf{U}\,r\,dr\,dz
\;+\; \oint_{\partial C_{ij}}\!\!\bigl(\mathbf{F}\,n_r + \mathbf{G}\,n_z\bigr)\,r\,d\ell
\;=\; \int_{C_{ij}}\!\mathbf{S}_p\,r\,dr\,dz.
$$

The boundary $$\partial C_{ij}$$ is the union of four faces; the unit
normal $$\hat{\mathbf{n}} = (n_r, n_z)$$ is purely radial on the
$$r$$-faces and purely axial on the $$z$$-faces.

Evaluating each face integral with the midpoint rule (which is
2nd-order accurate and natural for cell-centred schemes) gives

$$
\int_{r\text{-face at } r_{i \pm 1/2}}\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!
\mathbf{F}\,r\,d\ell
\;\approx\; \mathbf{F}^*_{i \pm 1/2,\,j}\,r_{i \pm 1/2}\,\Delta z,
$$

$$
\int_{z\text{-face at } z_{j \pm 1/2}}\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!\!
\mathbf{G}\,r\,d\ell
\;\approx\; \mathbf{G}^*_{i,\,j \pm 1/2}\,r_i\,\Delta r,
$$

where $$\mathbf{F}^*$$ and $$\mathbf{G}^*$$ are **numerical fluxes**
evaluated by the Riemann solver from the reconstructed left/right states
on each face — see [Riemann Solver and Reconstruction](riemann-reconstruction).

---

## Semi-discrete cell-update equation

Dividing by the cell volume $$V_{ij} = r_i\,\Delta r\,\Delta z$$ produces
the semi-discrete ODE for each cell average:

$$
\boxed{\;\;
\frac{d \mathbf{U}_{ij}}{dt}
\;=\;
- \frac{1}{V_{ij}}
  \Bigl[
    \bigl(\mathbf{F}^*_{i+1/2,j}\,r_{i+1/2} \;-\; \mathbf{F}^*_{i-1/2,j}\,r_{i-1/2}\bigr)\,\Delta z
    \;+\;
    \bigl(\mathbf{G}^*_{i,j+1/2} \;-\; \mathbf{G}^*_{i,j-1/2}\bigr)\,r_i\,\Delta r
  \Bigr]
\;+\; \mathbf{S}_{p,\,ij}
\;\;}
\tag{FV}
$$

with the cell-centred geometric pressure source

$$
\mathbf{S}_{p,\,ij} \;=\; \bigl(0,\;\; p_{ij}/r_i,\;\; 0,\;\; 0\bigr)^{T}.
$$

The right-hand side of (FV) is the **spatial residual**
$$\mathbf{L}(\mathbf{U})_{ij}$$ that the explicit time integrator
([Time Integration](time-integration)) advances. In code this assembly
appears verbatim in `compute_residual!` in
[`src/timeint.jl`](https://github.com/jman87/cosmo/blob/main/src/timeint.jl).

The form of (FV) makes three desirable properties manifest:

1. **Local conservation**: the flux $$\mathbf{F}^*_{i+1/2,\,j}$$ leaves
   cell $$(i, j)$$ and enters cell $$(i+1, j)$$ with opposite sign, so
   any quantity governed by a homogeneous balance law (mass, total
   energy, axial momentum) is conserved exactly to floating-point
   round-off across an interior face. Radial momentum has the geometric
   source term and is conserved only globally up to the source integral,
   which is the correct behaviour in cylindrical coordinates.
2. **Telescoping**: when (FV) is summed over a strip of cells along
   either direction the interior fluxes cancel pairwise, leaving only
   boundary contributions — an algebraic statement of the divergence
   theorem.
3. **Bounded by Riemann physics**: each face flux is the solution to a
   Riemann problem in the face-normal direction, so the cell update
   inherits the wave-propagation structure of the Euler system.

---

## Annular face areas

The radial flux $$\mathbf{F}^*$$ is multiplied by the **face radius**
$$r_{i \pm 1/2}$$, **not** the cell-centre radius. This choice is what
makes the discretisation faithful to the cylindrical metric:

* The continuum integrand $$r\,\mathbf{F}$$ in the $$r$$-direction surface
  integral evaluates the *radius at the face* — using $$r_i$$ here would
  introduce an $$\mathcal{O}(\Delta r^2)$$ inconsistency in the
  divergence and break a uniform-flow steady state at finite $$r$$.
* The $$z$$-direction integrand $$\mathbf{G}$$ has *no* $$r$$-derivative
  in the metric, so a single $$r_i$$ representative of the cell suffices
  for the axial-face integral.

Numerically, on the symmetry axis the leftmost interior face has
$$r_{n_g+1/2} = 0$$, which automatically zeroes the radial-flux
contribution there. The cell volume $$V_{n_g+1,\,j} = r_{n_g+1}\,\Delta r\,\Delta z$$
remains finite because the cell *centre* is at the half-step
$$r_{n_g+1} = \Delta r / 2 > 0$$. The discretisation is therefore
geometrically well-defined right up to the symmetry axis, with no need
for the $$\max(r,\,\varepsilon)$$ floor used in the Python sibling
solver's weak form.

---

## Geometric pressure source

Both the divergence form (FV) and the cell-centred source
$$\mathbf{S}_{p,\,ij} = (0,\,p_{ij}/r_i,\,0,\,0)^{T}$$ are needed because
the radial-momentum balance in cylindrical coordinates carries a
metric-induced source term:

$$
\frac{1}{r}\frac{\partial}{\partial r}\!\bigl[r\,(\rho u_r^2 + p)\bigr]
\;=\;
\frac{1}{r}\frac{\partial(r\,\rho u_r^2)}{\partial r}
\;+\; \frac{\partial p}{\partial r}
\;+\; \frac{p}{r}.
$$

The divergence form on the LHS absorbs the first two terms;
$$+p/r$$ remains as the **hoop-stress source** on the RHS. Physically it
is the net outward pressure force per unit volume on a wedge-shaped fluid
element. Numerically it must be discretised with the cell-centred
pressure $$p_{ij}$$ and the cell-centre radius $$r_i$$:

$$
S_{p,\,ij} \;=\; \frac{p_{ij}}{r_i}, \qquad p_{ij} = p\!\left(\rho_{ij},\, e_{ij};\, \lambda_{ij}\right).
$$

Other components of $$\mathbf{S}_p$$ are zero. The treatment is identical
to the textbook "balance method" for cylindrical Euler systems
(e.g. LeVeque 2002 §17.2; Toro 2009 §8.5).

---

## Reconstruction and numerical flux

The numerical flux $$\mathbf{F}^*_{i+1/2,\,j}$$ on the radial face
between cells $$(i, j)$$ and $$(i+1, j)$$ is obtained by

1. **MUSCL reconstruction** of one-sided face states $$\mathbf{U}^L,\,
   \mathbf{U}^R$$ from the four-cell stencil
   $$(\mathbf{U}_{i-1,j},\,\mathbf{U}_{i,j},\,\mathbf{U}_{i+1,j},\,\mathbf{U}_{i+2,j})$$
   with the **minmod** slope limiter;
2. **HLLC** approximate Riemann-solver evaluation
   $$\mathbf{F}^* = \mathcal{F}_\text{HLLC}(\mathbf{U}^L,\,\mathbf{U}^R,\,\lambda^L,\,\lambda^R,\,\hat{r})$$
   using the local cell-centred reaction-progress values
   $$\lambda^L,\lambda^R$$ to evaluate the EOS on each side;
3. The same procedure on $$z$$-faces, exchanging the role of the normal
   and tangential velocity components.

Both steps are documented at length in
[Riemann Solver and Reconstruction](riemann-reconstruction).

---

## Passive-scalar discretisation

The passive scalar $$\rho Y$$ ([Programmed Burn Model](burn-model)) obeys
the same flux-divergence form as the mass equation, so its semi-discrete
update is

$$
\frac{d (\rho Y)_{ij}}{dt}
\;=\;
- \frac{1}{V_{ij}}
  \Bigl[
    \bigl(\Phi^*_{i+1/2,j}\,r_{i+1/2} \;-\; \Phi^*_{i-1/2,j}\,r_{i-1/2}\bigr)\,\Delta z
    \;+\;
    \bigl(\Psi^*_{i,j+1/2} \;-\; \Psi^*_{i,j-1/2}\bigr)\,r_i\,\Delta r
  \Bigr],
$$

where the face-mass flux $$F^{(\rho)*}$$ from the HLLC solve is reused as
the upwind selector:

$$
\Phi^*_{i+1/2,j} \;=\;
\begin{cases}
  F^{(\rho)*}_{i+1/2,j}\,Y_{i,j},   & F^{(\rho)*}_{i+1/2,j} \;\geq\; 0, \\[2pt]
  F^{(\rho)*}_{i+1/2,j}\,Y_{i+1,j}, & F^{(\rho)*}_{i+1/2,j} \;<\; 0,
\end{cases}
$$

with $$Y_{i,j} = (\rho Y)_{ij} / \max(\rho_{ij},\,\rho_\text{floor})$$
and the analogous expression for $$\Psi^*$$ on $$z$$-faces. This is
1st-order donor-cell upwinding in conservative form. Reusing the HLLC
mass flux automatically consistency-couples the scalar transport to the
Euler update — the passive scalar advances at exactly the speed of the
mass flux that the Riemann solver produces.

---

## Order of accuracy

| Component | Local truncation error |
|-----------|------------------------|
| Cell-centred storage | $$\mathcal{O}(\Delta x^2)$$ for smooth flow |
| Midpoint-rule face integration | $$\mathcal{O}(\Delta x^2)$$ |
| MUSCL reconstruction (smooth) | $$\mathcal{O}(\Delta x^2)$$ |
| MUSCL at extrema (limiter active) | $$\mathcal{O}(\Delta x)$$ |
| HLLC numerical flux | exact for isolated waves; smooth-flow consistency $$\mathcal{O}(\Delta x^2)$$ |
| Donor-cell scalar advection | $$\mathcal{O}(\Delta x)$$ |
| SSP-RK2 time stepper | $$\mathcal{O}(\Delta t^2)$$ in smooth flow, $$\mathcal{O}(\Delta t)$$ at shocks |

The combined scheme is therefore formally **2nd-order accurate in
smooth regions** and degrades to 1st order at shocks and extrema, which
is the expected and desirable behaviour for a TVD finite-volume method.
The scalar tag is intentionally only 1st order — its sole purpose is
visualisation and material identification, and a more accurate
reconstruction would not change any physical quantity in the solution.
