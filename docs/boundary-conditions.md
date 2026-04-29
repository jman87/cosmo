---
title: Boundary Conditions
nav_order: 9
description: "Ghost-cell mirror reflection on the symmetry axis and the rigid ground, zero-gradient outflow elsewhere, and inter-rank halo seams."
---

# Boundary Conditions
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Ghost-cell formulation

The cell-centred finite-volume framework imposes boundary conditions
through values written into the **ghost cells** that border the
interior on every side. Each rank owns $$n_g = 2$$ ghost layers
([Mesh and Grid](mesh-discretization)). The MUSCL reconstruction and
HLLC flux at the first/last interior face read from the ghost cells in
exactly the same way they read from interior cells, so a BC is enforced
implicitly by populating the ghost values *before* the residual is
assembled.

Each ghost layer is filled in one of two ways:

1. By the **physical-boundary BC routine** `apply_bcs!` if the layer
   lies on the rank's exterior edge of the global domain;
2. By the **MPI halo-exchange** `exchange_halos!` if the layer abuts an
   internal partition seam.

The two cases are mutually exclusive on every ghost layer. A rank that
sits on the global $$r = 0$$ edge applies the symmetry-axis BC on its
left ghost layer; a rank deeper in the interior receives that ghost
layer from its neighbour. Edge tests (`at_r_min_edge`, `at_r_max_edge`,
etc.) on the rank's Cartesian coordinates select between the two paths.

The implementation lives in
[`src/boundary.jl`](https://github.com/jman87/cosmo/blob/main/src/boundary.jl)
(physical-edge BCs) and
[`src/parallel_halo.jl`](https://github.com/jman87/cosmo/blob/main/src/parallel_halo.jl)
(MPI halo exchange).

---

## Symmetry axis at $$r = 0$$

The axisymmetric formulation requires that the velocity field be
mirror-symmetric about the $$z$$-axis: at $$r = 0$$ the radial velocity
must vanish, and density / axial-velocity / total-energy fields must
have zero radial gradient. The discrete analogue is a **mirror
reflection** of the interior cells across the axis, with the radial
momentum negated:

$$
\rho^{\,g}_{i,j} \;=\; \rho^{\,\text{src}}_{i,j},
$$

$$
(\rho u_r)^{\,g}_{i,j} \;=\; -\,(\rho u_r)^{\,\text{src}}_{i,j},
$$

$$
(\rho u_z)^{\,g}_{i,j} \;=\; (\rho u_z)^{\,\text{src}}_{i,j},
$$

$$
(\rho E)^{\,g}_{i,j} \;=\; (\rho E)^{\,\text{src}}_{i,j},
$$

where the source-cell index for ghost cell $$g \in \{1,\,2\}$$ is the
mirror image across the boundary face:

$$
\text{src} \;=\; 2\,n_g + 1 - g.
$$

For $$n_g = 2$$ this gives $$g = 1 \;\to\; \text{src} = 4$$ and
$$g = 2 \;\to\; \text{src} = 3$$, i.e. the second-nearest interior cell
maps into the first ghost cell, etc. The mirror reflection produces a
state on the ghost side of the axis that is the exact reflection of the
state on the interior side, so the HLLC face flux at the axis face
recovers $$u_r = 0$$ to machine precision and the radial-momentum flux
through the axis face is exactly zero (modulo signs in symmetric pairs
that cancel).

The reflective BC is **active only when $$|r_\text{min}| < 10^{-12}$$**,
i.e. the user has placed the inner edge of the domain on the symmetry
axis. For an off-axis problem ($$r_\text{min} > 0$$) the inner edge
becomes a zero-gradient outflow boundary instead (see below). The
selection happens once at startup based on the JSON `domain.r_min`.

The same mirror+sign-flip rule is applied automatically by virtue of
the choice $$\text{src} = 2 n_g + 1 - g$$, regardless of the number of
ghost layers, so the BC scales naturally to any future increase in
$$n_g$$.

---

## Rigid ground at $$z = 0$$

When $$z_\text{min} = 0$$ the lower boundary represents a **rigid
horizontal symmetry plane** (the standard ground-burst assumption: the
solver evolves the upper half-space and the ground is modelled as a
perfect reflector). The required symmetry is

$$
u_z\big|_{z=0} \;=\; 0,
\qquad
\partial_z \rho \;=\;
\partial_z u_r \;=\;
\partial_z (\rho E) \;=\; 0.
$$

The discrete BC is the analogue of the axis case with the role of $$r$$
and $$z$$ swapped: mirror the interior across the ground plane and
negate the *axial* momentum component. As above, the BC is **active
only when $$|z_\text{min}| < 10^{-12}$$**; otherwise the lower edge is a
zero-gradient outflow.

For a free-air burst (no ground), set $$z_\text{min}$$ to a non-zero
value (typically a large negative number that places the charge well
above the lower edge) so that the ground-plane reflection does not
activate.

---

## Outer boundaries: $$r = r_\text{max}$$ and $$z = z_\text{max}$$

The two outer edges are always treated as **zero-gradient (Neumann)
outflow** boundaries:

$$
\rho^{\,g}_{i,j} \;=\; \rho^{\,\text{last interior}}_{i,j},
$$

$$
(\rho u_r)^{\,g}_{i,j} \;=\; (\rho u_r)^{\,\text{last interior}}_{i,j},
$$

$$
(\rho u_z)^{\,g}_{i,j} \;=\; (\rho u_z)^{\,\text{last interior}}_{i,j},
$$

$$
(\rho E)^{\,g}_{i,j} \;=\; (\rho E)^{\,\text{last interior}}_{i,j}.
$$

This is the simplest **non-reflecting** BC: each ghost cell mirrors the
state of the **last interior** cell on the same row/column, so any wave
travelling outward sees no reflecting jump. It is a 1st-order
approximation of a true characteristic outflow and reflects a small
fraction of the incident wave back into the domain (most visibly for
weak acoustic waves arriving nearly normal to the boundary). For blast
problems whose outer-domain shock arrival times are well past the time
of physical interest this is acceptable; for runs where the outer wave
does drive analysis, the user must extend the domain so the wave has
not yet reached it.

The outer-radial and outer-axial outflow BCs are **always active**,
independent of the corresponding inner-edge setting, because there is
no symmetric counterpart on the outside of the domain. There is no
"reflective outer wall" mode in the current implementation; if a hard
wall is needed (e.g. a confined explosion), the user must implement it.

---

## Off-axis inner-radial boundary

If the user places $$r_\text{min} > 0$$ (i.e. an annular sub-domain
that excludes the axis), the inner-radial edge is treated as **zero-
gradient outflow** with the same mirror-source rule used at $$r =
r_\text{max}$$. The mirror reflection at the axis is *not* applied in
this case — the symmetry-axis assumption is invalid for off-axis
problems and applying it would impose an unphysical reflective
constraint.

This option is rarely used in practice for blast applications but is
exposed for verification cases (e.g. comparing against a 1-D radial
similarity solution evaluated on $$[r_0,\,r_\text{max}]$$ with $$r_0 >
0$$).

---

## Passive-scalar BC

The passive scalar $$\rho Y$$ obeys a strictly hyperbolic transport
equation with no source term. Both the mirror-reflection BC at the
axis / ground and the zero-gradient outflow BC reduce to the same
scalar copy operation:

$$
(\rho Y)^{\,g}_{i,j} \;=\; (\rho Y)^{\,\text{src}}_{i,j},
$$

with `src` = mirror-image cell at axis/ground edges and `src` = last
interior cell at outflow edges. No sign flip is applied because the
scalar carries no signed direction. The implementation is
`apply_scalar_bcs!` in `src/boundary.jl`.

---

## Inter-rank halo seams

A ghost layer that abuts another rank is **never written by
`apply_bcs!`** — it is filled by the halo-exchange routine
`exchange_halos!` (or `exchange_halo_scalar!` for the burn fraction and
the passive scalar). The same code path applies the physical BC on
every edge that maps to `MPI.PROC_NULL` and the halo exchange on every
edge that does not, so the per-rank residual sees a consistent ghost
state regardless of where the rank sits in the global domain.

The two-step "exchange first, then physical BC" ordering used in
`rhs!` (`src/timeint.jl`) is critical: a rank in the bottom-left corner
of the global domain has both the symmetry-axis and ground-plane
physical edges *and* halo seams with neighbours on the right and
above, so both routines must run in sequence to fill the four ghost
slabs correctly.

The MPI side of the halo exchange is described in
[MPI Parallelisation](parallelisation).

---

## Summary

| Edge | Default condition | Activates when | Effect |
|------|-------------------|----------------|--------|
| $$r = r_\text{min}$$ | Mirror reflection (symmetry axis) | $$\|r_\text{min}\| < 10^{-12}$$ | $$u_r$$ flipped; other components copied |
| $$r = r_\text{min}$$ | Zero-gradient outflow | $$r_\text{min} > 0$$ | All components copied from last interior cell |
| $$r = r_\text{max}$$ | Zero-gradient outflow | always | All components copied from last interior cell |
| $$z = z_\text{min}$$ | Mirror reflection (rigid ground) | $$\|z_\text{min}\| < 10^{-12}$$ | $$u_z$$ flipped; other components copied |
| $$z = z_\text{min}$$ | Zero-gradient outflow | $$z_\text{min} \neq 0$$ | All components copied from last interior cell |
| $$z = z_\text{max}$$ | Zero-gradient outflow | always | All components copied from last interior cell |
| Internal seam | Halo exchange | rank has a non-NULL neighbour | Ghost slab copied from neighbour rank's interior |
