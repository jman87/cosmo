---
title: Boundary Conditions
nav_order: 9
description: "Symmetry axis, ground plane, and outer-wall boundary conditions."
---

# Boundary Conditions
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Weak-form boundary integrals

When the flux divergence is integrated by parts in the weak form, a surface integral
appears on the domain boundary $$\partial\Omega$$:

$$
\oint_{\partial\Omega} \mathbf{F}\cdot\hat{\mathbf{n}}\; v\; r\, d\ell,
$$

where $$\hat{\mathbf{n}}$$ is the outward unit normal and $$d\ell$$ is the arc-length element.
**COSMO drops this integral entirely**, which is equivalent to setting the normal
flux to zero on every boundary:

$$
\mathbf{F}\cdot\hat{\mathbf{n}} = 0 \quad \text{on } \partial\Omega.
$$

This is the **impermeable solid-wall** (zero-normal-flux) condition.  It is the natural
Neumann BC for the CG formulation and requires no explicit imposition.

---

## Symmetry axis ($$r = 0$$)

### Physical requirement

The axisymmetric formulation requires that the flow field be mirror-symmetric about the
$$z$$-axis.  At $$r = 0$$ the radial velocity must vanish:

$$
u_r\big|_{r=0} = 0.
$$

The zero-flux BC from the dropped boundary integral enforces zero normal flux $$\rho u_r$$
at $$r = 0$$.  However, accumulation of floating-point errors in the explicit update can
cause a small non-zero $$\rho u_r$$ to develop there.  The solver therefore **explicitly
zeroes** the radial momentum DOFs on the axis after every RK stage:

$$
(\rho u_r)_i = 0 \quad \forall\, \text{DOF } i \text{ with } r_i = 0.
$$

This is implemented in `_apply_axis_bc` by computing the flat indices of the $$r$$-component
DOFs at $$r = 0$$ once at startup (using `_find_symmetry_dofs`) and zeroing them each stage.

### Regularity on the axis

All other fields ($$\rho$$, $$\rho u_z$$, $$\rho E$$) are smooth and bounded on the axis by
physics; no special treatment is required beyond the floor on $$r_\text{safe}$$ used in
the weak form.  The lumped mass row-sum is zero at axis nodes (the cylindrical measure
$$r\,dr\,dz$$ vanishes there), but these entries are floored to prevent division by zero
while maintaining the expected near-zero updates (see [Time Integration](time-integration)).

---

## Ground-plane symmetry ($$z = 0$$)

When `domain.z_min = 0` (the default), the lower boundary of the domain represents a
**symmetry plane** (e.g., a ground burst with half-space modelling).  The axial velocity
must vanish there:

$$
u_z\big|_{z=0} = 0.
$$

The solver zeroes the axial momentum DOFs at $$z = 0$$ after every RK stage:

$$
(\rho u_z)_i = 0 \quad \forall\, \text{DOF } i \text{ with } z_i = 0.
$$

When `domain.z_min < 0` the lower boundary is treated as a standard solid wall (zero
normal flux only); the axial-velocity constraint is not applied there.

---

## Outer walls ($$r = r_\text{max}$$ and $$z = z_\text{max}$$)

The outer boundaries receive the same **solid-wall (zero-flux) condition** as all other
boundaries.  When a shock wave reaches the outer wall it reflects back into the domain —
which is **physically incorrect** for a free-air blast.  This is acceptable for short
simulations where the shock has not yet reached the outer boundary.

For longer simulations, the user should either:

* Extend the domain so the shock does not reach the boundary within the simulation time,
  or
* Implement a **non-reflecting outflow condition** (e.g., a characteristics-based
  absorbing BC or a Perfectly Matched Layer).

The domain should be sized so that the scaled standoff distance at the outer boundary
exceeds the range of interest by at least 10%.  The example input files are set up this
way.

---

## Summary

| Boundary | Condition | How enforced |
|----------|-----------|--------------|
| $$r = 0$$ (symmetry axis) | $$u_r = 0$$; zero normal flux | Dropped BC integral + explicit DOF zeroing |
| $$z = 0$$ (when $$z_\text{min} = 0$$) | $$u_z = 0$$; zero normal flux | Dropped BC integral + explicit DOF zeroing |
| $$r = r_\text{max}$$ (outer radial wall) | Zero normal flux; wall reflection | Dropped BC integral |
| $$z = z_\text{max}$$ (outer axial wall) | Zero normal flux; wall reflection | Dropped BC integral |
