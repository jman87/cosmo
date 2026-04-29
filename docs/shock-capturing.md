---
title: Shock Capturing
nav_order: 7
description: "Von Neumann–Richtmyer artificial bulk viscosity for shock smearing."
---

# Shock Capturing
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## The challenge of shock waves in finite elements

The compressible Euler equations admit discontinuous solutions — shock waves — even from
smooth initial data.  A continuous Galerkin (CG) scheme applied to a shock-bearing problem
without stabilization produces **spurious post-shock oscillations** (Gibbs phenomenon).
These oscillations arise because CG basis functions cannot represent jump discontinuities
exactly; the expansion coefficients attempt to approximate the jump with high-frequency
ripples that can amplify to the point of crashing the simulation.

Several remedies exist:

| Approach | Description |
|----------|-------------|
| **Artificial viscosity (AV)** | Add a dissipative term that smears the shock over a few elements |
| **Upwind/Godunov flux** | Use a Riemann solver at element faces to impose upwinding |
| **Discontinuous Galerkin (DG)** | Allow discontinuities between elements; use numerical fluxes |
| **Flux-corrected transport (FCT)** | Blend high- and low-order schemes element-by-element |

COSMO uses the **von Neumann–Richtmyer artificial bulk viscosity** because it is
simple to implement as an additive pressure term in an existing CG formulation and has
a long track record in weapons physics codes.

---

## Von Neumann–Richtmyer viscosity

### Historical context

Von Neumann and Richtmyer (1950) introduced an artificial viscosity pressure $$q$$ to
spread shock waves over a small number of mesh cells in finite-difference calculations
of 1-D implosions.  Wilkins (1980) extended the method to multi-dimensional Lagrangian
codes.  The form used here follows Wilkins.

### Formulation

The effective pressure in all flux expressions is augmented by $$q$$:

$$
p_\text{eff} = p + q.
$$

The artificial pressure $$q$$ is active **only in compressive zones** and is defined as

$$
\boxed{
q = \rho\, h\, s \left( c_q\, h\, s + c_l\, c_s \right), \quad
s = \max(-\nabla\cdot\mathbf{u},\; 0)
}
$$

where:

| Symbol | Meaning | Default |
|--------|---------|---------|
| $$h$$ | Local cell diameter (`ufl.CellDiameter`) | — |
| $$s$$ | Positive part of velocity divergence (compression rate) | — |
| $$c_q$$ | Quadratic VNR coefficient | 1.5 |
| $$c_l$$ | Linear damping coefficient | 0.06 |
| $$c_s$$ | Local sound-speed estimate $$\sqrt{\gamma(\gamma-1)e}$$ | — |

The condition $$s = \max(-\nabla\cdot\mathbf{u}, 0)$$ ensures that $$q > 0$$ only when the
flow is compressive ($$\nabla\cdot\mathbf{u} < 0$$).  In expansion zones $$q = 0$$ and the
physical pressure is unchanged.

---

## Physical interpretation of each term

### Quadratic term: $$c_q\, \rho\, h^2\, s^2$$

This term is proportional to the square of the compression rate and dominates at **strong
shocks**.  It spreads a strong shock over approximately $$1/c_q$$ cells.  The $$h^2$$
dependence means the viscosity vanishes as the mesh is refined (for a fixed shock jump),
so the scheme converges to the correct Rankine–Hugoniot solution in the limit
$$h \to 0$$.

The dimensional argument: at a shock with velocity jump $$\Delta u$$ across a cell of size
$$h$$, the compression rate is $$s \sim \Delta u / h$$.  The pressure jump due to $$q$$ is
$$\sim \rho h^2 (\Delta u/h)^2 = \rho (\Delta u)^2$$, which matches the order of
magnitude of the physical shock-pressure jump.

### Linear term: $$c_l\, \rho\, h\, s\, c_s$$

This term is proportional to $$h$$ rather than $$h^2$$ and dominates at **weak shocks and
acoustic waves**.  Its primary role is to damp high-frequency post-shock oscillations
(ringing) that the quadratic term alone cannot suppress.  Setting $$c_l = 0$$ disables
the linear term and may increase post-shock oscillations for weaker shocks.

---

## Effect on the conservation laws

Adding $$q$$ to the pressure modifies the flux vectors:

$$
\mathbf{F}_{\rho u_r}^{(r)} = \rho u_r^2 + p + q, \qquad
\mathbf{F}_{\rho E}^{(r)} = (\rho E + p + q)\,u_r,
$$

and similarly for other components.  The energy flux modification ensures that the
artificial viscosity does work on the fluid, converting kinetic energy into internal
energy at the rate $$q\,\nabla\cdot\mathbf{u}$$.  This heating is unphysical but
negligible when the mesh is fine enough to resolve the shock.

Mass conservation is not affected by $$q$$ (density flux $$\rho\mathbf{u}$$ does not
contain pressure).

---

## Coefficient recommendations

| Scenario | $$c_q$$ | $$c_l$$ | Notes |
|----------|-------|-------|-------|
| Default (blast) | 1.5 | 0.06 | Good balance; ~3-cell shock width |
| Strong shocks | 1.5–2.0 | 0.05–0.10 | Wider smearing prevents carbuncle modes |
| Convergence study | 0.5–1.0 | 0.02 | Reduced artificial heating |
| Linear term off | 1.5 | 0 | May show post-shock ringing |

---

## Disabling shock capturing

Setting `shock_capturing: false` in the input JSON (or `{"enabled": false}`) returns
$$q \equiv 0$$ identically.  This is useful for smooth-flow verification cases where
artificial dissipation would pollute the solution, but it will cause unphysical oscillations
or blow-up for any problem containing shocks.
