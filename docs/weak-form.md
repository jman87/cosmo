---
title: Weak Formulation
nav_order: 6
description: "Galerkin variational form of the axisymmetric Euler equations in cylindrical coordinates."
---

# Weak Formulation
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Overview

The finite-element method replaces the strong-form PDEs with an equivalent
**variational (weak) statement** that requires the residual to be orthogonal to a
suitable space of test functions.  This section derives the semi-discrete weak form
implemented in `solver/weak_form.py`.

The outcome is a system of ordinary differential equations (ODEs) in time:

$$
\mathbf{M}\,\dot{\mathbf{U}} = \mathbf{R}(\mathbf{U}),
$$

where $$\mathbf{M}$$ is the (lumped) mass matrix, $$\mathbf{U}$$ is the global DOF vector,
and $$\mathbf{R}$$ is the assembled spatial residual vector.  Time integration of this ODE
system is described in [Time Integration](time-integration).

---

## Cylindrical volume element

In 3-D cylindrical coordinates the volume element is $$dV = r\,dr\,d\theta\,dz$$.
Integrating out the azimuthal angle (all fields are $$\theta$$-independent) gives the
**effective 2-D measure**

$$
d\Omega = r\,dr\,dz.
$$

The factor $$r$$ multiplies every term in the weak form.  It vanishes on the symmetry axis
$$r = 0$$, which is handled numerically by a floor $$r_\text{safe} = \max(r, 10^{-15})$$.

---

## Weak form derivation

### General procedure

For a scalar conservation law

$$
\frac{\partial U}{\partial t} + \nabla_\text{cyl} \cdot \mathbf{F} = S,
$$

with cylindrical divergence $$\nabla_\text{cyl} \cdot \mathbf{F} = \frac{1}{r}\partial_r(rF_r) + \partial_z F_z$$:

1. Multiply by a test function $$v \in V_h$$ and integrate with measure $$r\,dr\,dz$$:

$$
\int_\Omega \frac{\partial U}{\partial t}\,v\,r\,dr\,dz
+ \int_\Omega \left[\frac{\partial(rF_r)}{\partial r} + r\frac{\partial F_z}{\partial z}\right] v\,dr\,dz
= \int_\Omega S\,v\,r\,dr\,dz.
$$

2. Integrate by parts on the flux terms (dropping the boundary integral, which imposes
   the zero-flux boundary condition — see [Boundary Conditions](boundary-conditions)):

$$
\int_\Omega \frac{\partial(rF_r)}{\partial r}\,v\,dr\,dz
= -\int_\Omega r\,F_r\,\frac{\partial v}{\partial r}\,dr\,dz + \underbrace{\oint_{\partial\Omega} r\,F_r\,v\,dz}_{=\,0}.
$$

3. Combine radial and axial contributions:

$$
\int_\Omega \frac{\partial U}{\partial t}\,v\,r\,dr\,dz
= \int_\Omega (F_r,\,F_z)\cdot\nabla v\,r\,dr\,dz + \int_\Omega S\,v\,r\,dr\,dz.
$$

The right-hand side is the **residual** $$R(v)$$.

### Continuity equation

With $$U = \rho$$, $$\mathbf{F}_\rho = (\rho u_r,\, \rho u_z)$$, and $$S = 0$$:

$$
\int_\Omega \dot\rho\,v\,r\,dr\,dz
= \int_\Omega \mathbf{F}_\rho \cdot \nabla v\,r\,dr\,dz.
$$

In code notation:

```python
R_rho = ufl.inner(F_rho, ufl.grad(v_rho)) * r_safe * dx
```

### Radial momentum equation

With $$U = \rho u_r$$, radial flux $$\mathbf{F}_{\rho u_r} = (\rho u_r^2 + p,\; \rho u_r u_z)$$,
and the geometric hoop-stress source $$S = p/r$$:

$$
\int_\Omega \dot{(\rho u_r)}\,v_r\,r\,dr\,dz
= \int_\Omega \mathbf{F}_{\rho u_r} \cdot \nabla v_r\,r\,dr\,dz
+ \underbrace{\int_\Omega \frac{p}{r}\,v_r\,r\,dr\,dz}_{\displaystyle=\,\int_\Omega p\,v_r\,dr\,dz}.
$$

The hoop-stress term $$\int_\Omega p\,v_r\,dr\,dz$$ has the $$r$$ factors cancel, leaving a
simple unweighted integral.  In code notation:

```python
R_rhou = (
    ufl.inner(ufl.as_tensor([F_rhou_r, F_rhou_z]), ufl.grad(v_rhou)) * r_safe * dx
    + p_eff * v_rhou[0] * dx
)
```

The positive sign is correct: $$+p/r$$ is the geometric source that drives radial expansion
of a cylindrical fluid element.

### Axial momentum equation

With $$\mathbf{F}_{\rho u_z} = (\rho u_r u_z,\; \rho u_z^2 + p)$$ and $$S = 0$$:

$$
\int_\Omega \dot{(\rho u_z)}\,v_z\,r\,dr\,dz
= \int_\Omega \mathbf{F}_{\rho u_z} \cdot \nabla v_z\,r\,dr\,dz.
$$

This term is the axial component of `R_rhou` (component index 1 of the vector test
function).

### Energy equation

With $$\mathbf{F}_{\rho E} = ((\rho E + p)\,u_r,\; (\rho E + p)\,u_z)$$ and $$S = 0$$:

$$
\int_\Omega \dot{(\rho E)}\,v\,r\,dr\,dz
= \int_\Omega \mathbf{F}_{\rho E} \cdot \nabla v\,r\,dr\,dz.
$$

In code notation:

```python
R_rhoe = ufl.inner(F_rhoe, ufl.grad(v_rhoe)) * r_safe * dx
```

---

## Effective pressure

The pressure used in all flux expressions is the **effective pressure** including the
artificial bulk viscosity term $$q$$:

$$
p_\text{eff} = p(\rho, e) + q(\rho, \mathbf{u}),
$$

where $$p$$ is the blended JWL/ideal-gas pressure and $$q$$ is the von Neumann–Richtmyer
viscosity (see [Shock Capturing](shock-capturing)).

---

## Flux tensors in full

For reference, the complete flux tensors used in the code are:

$$
\mathbf{F}_\rho = \begin{pmatrix}\rho u_r \\ \rho u_z\end{pmatrix},
$$

$$
\underline{\mathbf{F}}_{\rho\mathbf{u}} =
\begin{pmatrix}
  \mathbf{F}_{\rho u_r}^T \\[4pt]
  \mathbf{F}_{\rho u_z}^T
\end{pmatrix}
=
\begin{pmatrix}
  \rho u_r^2 + p_\text{eff} & \rho u_r u_z \\
  \rho u_r u_z & \rho u_z^2 + p_\text{eff}
\end{pmatrix}
= \rho\,\mathbf{u}\otimes\mathbf{u} + p_\text{eff}\,\mathbf{I},
$$

$$
\mathbf{F}_{\rho E} = \begin{pmatrix}(\rho E + p_\text{eff})\,u_r \\ (\rho E + p_\text{eff})\,u_z\end{pmatrix}.
$$

The momentum residual uses the double contraction
$$\underline{\mathbf{F}}_{\rho\mathbf{u}} : \nabla\mathbf{v} = F_{ij}\,\partial_j v_i$$,
which is `ufl.inner(tensor, grad(v))` in UFL.

---

## Numerical floors and stability

| Floor | Value | Purpose |
|-------|-------|---------|
| $$\rho_\text{floor}$$ | $$10^{-30}$$ | Prevent division by zero in $$u = (\rho u)/\rho$$ |
| $$e_\text{floor}$$ | $$10^{-30}$$ | Prevent negative sound speed in $$c^2 = \gamma(\gamma-1)e$$ |
| $$r_\text{floor}$$ | $$10^{-15}$$ | Prevent division by zero in $$r$$-weighted integrals on axis |
| $$p \geq 0$$ | — | Physical constraint; enforced via `ufl.max_value(p, 0)` |
