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

COSMO models a detonation in a fluid domain that is rotationally symmetric about the
$$z$$-axis.  It is therefore sufficient to solve on the 2-D meridional half-plane
$$\Omega \subset \{(r,z) \mid r \geq 0\}$$, where

* $$r$$ is the radial distance from the symmetry axis (horizontal in the mesh), and
* $$z$$ is the axial coordinate (vertical in the mesh).

The full 3-D solution is recovered by rotating the 2-D solution through $$2\pi$$ radians.
No azimuthal ($$\theta$$) derivatives appear because all fields are independent of $$\theta$$.

---

## Conserved variables

The solver advances four conserved quantities per spatial point:

$$
\mathbf{U} =
\begin{pmatrix}
  \rho \\
  \rho u_r \\
  \rho u_z \\
  \rho E
\end{pmatrix}
$$

| Symbol | Meaning | Units |
|--------|---------|-------|
| $$\rho$$ | Mass density | lbf·s²/in⁴ |
| $$\rho u_r$$ | Radial momentum density | lbf·s/in³ |
| $$\rho u_z$$ | Axial momentum density | lbf·s/in³ |
| $$\rho E$$ | Total energy density | lbf/in² |

where $$E = e + \tfrac{1}{2}(u_r^2 + u_z^2)$$ is the specific total energy, $$e$$ is the
specific internal energy (in²/s²), and $$u_r$$, $$u_z$$ are the velocity components (in/s).

---

## Strong conservative form

The inviscid, compressible Euler equations in 2-D axisymmetric coordinates are

$$
\frac{\partial \rho}{\partial t}
  + \frac{1}{r}\frac{\partial (r\,\rho u_r)}{\partial r}
  + \frac{\partial (\rho u_z)}{\partial z}
  = 0
\tag{mass}
$$

$$
\frac{\partial (\rho u_r)}{\partial t}
  + \frac{1}{r}\frac{\partial\bigl(r(\rho u_r^2 + p)\bigr)}{\partial r}
  + \frac{\partial (\rho u_r u_z)}{\partial z}
  - \frac{p}{r}
  = 0
\tag{$r$-momentum}
$$

$$
\frac{\partial (\rho u_z)}{\partial t}
  + \frac{1}{r}\frac{\partial (r\,\rho u_r u_z)}{\partial r}
  + \frac{\partial (\rho u_z^2 + p)}{\partial z}
  = 0
\tag{$z$-momentum}
$$

$$
\frac{\partial (\rho E)}{\partial t}
  + \frac{1}{r}\frac{\partial\bigl(r\,(\rho E + p)\,u_r\bigr)}{\partial r}
  + \frac{\partial\bigl((\rho E + p)\,u_z\bigr)}{\partial z}
  = 0
\tag{energy}
$$

---

## Compact notation

Defining the radial and axial flux vectors

$$
\mathbf{F}({\mathbf{U}}) =
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

and the geometric source vector

$$
\mathbf{S}(\mathbf{U}) =
\begin{pmatrix} 0 \\ p/r \\ 0 \\ 0 \end{pmatrix}
$$

the system is written as

$$
\frac{\partial \mathbf{U}}{\partial t}
  + \frac{1}{r}\frac{\partial (r\,\mathbf{F})}{\partial r}
  + \frac{\partial \mathbf{G}}{\partial z}
  = \mathbf{S}.
$$

The geometric source $$p/r$$ in the radial-momentum equation is the **hoop stress** (or
centripetal pressure) term.  It arises because the coordinate divergence $$\nabla \cdot \mathbf{F}$$
in cylindrical coordinates contains a $$1/r$$ metric factor: expanding
$$\frac{1}{r}\partial_r(r(\rho u_r^2 + p)) = \frac{1}{r}\partial_r(r\,\rho u_r^2) + \partial_r p + p/r$$.
The divergence form on the left-hand side absorbs the first two terms, leaving $$-p/r$$ as
an apparent source (equivalently $$+p/r$$ on the right-hand side).  Physically it represents
the net outward pressure force per unit volume on a cylindrical fluid element due to its
finite angular extent.

---

## Primitive variable recovery

The weak form and the EOS require *primitive* variables $$(\rho, u_r, u_z, p)$$ recovered
from the conserved state:

$$
u_r = \frac{\rho u_r}{\rho}, \qquad
u_z = \frac{\rho u_z}{\rho},
$$

$$
e = \frac{\rho E}{\rho} - \frac{1}{2}(u_r^2 + u_z^2),
$$

$$
p = p(\rho,\, e) \quad \text{via the equation of state.}
$$

Division by $$\rho$$ is protected by a numerical floor $$\rho \geq \rho_\text{floor} = 10^{-30}$$
to avoid NaN propagation in near-vacuum cells at the domain boundary.  Similarly, the
specific internal energy is clamped: $$e \geq e_\text{floor} = 10^{-30}$$.

---

## Isentropic sound speed

The CFL condition (see [Time Integration](time-integration)) requires a local wave-speed
estimate.  For a general equation of state the squared isentropic sound speed is

$$
c^2 = \left.\frac{\partial p}{\partial \rho}\right|_s
    = \frac{\partial p}{\partial \rho}\bigg|_e
      + \frac{p}{\rho^2} \frac{\partial p}{\partial e}\bigg|_\rho.
$$

For the **ideal-gas** EOS this simplifies to $$c^2 = \gamma(\gamma-1)e = \gamma p/\rho$$.
Inside the explosive charge the **JWL** sound speed is larger; the solver additionally
bounds the time step by the user-supplied detonation velocity $$D$$ whenever any material
indicator $$\lambda > 0$$ is present (see [Time Integration](time-integration)).

---

## Assumptions and limitations

| Assumption | Implication |
|------------|-------------|
| Inviscid flow | No viscous stress; no thermal conduction |
| No body forces | Gravity and magnetic forces neglected |
| Axisymmetry | No azimuthal velocity $$u_\theta$$; all $$\partial/\partial\theta = 0$$ |
| Single-material cells | Pressure blended by scalar burn fraction; no sharp interface tracking |
| Fixed Eulerian mesh | No mesh motion; material passes through cell faces (Eulerian description) |
| Symmetry axis at $$r = 0$$ | Domain must include $$r = 0$$; no off-axis problems |
