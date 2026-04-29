---
title: Equations of State
nav_order: 3
description: "JWL EOS for TNT detonation products and ideal-gas EOS for ambient air."
---

# Equations of State
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Role of the EOS

The Euler equations contain five unknowns $$(\rho, u_r, u_z, E, p)$$ but only four
conservation laws.  The **equation of state** (EOS) closes the system by relating
pressure to the two thermodynamic state variables — density and specific internal energy:

$$
p = p(\rho,\, e).
$$

COSMO uses two EOS models:

1. **Jones–Wilkins–Lee (JWL)** for high-explosive (TNT) detonation products.
2. **Calorically perfect ideal gas** ($$\gamma = 1.4$$) for ambient air.

A scalar material indicator $$\lambda$$ (the burn fraction, $$\in [0,1]$$) linearly blends
the two pressures across partially-detonated cells:

$$
p = \lambda\, p_\text{JWL}(\rho, e) + (1 - \lambda)\, p_\text{air}(\rho, e).
\tag{blended}
$$

---

## JWL equation of state

### Physical basis

The JWL form was introduced by Lee, Finger, and Collins (1973) and is the industry
standard for high-explosive product gases.  It is calibrated to cylinder-expansion
test data and captures the two primary features of detonation products:

* a stiff, exponentially-varying contribution at high density (near the Chapman–Jouguet
  state), and
* a Grüneisen-type term $$\omega \rho e$$ that dominates at large expansion.

### Formula

$$
\boxed{
p_\text{JWL}(\rho, e) =
  A\!\left(1 - \frac{\omega}{\mathcal{R}_1\,\eta}\right)\exp\!\left(-\frac{\mathcal{R}_1}{\eta}\right)
+ B\!\left(1 - \frac{\omega}{\mathcal{R}_2\,\eta}\right)\exp\!\left(-\frac{\mathcal{R}_2}{\eta}\right)
+ \omega\,\rho\, e
}
$$

where

$$
\eta = \frac{\rho}{\rho_0}
$$

is the **compression ratio** (ratio of current density to the unreacted explosive density
$$\rho_0$$).  The first two terms account for the repulsive high-density behavior; the third
is the Mie–Grüneisen (thermal) term with Grüneisen parameter $$\omega$$.

### TNT parameters

The following constants are from Dobratz & Crawford (1985) / LLNL EOS Handbook, converted
from SI to imperial units using

$$
1\ \text{Pa} = 1.450\,377 \times 10^{-4}\ \text{psi},
\qquad
1\ \text{kg/m}^3 = 9.357 \times 10^{-8}\ \text{lbf}\cdot\text{s}^2/\text{in}^4,
\qquad
1\ \text{J/kg} = 1550.003\ \text{in}^2/\text{s}^2.
$$

| Parameter | SI value | Imperial value | Units (imperial) |
|-----------|----------|----------------|-----------------|
| $$A$$ | 373.77 GPa | 5.420 × 10⁷ | psi |
| $$B$$ | 3.747 GPa | 5.435 × 10⁵ | psi |
| $$\mathcal{R}_1$$ | 4.15 | 4.15 | — |
| $$\mathcal{R}_2$$ | 0.90 | 0.90 | — |
| $$\omega$$ | 0.35 | 0.35 | — |
| $$\rho_0$$ | 1630 kg/m³ | 1.526 × 10⁻⁴ | lbf·s²/in⁴ |
| $$p_\text{CJ}$$ | 21.0 GPa | 3.046 × 10⁶ | psi |
| $$e_\text{CJ}$$ | 6.0 MJ/kg | 9.300 × 10⁹ | in²/s² |

$$p_\text{CJ}$$ and $$e_\text{CJ}$$ are the Chapman–Jouguet detonation pressure and specific
energy; they set the initial thermodynamic state inside the charge (see
[Programmed Burn Model](burn-model)).

### JWL sound speed

The squared isentropic sound speed for JWL is

$$
c^2 = \frac{\partial p}{\partial \rho}\bigg|_e + \frac{p}{\rho^2}\frac{\partial p}{\partial e}\bigg|_\rho.
$$

The two partial derivatives are

$$
\frac{\partial p}{\partial e}\bigg|_\rho = \omega\,\rho,
$$

$$
\frac{\partial p}{\partial \rho}\bigg|_e =
  \frac{1}{\rho_0}\left[
    A\left(\frac{\mathcal{R}_1}{\eta} - \frac{\omega}{\mathcal{R}_1\,\eta^2}\right)e^{-\mathcal{R}_1/\eta}
  + B\left(\frac{\mathcal{R}_2}{\eta} - \frac{\omega}{\mathcal{R}_2\,\eta^2}\right)e^{-\mathcal{R}_2/\eta}
  \right].
$$

In practice, the CFL time step uses the ideal-gas sound speed formula for simplicity but
supplements it with the detonation velocity $$D$$ as a safety bound whenever explosive
products are present (see [Time Integration](time-integration)).

---

## Ideal-gas equation of state

### Formula

For a **calorically perfect ideal gas** the pressure is

$$
\boxed{p_\text{air}(\rho, e) = (\gamma - 1)\,\rho\, e}
$$

where $$\gamma = C_p/C_v = 1.4$$ for diatomic air.  This follows from the thermal equation
of state $$p = \rho R T$$ combined with $$e = c_v T = R T/(\gamma - 1)$$.

The squared isentropic sound speed simplifies to

$$
c^2 = \gamma(\gamma - 1)\, e = \frac{\gamma\, p}{\rho}.
$$

### Ambient air initial conditions

Standard sea-level atmosphere in imperial units:

| Quantity | Value | Units |
|----------|-------|-------|
| $$p_\text{atm}$$ | 14.696 | psi |
| $$\rho_\text{air}$$ | 1.146 × 10⁻⁷ | lbf·s²/in⁴ |
| $$e_\text{air} = p_\text{atm}/[(\gamma-1)\rho_\text{air}]$$ | 3.204 × 10⁸ | in²/s² |

---

## Pressure blending

Because COSMO uses a single-fluid Eulerian description, a cell may contain a
mixture of explosive products and air as the detonation front passes.  The effective
pressure is the linear blend

$$
p_\text{eff} = \lambda\, p_\text{JWL}(\rho, e) + (1 - \lambda)\, p_\text{air}(\rho, e),
$$

where $$\lambda \in [0,1]$$ is the programmed-burn fraction (see [Programmed Burn Model](burn-model)).
For fully-burned cells $$\lambda = 1$$ (JWL only); for pure air cells $$\lambda = 0$$
(ideal gas only).

The blended pressure is additionally clamped non-negative:

$$
p \leftarrow \max(p_\text{eff},\, 0)
$$

because negative (tensile) pressures are physically inadmissible for this problem and
would crash the explicit time integrator.

---

## UFL / NumPy dispatch

Both EOS classes expose a `pressure(rho, e)` method that works with plain NumPy arrays
(used during initial-condition setup and unit testing) and with UFL symbolic expressions
(used inside the DOLFINx weak form, where automatic differentiation generates the
element-level Jacobian).  The dispatch is handled by a small `_exp()` helper that calls
`ufl.exp()` when the argument is a UFL expression and `numpy.exp()` otherwise.
