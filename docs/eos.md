---
title: Equations of State
nav_order: 3
description: "JWL EOS for TNT detonation products and ideal-gas EOS for ambient air; reaction-progress pressure blend."
---

# Equations of State
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Role of the EOS

The compressible Euler equations contain five unknowns
$$(\rho,\, u_r,\, u_z,\, E,\, p)$$ but only four conservation laws. The
**equation of state** (EOS) closes the system by relating pressure to two
thermodynamic state variables — density and specific internal energy:

$$
p = p(\rho,\, e).
$$

COSMO uses two EOS models in tandem:

1. **Jones–Wilkins–Lee (JWL)** for high-explosive (TNT) detonation products.
2. **Calorically perfect ideal gas** ($$\gamma = 1.4$$) for ambient air.

A scalar reaction-progress variable $$\lambda \in [0, 1]$$ — the
**burn fraction** — linearly blends the two pressures across cells whose
state lies between unreacted explosive ($$\lambda = 0$$) and fully-reacted
detonation products ($$\lambda = 1$$):

$$
p(\rho, e;\, \lambda)
\;=\; (1 - \lambda)\,p_\text{air}(\rho, e)
\;+\; \lambda\,p_\text{JWL}(\rho, e).
\tag{blend}
$$

The implementation of the blend lives in `mix_pressure` in
[`src/eos.jl`](https://github.com/jman87/cosmo/blob/main/src/eos.jl).

---

## Ideal-gas equation of state

### Formula

For a calorically perfect ideal gas the pressure is

$$
\boxed{\,p_\text{air}(\rho,\, e) \;=\; (\gamma - 1)\,\rho\, e\,}
$$

with $$\gamma = C_p / C_v = 1.4$$ for diatomic air. This follows from the
thermal EOS $$p = \rho R T$$ combined with the caloric relation
$$e = c_v T = R T / (\gamma - 1)$$.

The squared isentropic sound speed simplifies to

$$
c_\text{air}^2 \;=\; \gamma\,(\gamma - 1)\,e \;=\; \frac{\gamma\,p}{\rho}.
$$

Both expressions appear in the CFL estimator (see
[Time Integration](time-integration)) and, indirectly, in the HLLC wave-speed
estimates (see [Riemann Solver and Reconstruction](riemann-reconstruction)).

### Ambient air at sea level

The standard sea-level atmosphere in imperial units is

| Quantity | Value | Units |
|----------|-------|-------|
| $$p_\text{atm}$$ | $$14.696$$ | psi |
| $$\rho_\text{air}$$ | $$1.146 \times 10^{-7}$$ | $$\text{lbf}\cdot\text{s}^2/\text{in}^4$$ |
| $$e_\text{air} = p_\text{atm} / [(\gamma-1)\rho_\text{air}]$$ | $$3.204 \times 10^8$$ | $$\text{in}^2/\text{s}^2$$ |
| $$c_\text{air} = \sqrt{\gamma p_\text{atm} / \rho_\text{air}}$$ | $$\approx 1.34 \times 10^4$$ | $$\text{in}/\text{s}$$ |

These values populate every air cell at $$t = 0$$ (the velocity is zero, so
$$E_\text{air} = e_\text{air}$$).

---

## JWL equation of state

### Physical basis

The Jones–Wilkins–Lee form was introduced by Lee, Finger and Collins (1973)
and is the industry-standard EOS for high-explosive detonation products.
Calibrated to cylinder-expansion tests, it captures the two dominant
features of expanding products:

* a stiff, exponentially-varying contribution active at high density (near
  the Chapman–Jouguet state and above);
* a Mie–Grüneisen term $$\omega \rho e$$ that dominates at large expansion
  ratios.

### Formula

$$
\boxed{
p_\text{JWL}(\rho, e)
\;=\;
A\!\left(1 - \frac{\omega}{\mathcal{R}_1\,\eta}\right)\exp\!\left(-\frac{\mathcal{R}_1}{\eta}\right)
\;+\;
B\!\left(1 - \frac{\omega}{\mathcal{R}_2\,\eta}\right)\exp\!\left(-\frac{\mathcal{R}_2}{\eta}\right)
\;+\;
\omega\,\rho\,e
}
$$

where

$$
\eta \;=\; \frac{\rho}{\rho_0}
$$

is the **compression ratio** relative to the unreacted-explosive density
$$\rho_0$$. The first two exponential terms model the cold-curve (repulsive)
behaviour at high compression; the third is the Grüneisen (thermal) term
with Grüneisen coefficient $$\omega$$. The implementation lives in
`pressure_jwl` in `src/eos.jl`.

### TNT parameters

The constants below are the Dobratz & Crawford (1985) / LLNL EOS Handbook
values for TNT, converted from SI to imperial inside the `TNT = let ... end`
block in `src/eos.jl`:

| Parameter | SI value | Imperial value | Units (imperial) |
|-----------|----------|----------------|------------------|
| $$A$$ | $$373.77$$ GPa | $$5.420 \times 10^{7}$$ | psi |
| $$B$$ | $$3.747$$ GPa | $$5.435 \times 10^{5}$$ | psi |
| $$\mathcal{R}_1$$ | 4.15 | 4.15 | — |
| $$\mathcal{R}_2$$ | 0.90 | 0.90 | — |
| $$\omega$$ | 0.35 | 0.35 | — |
| $$\rho_0$$ | $$1630$$ kg/m³ | $$1.526 \times 10^{-4}$$ | $$\text{lbf}\cdot\text{s}^2/\text{in}^4$$ |
| $$p_\text{CJ}$$ | $$21.0$$ GPa | $$3.046 \times 10^{6}$$ | psi |
| $$e_\text{CJ}$$ | $$6.0$$ MJ/kg | $$9.300 \times 10^{9}$$ | $$\text{in}^2/\text{s}^2$$ |
| $$D$$ (CJ velocity) | $$6930$$ m/s | $$2.728 \times 10^{5}$$ | $$\text{in}/\text{s}$$ |

$$p_\text{CJ}$$ and $$e_\text{CJ}$$ are the Chapman–Jouguet detonation
pressure and specific internal energy of TNT; the latter sets the
in-charge thermodynamic state in the programmed-burn IC (see
[Programmed Burn Model](burn-model)). $$D$$ is the CJ detonation velocity
and controls the kinematics of the burn front.

The unit conversions used (constants `PA_TO_PSI`, `KG_M3_TO_LBF_S2_IN4`,
`J_KG_TO_IN2_S2` in code) are

$$
1\ \text{Pa} = 1.450\,377 \times 10^{-4}\ \text{psi},
$$

$$
1\ \text{kg/m}^3 = 9.357 \times 10^{-8}\ \text{lbf}\cdot\text{s}^2/\text{in}^4,
$$

$$
1\ \text{J/kg} = 1550.003\ \text{in}^2/\text{s}^2.
$$

### JWL sound speed used in the solver

The exact squared isentropic sound speed of JWL is

$$
c_\text{JWL}^2 \;=\;
\left.\frac{\partial p}{\partial \rho}\right|_e
\;+\; \frac{p}{\rho^2}\,\left.\frac{\partial p}{\partial e}\right|_\rho,
$$

with

$$
\left.\frac{\partial p}{\partial e}\right|_\rho \;=\; \omega\,\rho,
$$

$$
\left.\frac{\partial p}{\partial \rho}\right|_e
\;=\; \frac{1}{\rho_0}\Bigg[\,
A\!\left(\frac{\mathcal{R}_1}{\eta} - \frac{\omega}{\mathcal{R}_1\,\eta^2}\right)\!e^{-\mathcal{R}_1/\eta}
\;+\; B\!\left(\frac{\mathcal{R}_2}{\eta} - \frac{\omega}{\mathcal{R}_2\,\eta^2}\right)\!e^{-\mathcal{R}_2/\eta}
\,\Bigg].
$$

For the CFL bound and the HLLC wave-speed estimates the solver replaces this
with the cheaper Mie–Grüneisen-only over-estimate

$$
c_\text{JWL}^2 \;\approx\; (\omega + 1)\,\frac{p}{\rho},
$$

which bounds the true value from above for the TNT parameter ranges
encountered in practice and is the standard formula used in production blast
codes. Because both the CFL and HLLC use only an upper bound on the wave
speed, this approximation introduces only a marginal time-step penalty
(via a slightly smaller stable $$\Delta t$$) and never compromises stability.

---

## Pressure blending and clamping

The cell pressure is the linear blend (blend) of the air and JWL pressures
weighted by the local burn fraction $$\lambda$$. To absorb tiny round-off
excursions the result is clamped to a positive floor:

$$
p_\text{cell} \;\leftarrow\; \max\!\left(\,
(1 - \lambda)\,p_\text{air}(\rho, e) \;+\; \lambda\,p_\text{JWL}(\rho, e),
\;\; p_\text{floor}
\right),
\quad p_\text{floor} = 10^{-6}\ \text{psi}.
$$

For pure-air cells ($$\lambda = 0$$) the JWL term is skipped entirely and
only $$p_\text{air}$$ is evaluated. The blended sound speed used by HLLC and
the CFL bound is the same convex combination:

$$
c_\text{cell} \;=\; (1 - \lambda)\,c_\text{air}(\rho, p)
              \;+\; \lambda\,c_\text{JWL}(\rho, p).
$$

---

## Inversion: primitives → conservatives

Setting initial conditions or defining inflow states (not used currently)
requires the inverse map: given primitives $$(\rho,\, u_r,\, u_z,\, p)$$ and
$$\lambda$$, compute $$\mathbf{U} = (\rho,\, \rho u_r,\, \rho u_z,\, \rho E)$$.
The total energy density is

$$
\rho E \;=\; \rho\,e \;+\; \tfrac{1}{2}\rho\,(u_r^2 + u_z^2),
$$

so the only nontrivial step is recovering $$e$$ from the pressure–density
pair at the prescribed $$\lambda$$. For pure air this is direct:

$$
e \;=\; \frac{p}{(\gamma - 1)\,\rho}.
$$

For $$\lambda > 0$$ the implicit relation $$p = p_\text{mix}(\rho, e;\lambda)$$
is solved by a Newton iteration starting from the ideal-gas guess. The
Jacobian is constant in $$e$$ at fixed $$\rho$$:

$$
\left.\frac{\partial p_\text{mix}}{\partial e}\right|_\rho
\;=\; (1 - \lambda)\,(\gamma - 1)\,\rho \;+\; \lambda\,\omega\,\rho,
$$

so the iteration

$$
e^{(k+1)}
\;=\; e^{(k)} \;+\; \frac{p_\text{target} - p_\text{mix}\!\left(\rho,\, e^{(k)};\, \lambda\right)}
                          {\left.\partial p_\text{mix}/\partial e\right|_\rho}
$$

converges to relative tolerance $$10^{-10}$$ in two to four iterations. The
implementation is `_invert_mix` in `src/eos.jl`.

---

## Numerical floors

| Floor (code constant) | Value | Role |
|-----------------------|-------|------|
| `RHO_FLOOR` | $$10^{-12}$$ | Lower bound on $$\rho$$ used in primitive recovery and sound speed |
| `E_FLOOR` | $$10^{-6}$$ | Lower bound on specific internal energy entering the EOS |
| `P_FLOOR` | $$10^{-6}$$ | Lower bound on the blended pressure returned by `mix_pressure` |

The floors do not represent physical values but exist solely to preserve
finite arithmetic in unphysical or pathologically-coarse states. Their
chosen values are many orders of magnitude smaller than any flow scale
encountered in a properly-resolved blast simulation, so they are inactive
on the interior mesh during normal operation.
