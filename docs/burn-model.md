---
title: Programmed Burn Model
nav_order: 4
description: "Detonation model, burn fraction, and initial condition assignment."
---

# Programmed Burn Model
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Detonation physics background

A high explosive detonates when a supersonic shock wave — the **detonation front** —
propagates through the unreacted material and instantaneously converts it to hot, high-
pressure product gases.  The steady propagation speed $$D$$ and the thermodynamic state
immediately behind the front (the **Chapman–Jouguet (CJ) state**) are determined by the
Rankine–Hugoniot conditions applied to the detonation wave.

For TNT the CJ detonation velocity is $$D \approx 6.86\text{–}8.31\ \text{km/s}$$
$$(2.70\text{–}3.27)\times10^5\ \text{in/s}$$, and the CJ pressure is
$$p_\text{CJ} \approx 21\ \text{GPa}$$ $$(3.05 \times 10^6\ \text{psi})$$.

Rigorous treatment of detonation requires a reactive flow model with one or more progress
variables and an Arrhenius-type reaction rate.  The **programmed burn** approximation
replaces this reactive system with a purely geometric rule: the detonation front is
prescribed to travel at velocity $$D$$ from the initiation point, eliminating the need to
solve stiff reaction-rate equations.

---

## Programmed burn algorithm

### Burn fraction

At simulation time $$t$$, a material point at position $$(r, z)$$ has been reached by the
detonation front if the straight-line distance $$d$$ from that point to the charge center
$$(r_c, z_c)$$ satisfies

$$
d \leq D \cdot t, \qquad d = \sqrt{(r - r_c)^2 + (z - z_c)^2}.
$$

The **burn fraction** is

$$
\lambda(r, z, t) =
\begin{cases}
  1 & \text{if inside charge and } D\,t \geq d, \\
  0 & \text{otherwise.}
\end{cases}
\tag{burn fraction}
$$

$$\lambda = 1$$ means fully converted to detonation products (JWL EOS active);
$$\lambda = 0$$ means air or unreacted explosive (ideal-gas EOS active).

The burn fraction is stored in a piecewise-constant DG0 function (`mat`) and is updated
once per time step — before the RK stage evaluations — from the analytic formula above.
Holding `mat` fixed within an RK sub-step is consistent because the programmed-burn front
is not influenced by the computed flow field.

### Charge geometry

The solver computes the charge radius from the user-supplied weight (in lbm) and the TNT
reference density $$\rho_0$$:

1. Convert weight to mass:
$$
m = \frac{W_\text{lbm}}{g_c}, \qquad g_c = 386.088\ \text{in/s}^2.
$$

2. Compute volume:
$$
V = \frac{m}{\rho_0}.
$$

3. Extract radius by shape:

   * **Sphere:** $$R = \left(\dfrac{3V}{4\pi}\right)^{1/3}$$
   * **Right-circular cylinder** (half-height $$H$$ user-supplied):
     $$R = \sqrt{\dfrac{V}{2\pi H}}$$

The charge occupies the region

$$
\text{sphere:} \quad (r-r_c)^2 + (z-z_c)^2 \leq R^2,
$$

$$
\text{cylinder:} \quad |z - z_c| \leq H \;\text{ and }\; (r - r_c)^2 \leq R^2.
$$

Note that the axisymmetric cylinder is a right-circular cylinder in 3-D (no eccentricity).

---

## Initial conditions

At $$t = 0$$ all material is at rest ($$\mathbf{u} = \mathbf{0}$$).  Density and energy are
assigned by material region:

### Inside the charge

The explosive is initialised at the CJ state with the full chemical energy pre-loaded:

$$
\rho = \rho_0 = 1.526 \times 10^{-4}\ \text{lbf·s}^2/\text{in}^4,
$$

$$
e = e_\text{CJ} = 9.30 \times 10^9\ \text{in}^2/\text{s}^2,
$$

$$
\rho E = \rho\, e \quad (\text{no kinetic energy at rest}).
$$

This pre-loading of $$e_\text{CJ}$$ is the key simplification of the programmed-burn
approach: the chemical energy is treated as already released, and the detonation model
merely controls *when* each point transitions from the unreacted state to the high-energy
product state.  A more sophisticated treatment would use a reaction progress variable and
release the energy incrementally as the burn fraction increases.

### Outside the charge

Ambient air at standard sea-level conditions:

$$
\rho = \rho_\text{air} = 1.146 \times 10^{-7}\ \text{lbf·s}^2/\text{in}^4,
$$

$$
e = e_\text{air} = \frac{p_\text{atm}}{(\gamma - 1)\,\rho_\text{air}}
  = 3.204 \times 10^8\ \text{in}^2/\text{s}^2.
$$

### Material indicator

The DG0 material field is set to $$\lambda = 1$$ inside the charge and $$\lambda = 0$$
outside.  Because DG0 DOFs are cell-centered, the indicator uses the cell-centroid
coordinates to determine which region each cell falls in.

---

## Limitations and potential improvements

| Limitation | Effect | Possible remedy |
|------------|--------|-----------------|
| Instantaneous burn at front | No partial-reaction zone; detonation products at full CJ state | History-variable programmed burn; Lee–Tarver ignition & growth |
| Burn direction is spherical from center | Cylindrical charges ignited end-on have a different wave shape | User-prescribed detonation vector or fit-to-test front shape |
| No afterburn | Air–product mixing and late-time oxidation not modelled | Two-phase or afterburn source terms |
| Single material per cell | No sharp interface; pressure blended by $$\lambda$$ | Level-set or VOF interface tracking |
| Pre-loaded CJ energy | Instantaneous energy release ignores the wave structure | Reactive-flow EOS with progress variable |
