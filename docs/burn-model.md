---
title: Programmed Burn Model
nav_order: 4
description: "Detonation kinematics, burn fraction, charge geometry, initial conditions, and the explosive-material passive scalar."
---

# Programmed Burn Model
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Detonation physics background

A high explosive detonates when a supersonic shock — the **detonation
front** — propagates through the unreacted material and converts it
essentially instantaneously into hot, high-pressure product gases. The
steady propagation speed $$D$$ and the thermodynamic state immediately
behind the front (the **Chapman–Jouguet** or CJ state) are determined by
the Rankine–Hugoniot conditions applied to the detonation wave.

For TNT, the reference values built into `src/eos.jl` are

$$
D_\text{CJ} \;\approx\; 6930\ \text{m/s} \;\approx\; 2.728 \times 10^5\ \text{in/s},
\qquad
p_\text{CJ} \;\approx\; 21\ \text{GPa} \;\approx\; 3.05 \times 10^6\ \text{psi},
$$

$$
e_\text{CJ} \;\approx\; 6.0\ \text{MJ/kg} \;\approx\; 9.30 \times 10^9\ \text{in}^2/\text{s}^2.
$$

A rigorous reactive-flow treatment of detonation requires one or more
progress variables, an Arrhenius-type rate law, and very fine resolution of
the reaction zone. The **programmed-burn** approximation replaces this
stiff system with a purely *geometric* rule: the detonation front is
prescribed to travel at velocity $$D$$ from the initiation point, and the
high-pressure CJ-state energy is pre-loaded into every charge cell at
$$t = 0$$. This eliminates the reaction kinetics, lets the explicit time
step reflect the flow CFL alone, and produces excellent far-field blast
predictions for engineering purposes.

---

## Initial-condition modes

Two IC modes are supported, selected by the JSON key
`charge.initial_condition`. Both are implemented in
[`src/charge.jl`](https://github.com/jman87/cosmo/blob/main/src/charge.jl).

| Mode | Key | In-charge state at $$t = 0$$ | EOS evolution |
|------|-----|------------------------------|---------------|
| Programmed burn (default) | `"programmed_burn"` | $$\rho = \rho_0,\;\; e = e_\text{CJ},\;\; \mathbf{u} = \mathbf{0}$$, $$\lambda(r, z, 0) = 0$$ | $$\lambda$$ is updated kinematically each step until the front sweeps past every charge cell |
| Brode compressed-gas balloon | `"brode"` (alias `"balloon"`) | $$\rho = \rho_0,\;\; e = Q_\text{TNT},\;\; \mathbf{u} = \mathbf{0}$$, $$\lambda = 1$$ | $$\lambda$$ is frozen at $$t = 0$$; full JWL behaviour from the first step |

The Brode IC uses the **heat of detonation**

$$
Q_\text{TNT} \;\approx\; 4.184 \times 10^6\ \text{J/kg}
            \;\approx\; 6.486 \times 10^9\ \text{in}^2/\text{s}^2
$$

(coded as `4.184e6 * J_KG_TO_IN2_S2` in `build_charge`) instead of
$$e_\text{CJ}$$. The Brode model is unconditionally robust on any mesh
because there is no detonation front to resolve, but its initial pressure
field is unphysical (a uniform high-pressure sphere) and its near-field
shock arrival times are off by $$\mathcal{O}(R / D)$$. It is useful for
very coarse meshes and for verification cases where the only requirement
is mass and energy conservation.

The **programmed-burn** mode is the production setting and matches the
behaviour of the original Python sibling solver.

---

## Burn fraction

At simulation time $$t$$, a point at position $$(r, z)$$ inside the charge
is considered to have been reached by the detonation front if the
straight-line distance from that point to the charge centre satisfies

$$
d(r, z) \;\leq\; D\,t,
\qquad
d(r, z) \;=\; \sqrt{(r - r_c)^2 + (z - z_c)^2}.
$$

The discrete cell-centred burn fraction is therefore

$$
\lambda_{ij}(t) \;=\;
\begin{cases}
  1, & \text{if cell } (i, j) \text{ is inside the charge and } D\,t \;\geq\; d(r_i, z_j), \\
  \lambda_{ij}(t - \Delta t), & \text{otherwise.}
\end{cases}
\tag{burn rule}
$$

Two consequences are worth noting:

1. The rule is **monotone**: once $$\lambda_{ij} = 1$$ the cell never
   reverts to unburned, so `update_burn!` only needs to walk over cells
   with $$\lambda_{ij} = 0$$ each step. Cells outside the charge never
   reach $$\lambda = 1$$.
2. The transition is **piecewise binary**: an individual cell's
   $$\lambda$$ jumps from 0 to 1 within a single step. The coarse-grained
   smoothness of the burn surface comes from neighbouring cells flipping
   on adjacent steps as the front sweeps across them, not from a smooth
   sub-cell ramp.

The burn rule treats the detonation as a spherical wave from the charge
centre regardless of the charge shape — a simplification that is exact
for a point-initiated sphere and a reasonable approximation for a
centre-initiated cylinder. End-initiated geometries with a planar wave
front would require user-supplied burn-tables or a level-set initialisation
of the front.

The burn fraction is **constant within an SSP-RK substep**: it is updated
at the start of each outer step and then frozen for both stage residuals.
Holding $$\lambda$$ fixed within the RK is consistent because the
programmed-burn front is purely kinematic — it does not depend on the
flow field — and the resulting two-stage update is genuinely 2nd-order
accurate everywhere except in the handful of cells the front passes
through during that step.

---

## Charge geometry

The solver computes the charge **radius** from a user-supplied charge
size and the unreacted-explosive density $$\rho_0$$ from the TNT JWL
parameters. The size is given in the input file as **either** a mass

$$
m \;\bigl[\text{lbf}\cdot\text{s}^2/\text{in}\bigr] \;\;\text{(\texttt{charge.mass})},
$$

**or** a weight in pound-force

$$
W\;[\text{lbf}] \;\;\text{(\texttt{charge.weight})},
$$

but not both — the loader rejects files with both keys present and
files using the legacy `weight_lbm` key. When a weight is supplied, the
solver converts to mass via the consistent-unit form of Newton's second
law

$$
m \;=\; \frac{W}{g_c}, \qquad g_c = 386.088\ \text{in}/\text{s}^2,
$$

so that mass $$\times$$ acceleration ($$\text{lbf}\cdot\text{s}^2/\text{in}$$
$$\times$$ $$\text{in}/\text{s}^2$$) is directly in lbf. Internally the
solver only ever uses the mass form; the pound-mass (lbm) unit is not
used or accepted anywhere in the codebase. The volume is then

$$
V \;=\; \frac{m}{\rho_0}.
$$

The radius is then extracted by shape:

* **Sphere** (`charge.shape = "sphere"`):

$$
R \;=\; \left(\frac{3\,V}{4\pi}\right)^{1/3}.
$$

* **Right-circular cylinder** (`charge.shape = "cylinder"`, half-height
  $$H$$ supplied as `charge.cylinder_half_height`):

$$
R \;=\; \sqrt{\frac{V}{2\pi\,H}}.
$$

The charge occupies the region

$$
\text{sphere:} \qquad (r - r_c)^2 + (z - z_c)^2 \;\leq\; R^2,
$$

$$
\text{cylinder:} \qquad |z - z_c| \;\leq\; H \;\;\text{and}\;\; (r - r_c)^2 \;\leq\; R^2.
$$

The cylinder is right-circular in 3-D (no eccentricity); its
azimuthal extent is automatic from the axisymmetric formulation. The
implementation is `inside_charge` in `src/charge.jl`.

---

## Initial conservative state

At $$t = 0$$ the fluid is at rest ($$\mathbf{u} = \mathbf{0}$$) and the
density and energy are blended on a per-cell basis using a sub-cell
volume-fraction sampler. For each cell the solver evaluates a
$$5 \times 5$$ array of sub-points, computes the fraction of those points
that lie inside the charge geometry, and uses that fraction $$\phi$$ to
linearly blend the in-charge and ambient states:

$$
\rho_{ij} \;=\; \phi_{ij}\,\rho_0 \;+\; (1 - \phi_{ij})\,\rho_\text{air},
$$

$$
(\rho E)_{ij} \;=\; \phi_{ij}\,\rho_0\,e_\text{in} \;+\; (1 - \phi_{ij})\,\rho_\text{air}\,e_\text{air},
$$

with $$e_\text{in} = e_\text{CJ}$$ for `programmed_burn` and
$$e_\text{in} = Q_\text{TNT}$$ for `brode`. The momentum components are
zero. The $$5 \times 5$$ sub-sampling resolves cells the geometric
boundary cuts through to roughly 4 % accuracy, which is sufficient given
the dominant near-field error sources (mesh resolution and the
binary-burn rule). The implementation is `set_initial_state!` in
`src/charge.jl`.

The initial reaction-progress field is then

$$
\lambda_{ij}(0) \;=\;
\begin{cases}
  \phi_{ij}, & \text{Brode mode (frozen for the run)}, \\
  0, & \text{programmed-burn mode (updated each step)}.
\end{cases}
$$

---

## Material tag (passive scalar)

In addition to the four-component conservative state, the solver carries
a passive scalar to mark the **explosive-origin mass fraction**
$$Y \in [0, 1]$$:

* $$Y = 1$$ inside the original charge geometry;
* $$Y = 0$$ in the ambient air at $$t = 0$$;
* $$Y$$ is conserved under flow (no source term).

The conservative form $$\rho Y$$ is advected by

$$
\frac{\partial(\rho Y)}{\partial t}
  \;+\; \frac{1}{r}\frac{\partial(r\,\rho u_r Y)}{\partial r}
  \;+\; \frac{\partial(\rho u_z Y)}{\partial z}
  \;=\; 0.
$$

The face-fluxes for $$\rho Y$$ reuse the **HLLC mass flux** $$F^{(\rho)}$$
already computed for the Euler update. Donor-cell upwinding on $$Y$$ —
i.e.

$$
(\rho u Y)_\text{face}
\;=\;
\begin{cases}
  F^{(\rho)}_\text{face}\,Y_L, & F^{(\rho)}_\text{face} \geq 0, \\[2pt]
  F^{(\rho)}_\text{face}\,Y_R, & F^{(\rho)}_\text{face} < 0,
\end{cases}
$$

— is the simplest 1st-order monotonic discretisation. It is sufficient
for what $$Y$$ is used for: a visualisation tag that follows the
products outward as they expand into the air. A higher-order
reconstruction on $$Y$$ would not change any physical quantity in the
solution.

The material tag is purely diagnostic — it has **no effect on the EOS or
on the burn fraction**. The reaction-progress variable $$\lambda$$ controls
the EOS blend; once a cell has finished burning ($$\lambda = 1$$) the
products' subsequent expansion into ambient air is captured by the JWL
expansion tail and the velocity field, while $$Y$$ records the material
provenance for post-processing.

The post-processed `material_id` field written to VTK output is
$$Y = \rho Y / \max(\rho,\, \rho_\text{floor})$$, clipped to
$$[0, 1]$$ to absorb tiny round-off excursions from the donor-cell update.

---

## Limitations and possible extensions

| Limitation | Effect | Possible remedy |
|------------|--------|-----------------|
| Instantaneous burn at the front | No reaction-zone width; cell pressure jumps from air-state to JWL within one cell | History-variable programmed burn or Lee–Tarver ignition-and-growth |
| Burn front spherical from the centre | Cylindrical charges initiated end-on have a different wave shape | User-prescribed burn-table or level-set initialisation of the front |
| Pre-loaded CJ energy | Chemical energy released at $$t = 0$$, not progressively as the front sweeps | Energy-release source $$\dot{q} = \rho\,Q\,\dot{\lambda}$$ tied to a non-binary $$\lambda$$ |
| No afterburn / post-detonation oxidation | Late-time reaction with ambient oxygen ignored | Two-phase or afterburn source terms |
| Single material per cell | Pressure linearly blended via $$\lambda$$; no sharp interface | Level-set or VOF interface tracking |
| TNT-only built-in JWL | Other explosives need user-supplied parameters | Parameter table read from the JSON input |
