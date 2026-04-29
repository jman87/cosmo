---
title: Home
nav_order: 1
description: "Theory and numerical methods for the COSMO 2-D axisymmetric finite-volume air-blast solver."
permalink: /
---

# COSMO Theory Manual
{: .no_toc }

**COSMO** (**C**omputational **O**verpressure **S**hock **M**odelling **O**perator)
is a 2-D axisymmetric, MPI-parallel **finite-volume** solver for near-field
air-blast problems written in Julia. It integrates the inviscid compressible
Euler equations on a uniform structured Cartesian $$(r, z)$$ mesh using the
**HLLC** approximate Riemann solver, **MUSCL** spatial reconstruction with the
**minmod** slope limiter, and **strong-stability-preserving** (SSP) explicit
Runge–Kutta time stepping.

---

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Physics and capabilities

| Feature | Description |
|---------|-------------|
| **Governing equations** | Inviscid compressible Euler equations in 2-D axisymmetric form |
| **Domain** | Structured Cartesian $$(r, z)$$ quarter-plane; symmetry axis at $$r = 0$$ |
| **Charge shapes** | Sphere or right-circular cylinder, sized from user-supplied weight |
| **Explosive model** | Programmed-burn detonation front from the initiation point at velocity $$D$$, or Brode compressed-gas balloon |
| **Explosive EOS** | Jones–Wilkins–Lee (JWL) for TNT detonation products |
| **Ambient EOS** | Calorically perfect ideal gas ($$\gamma = 1.4$$) for air |
| **Spatial discretisation** | Cell-centred finite volume on a structured $$(r, z)$$ mesh with annular cell volumes $$V = r_c\,\Delta r\,\Delta z$$ |
| **Reconstruction** | MUSCL with the minmod slope limiter (TVD; 2nd-order in smooth flow, 1st-order at extrema) |
| **Riemann solver** | HLLC (Toro, Spruce & Speares 1994), evaluated on r- and z-faces with EOS-blended sound speeds |
| **Time integration** | Forward Euler (debug) or SSP-RK2 / Heun (default), CFL-controlled |
| **Material tracking** | Reaction-progress variable $$\lambda \in [0,1]$$ blends air and JWL pressures; passive scalar $$\rho Y$$ tags explosive-origin material |
| **Parallelism** | MPI Cartesian topology with geometric block decomposition and 2-cell ghost halo exchange |
| **Output** | Per-step `.vtr` rectilinear-grid snapshots gathered onto rank 0 + master `.pvd` collection for ParaView |
| **Units** | Imperial gravitational (lbf, s, in) throughout |

---

## Solution process

```
Load JSON input ─ config.jl
        │
        ▼
Build MPI Cartesian topology and per-rank structured grid with
ng = 2 ghost layers ─ parallel.jl, grid.jl
        │
        ▼
Construct charge geometry and assign t = 0 conservative state
on each cell using a 5x5 sub-cell volume-fraction sampler
─ charge.jl
        │
        ▼
 ╔══════════════════════════════════════════════════════════════════╗
 ║  Explicit time loop (solver.jl, timeint.jl)                       ║
 ║                                                                   ║
 ║  while t < t_end:                                                 ║
 ║    1. Update programmed-burn fraction λ(r, z, t)                  ║
 ║       and exchange its halos                                       ║
 ║    2. Compute MPI-reduced CFL time step Δt                        ║
 ║    3. Repeat for each SSP-RK stage:                               ║
 ║       a. exchange_halos! on conservative state U and ρY           ║
 ║       b. apply_bcs! on physical edges (mirror or zero-gradient)   ║
 ║       c. MUSCL-reconstruct UL/UR at every interior r- and z-face  ║
 ║       d. Evaluate HLLC numerical flux F, G with EOS-blended c     ║
 ║       e. Assemble residual L(U) from the discrete divergence      ║
 ║          theorem with annular face areas; add geometric source    ║
 ║          (0, p/r_c, 0, 0) for r-momentum                          ║
 ║       f. Donor-cell upwind update of passive scalar ρY            ║
 ║       g. Convex SSP combination → next stage                      ║
 ║    4. Write .vtr snapshot every output_frequency steps            ║
 ╚══════════════════════════════════════════════════════════════════╝
        │
        ▼
Write final snapshot; finalise the .pvd time series ─ output.jl
```

---

## Document structure

| Section | Topic |
|---------|-------|
| [Governing Equations](governing-equations) | Axisymmetric Euler equations, conservative state, primitives, sound speed |
| [Equations of State](eos) | JWL (TNT) and ideal-gas models; reaction-progress pressure blend; Newton inversion |
| [Programmed Burn Model](burn-model) | Charge geometry, detonation front, burn fraction, Brode IC, material tag |
| [Mesh and Grid](mesh-discretization) | Structured $$(r, z)$$ grid, ghost cells, indexing convention |
| [Finite-Volume Discretisation](finite-volume) | Cell averages, discrete divergence theorem in cylindrical coordinates, residual assembly |
| [Riemann Solver and Reconstruction](riemann-reconstruction) | HLLC numerical flux + MUSCL minmod reconstruction |
| [Time Integration](time-integration) | SSP-RK2 (Heun) and forward Euler; CFL bound; detonation-velocity safeguard |
| [Boundary Conditions](boundary-conditions) | Ghost-cell mirror reflection, zero-gradient outflow, halo seams |
| [MPI Parallelisation](parallelisation) | Cartesian topology, block partitioning, halo exchange, gather-based output |
| [Output](output) | Cell-centred fields, MPI gather, `.vtr` + `.pvd` for ParaView |
| [References](references) | Bibliography |

---

## Unit system

All physical quantities use the **imperial gravitational (consistent)
unit system**: forces are in pound-force (lbf) and masses are in
$$\text{lbf}\cdot\text{s}^2 / \text{in}$$. The pound-mass (lbm) is **never
used** anywhere in the solver — neither in code, nor in input files, nor
in printed output — to eliminate the recurring lbm/lbf ambiguity that
plagues mixed-unit weapons-physics literature.

| Quantity | Unit | Symbol |
|----------|------|--------|
| Length | inch | $$\text{in}$$ |
| Time | second | $$\text{s}$$ |
| Force / weight | pound-force | $$\text{lbf}$$ |
| Mass | pound second squared per inch | $$\text{lbf} \cdot \text{s}^2 / \text{in}$$ |
| Pressure | pound per square inch | $$\text{psi} = \text{lbf}/\text{in}^2$$ |
| Density | pound second squared per inch^4 | $$\text{lbf} \cdot \text{s}^2 / \text{in}^4$$ |
| Specific energy | square inch per second squared | $$\text{in}^2/\text{s}^2$$ |

In this **consistent system** Newton's second law is simply $$F = m a$$
(no $$g_c$$ correction factor) with $$F$$ in lbf, $$m$$ in
$$\text{lbf}\cdot\text{s}^2/\text{in}$$, and $$a$$ in
$$\text{in}/\text{s}^2$$. The standard acceleration of gravity
$$g_c = 386.088\ \text{in}/\text{s}^2$$ enters only as the
**weight↔mass conversion** at the input boundary:

$$
m\;\bigl[\text{lbf}\cdot\text{s}^2/\text{in}\bigr]
\;=\;\frac{W\;[\text{lbf}]}{g_c}.
$$

The user may supply the explosive charge as either a **mass**
(`charge.mass`, in $$\text{lbf}\cdot\text{s}^2/\text{in}$$) or a
**weight** (`charge.weight`, in lbf) — but not both. The solver stores
the value internally as a mass and never reintroduces lbm.

Conversion of the JWL parameters from SI uses

$$
1\ \text{Pa} = 1.450\,377 \times 10^{-4}\ \text{psi}, \quad
1\ \text{kg/m}^3 = 9.357 \times 10^{-8}\ \text{lbf}\cdot\text{s}^2/\text{in}^4, \quad
1\ \text{J/kg} = 1550.003\ \text{in}^2/\text{s}^2.
$$

### Input-file format

Input files are parsed by `read_json_with_comments`, which strips
`//` line comments and `/* ... */` block comments before invoking
`JSON.parse`. Files may be named `.json` (strict) or `.jsonc`
(JSON-with-comments) — both extensions are accepted.

---

## Software dependencies

| Package | Role |
|---------|------|
| [Julia ≥ 1.10](https://julialang.org) | Implementation language |
| [`MPI.jl`](https://juliaparallel.org/MPI.jl/stable/) | Cartesian topology, halo exchange, allreduce, gatherv |
| [`WriteVTK.jl`](https://github.com/JuliaIO/WriteVTK.jl) | Rectilinear `.vtr` snapshots and `.pvd` collection |
| [`JSON.jl`](https://github.com/JuliaIO/JSON.jl) | Input-file parsing |

The MPI runtime is supplied by the `MPI.jl` JLL build by default; the user may
optionally bind to a system MPI (e.g. OpenMPI / MPICH) via
`MPIPreferences.use_system_binary()`.
