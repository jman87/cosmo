---
title: Home
nav_order: 1
description: "Theory and numerical methods for the COSMO 2-D axisymmetric air-blast solver."
permalink: /
---

# COSMO Theory Manual
{: .no_toc }

**COSMO** is a 2-D axisymmetric, MPI-parallel finite-volume solver for near-field
air-blast problems. The solver integrates the compressible Euler equations on a structured
Cartesian $$(r, z)$$ mesh using a continuous Galerkin formulation with explicit time stepping.

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
| **Charge shapes** | Sphere or right-circular cylinder |
| **Explosive model** | Programmed-burn detonation front at user-specified velocity |
| **Explosive EOS** | Jones–Wilkins–Lee (JWL) for TNT detonation products |
| **Ambient EOS** | Calorically perfect ideal gas ($$\gamma = 1.4$$) for air |
| **Elements** | Degree-2 Lagrange quadrilaterals (Q2, 9-node) on all state variables |
| **Time integration** | Forward Euler · SSP-RK2 · SSP-RK3 (CFL-controlled) |
| **Shock capturing** | Von Neumann–Richtmyer artificial bulk viscosity |
| **Parallelism** | MPI domain decomposition via DOLFINx / PETSc |
| **Output** | Per-step `.vtu` snapshots + master `.pvd` for ParaView |
| **Units** | Imperial gravitational (lbf, s, in) throughout |

---

## Solution process

The overall solution algorithm proceeds as follows.

```
Load & validate JSON input (config.py)
        │
        ▼
Build structured quad mesh; allocate Q2, Q2v, DG0, Q1 function spaces (mesh.py)
        │
        ▼
Instantiate JWL (TNT) and ideal-gas EOS objects (eos.py)
        │
        ▼
Compute charge geometry; set initial conditions (burn.py)
   ┌─ Inside charge: ρ = ρ₀_TNT,  e = e_CJ,  u = 0
   └─ Outside charge: ρ = ρ_air,  e = e_air, u = 0
        │
        ▼
 ╔══════════════════════════════════════════════════════════╗
 ║  Explicit time loop  (time_integration.py)               ║
 ║                                                          ║
 ║  while t < t_end:                                        ║
 ║    1. Update burn fraction λ(r,z,t) → mat field          ║
 ║    2. Compute CFL time step Δt                           ║
 ║    3. SSP-RK stage(s):                                   ║
 ║       a. Assemble residual R(U) via UFL (weak_form.py)   ║
 ║       b. Update:  U ← U + Δt M⁻¹ R(U)                   ║
 ║       c. Enforce symmetry BCs on momentum                ║
 ║    4. Write .vtu snapshot every output_frequency steps   ║
 ╚══════════════════════════════════════════════════════════╝
        │
        ▼
Write final snapshot; close PVD file (output.py)
```

---

## Document structure

| Section | Topic |
|---------|-------|
| [Governing Equations](governing-equations) | Axisymmetric Euler equations, conserved variables |
| [Equations of State](eos) | JWL (TNT) and ideal-gas models; pressure blending |
| [Programmed Burn Model](burn-model) | Detonation front, burn fraction, initial conditions |
| [Mesh & Discretization](mesh-discretization) | Structured mesh, Q2 elements, function spaces |
| [Weak Formulation](weak-form) | Variational form in cylindrical coordinates, IBP |
| [Shock Capturing](shock-capturing) | Von Neumann–Richtmyer artificial bulk viscosity |
| [Time Integration](time-integration) | SSP-RK schemes, lumped mass matrix, CFL condition |
| [Boundary Conditions](boundary-conditions) | Symmetry axis, ground plane, outer walls |
| [Output](output) | Field interpolation, VTK format, ParaView |
| [References](references) | Bibliography |

---

## Unit system

All physical quantities use the **imperial gravitational unit system**:

| Quantity | Unit | Symbol |
|----------|------|--------|
| Length | inch | in |
| Time | second | s |
| Force | pound-force | lbf |
| Mass | lbf·s²/in | — |
| Pressure | pound per square inch | psi = lbf/in² |
| Density | lbf·s²/in⁴ | — |
| Specific energy | square inch per second squared | in²/s² |

Newton's second law reads $$F = ma$$ with $$F$$ in lbf, $$m$$ in lbf·s²/in, and $$a$$ in in/s².
A body with weight $$W$$ lbm has mass $$m = W / g_c$$ where $$g_c = 386.088$$ in/s² is the
standard acceleration of gravity.

---

## Dependencies

| Package | Role |
|---------|------|
| [DOLFINx (FEniCSx)](https://fenicsproject.org) | Mesh, function spaces, UFL form assembly |
| [mpi4py](https://mpi4py.readthedocs.io) | MPI parallelism |
| [PETSc](https://petsc.org) | Vector/matrix assembly (via DOLFINx) |
| [NumPy](https://numpy.org) | Array operations in Python |
| [jsonschema](https://python-jsonschema.readthedocs.io) | Input file validation |
