---
title: Mesh & Discretization
nav_order: 5
description: "Structured Cartesian mesh, Q2 Lagrange quadrilateral elements, and function spaces."
---

# Mesh & Discretization
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Eulerian description

COSMO uses an **Eulerian** (fixed-mesh) description: the computational mesh is static
and fluid material passes through cell faces.  This contrasts with a Lagrangian approach
in which mesh nodes move with the fluid.  For blast problems the large deformations and
material interactions favor the Eulerian frame despite its higher numerical diffusion at
material interfaces.

---

## Mesh construction

The computational domain is the rectangle

$$
\Omega = [r_\text{min},\, r_\text{max}] \times [z_\text{min},\, z_\text{max}]
$$

partitioned into a uniform structured grid of $$n_r \times n_z$$ **quadrilateral cells**.
DOLFINx constructs this mesh via `create_rectangle` with `CellType.quadrilateral`.

Mesh coordinates map directly to physical coordinates:

* coordinate axis 0 → $$r$$ (radial),
* coordinate axis 1 → $$z$$ (axial).

The symmetry axis sits at $$r = r_\text{min} = 0$$ for standard axisymmetric problems.

### Cell size and resolution guidelines

For a charge of radius $$R$$, the cell dimensions are

$$
\Delta r = \frac{r_\text{max} - r_\text{min}}{n_r}, \qquad
\Delta z = \frac{z_\text{max} - z_\text{min}}{n_z}.
$$

A minimum of **~5 cells across the charge radius** is required to capture the initial
blast; 20–30 cells provides good resolution of the peak overpressure.  The examples
provided target 30 cells across the charge radius for production runs.

---

## Finite element spaces

Four function spaces are defined on the mesh.

### Q2 — degree-2 Lagrange quadrilaterals (scalar)

Used for density $$\rho$$, total energy density $$\rho E$$, and pressure output.

A 9-node (full tensor-product) biquadratic basis is used.  On the reference element
$$\hat{\Omega} = [-1,1]^2$$ the 9 basis functions are the tensor products of the
one-dimensional Lagrange polynomials through the nodes $$\{-1, 0, 1\}$$:

$$
\hat{\phi}_{ij}(\xi,\eta) = L_i(\xi)\, L_j(\eta), \quad i,j \in \{0,1,2\},
$$

where $$L_0(\xi) = \tfrac{1}{2}\xi(\xi-1)$$, $$L_1(\xi) = 1 - \xi^2$$,
$$L_2(\xi) = \tfrac{1}{2}\xi(\xi+1)$$.

The 9 local nodes are:

```
(-1,+1) ---- (0,+1) ---- (+1,+1)
   |              |              |
(-1, 0) ---- (0, 0) ---- (+1, 0)
   |              |              |
(-1,-1) ---- (0,-1) ---- (+1,-1)
```

(4 corners, 4 edge midpoints, 1 cell center).

**Global DOF count** (uniform mesh, continuous): approximately
$$(2n_r + 1)(2n_z + 1)$$ for a mesh with $$n_r \times n_z$$ cells.

### Q2v — degree-2 Lagrange quadrilaterals (vector)

Used for the momentum vector $$(\rho u_r, \rho u_z)$$.  This is the vector extension of
Q2: each DOF carries two components, so the global DOF array is twice as long as the
scalar Q2 array.  Components are interleaved (block size 2):
`array = [u_r_0, u_z_0, u_r_1, u_z_1, ...]`.

### DG0 — piecewise-constant discontinuous Galerkin

Used for the material indicator $$\lambda$$ (burn fraction).  Each cell has exactly one
DOF located at the cell centroid; no inter-cell continuity is imposed.  This is
appropriate because $$\lambda$$ jumps discontinuously at the detonation front.

**Global DOF count**: $$n_r \times n_z$$ (one per cell).

### Q1 — degree-1 Lagrange quadrilaterals (scalar)

Used for visualization output only.  The 4-node bilinear basis has one DOF per vertex.
Q2 state fields are interpolated to Q1 before writing to VTK files (see [Output](output)).

---

## MPI domain decomposition

When launched with multiple MPI ranks, DOLFINx automatically partitions the mesh by
distributing contiguous groups of cells to each rank using a graph partitioner (SCOTCH or
ParMETIS via PETSc).  Each rank owns a subset of cells and a layer of **ghost cells**
(cells owned by neighboring ranks that share a face with an owned cell).

Ghost DOF values are synchronized via `ghostUpdate()` calls after each vector assembly,
ensuring that the explicit update $$\mathbf{U} \leftarrow \mathbf{U} + \Delta t\,\mathbf{M}^{-1}\mathbf{R}$$
is consistent across ranks.

The CFL time step computation and the minimum cell size are reduced across all ranks
using `MPI.allreduce` with the `MPI.MIN` / `MPI.MAX` operators.

---

## Quadrature

All integrals in the weak form are evaluated with **Gauss–Legendre quadrature** on the
reference quadrilateral.  The mass-matrix form involves Q2 × Q2 × $$r$$ (polynomial
degree $$\leq 5$$ in each variable on a uniform mesh), which requires at least a
$$3 \times 3$$ Gauss rule to integrate exactly.  DOLFINx is instructed to use
`quadrature_degree = 5` for the lumped-mass assembly (see [Time Integration](time-integration)).
