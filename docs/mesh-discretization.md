---
title: Mesh and Grid
nav_order: 5
description: "Structured Cartesian (r, z) finite-volume mesh, ghost cells, and indexing convention."
---

# Mesh and Grid
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Eulerian, cell-centred description

COSMO uses a fixed (Eulerian) **cell-centred finite-volume** description.
The computational mesh is a uniform structured grid; fluid material passes
through cell faces and the conserved quantities $$(\rho,\, \rho u_r,\,
\rho u_z,\, \rho E)$$ are stored as cell-volume averages at the cell
centroids. This contrasts with a Lagrangian description in which mesh
nodes move with the fluid, and with a finite-element description in which
the unknowns are nodal coefficients of a continuous polynomial expansion.
For blast problems the very large deformations and material interactions
favour the Eulerian frame even though it incurs more numerical diffusion at
material interfaces than a Lagrangian or interface-tracking method.

---

## Global computational domain

The global mesh covers the rectangle

$$
\Omega \;=\; [r_\text{min},\, r_\text{max}] \;\times\; [z_\text{min},\, z_\text{max}]
$$

partitioned into $$n_r \times n_z$$ uniform rectangular cells:

$$
\Delta r \;=\; \frac{r_\text{max} - r_\text{min}}{n_r},
\qquad
\Delta z \;=\; \frac{z_\text{max} - z_\text{min}}{n_z}.
$$

By convention, axis 0 of the mesh is the **radial** direction $$r$$ and
axis 1 is the **axial** direction $$z$$. The symmetry axis sits at
$$r = r_\text{min} = 0$$ for standard axisymmetric problems.

### Cell centres and faces

Cell $$(i, j)$$ has

* **centre** at $$(r_i,\, z_j)$$ with
  $$r_i = r_\text{min} + (i - \tfrac{1}{2})\,\Delta r$$ and
  $$z_j = z_\text{min} + (j - \tfrac{1}{2})\,\Delta z$$,
* **radial faces** at $$r = r_{i-1/2}$$ (left) and $$r = r_{i+1/2}$$
  (right) with $$r_{i \pm 1/2} = r_\text{min} + (i \mp \tfrac{1}{2})\,\Delta r$$
  modulo the offset by one,
* **axial faces** at $$z = z_{j-1/2}$$ (bottom) and $$z = z_{j+1/2}$$ (top).

Because all azimuthal integrals contribute a uniform factor $$2\pi$$, the
solver works with the **annular** form of the cell volume and face areas
(the $$2\pi$$ cancels everywhere):

$$
V_{ij} \;=\; r_i\,\Delta r\,\Delta z,
$$

$$
A^{(r)}_{i\pm 1/2,\,j} \;=\; r_{i\pm 1/2}\,\Delta z,
\qquad
A^{(z)}_{i,\,j\pm 1/2} \;=\; r_i\,\Delta r.
$$

These appear directly in the discrete divergence theorem
([Finite-Volume Discretisation](finite-volume)).

### Resolution guidelines

Empirical practice for free-air blast problems is

* a **minimum of 5 cells across the charge radius** to capture the
  initial expansion;
* **20–30 cells across the charge radius** for converged peak overpressure
  in the near field;
* outer-domain mesh size matched to the inner mesh — the structured
  uniform grid does not stretch.

The bundled `examples/sphere/input.jsonc` uses 96 × 96 cells on a
$$60 \times 60$$ in domain, giving $$\Delta r = \Delta z = 0.625$$ in and
roughly 2.5 cells across the 1.59 in TNT-sphere radius. This is a
deliberately coarse smoke-test mesh; the `coarse_sphere` example is even
coarser. Production runs typically use $$256 \times 256$$ or larger.

---

## Ghost cells

Each rank stores its block of interior cells **plus** $$n_g$$ layers of
**ghost cells** on every side. The ghost layers serve two roles:

1. They provide the four-cell stencil needed for MUSCL reconstruction
   (`U_LL, U_L, U_R, U_RR`) at the first and last interior face;
2. They store either physical-boundary BC values (filled by `apply_bcs!`
   on edge ranks) or off-rank neighbour values (filled by
   `exchange_halos!` on internal seams).

The solver uses

$$
n_g \;=\; \texttt{NGHOST} \;=\; 2,
$$

which is the minimum required for the 4-point MUSCL stencil with minmod.
Increasing $$n_g$$ would only be necessary for a wider stencil
(e.g. WENO-5 or higher).

### Indexing convention

Per-rank arrays carry the ghost layers explicitly:

* total radial width: $$n_r^\text{tot} \;=\; n_r + 2 n_g$$
  (`grid.nri` in code),
* total axial width: $$n_z^\text{tot} \;=\; n_z + 2 n_g$$
  (`grid.nzj` in code),
* the **interior** cell index range is
  $$i \in [n_g + 1,\; n_g + n_r]$$,
  $$j \in [n_g + 1,\; n_g + n_z]$$.

In source code these ranges appear as `interior_i(grid)` /
`interior_j(grid)` (defined in
[`src/grid.jl`](https://github.com/jman87/cosmo/blob/main/src/grid.jl)).
The face index `iface` lies between cells `iface - 1` and `iface`, so the
face arrays in `src/timeint.jl` carry one more entry than the cell arrays
in each direction.

The cell-centre and face coordinates are stored at their **physical**
locations, so a ghost cell that sits adjacent to the symmetry axis has a
mathematically negative $$r$$. This negative-$$r$$ value is not used as a
geometry input — the volume and face areas in the residual always use
interior cell-centres — but it preserves the natural indexing for the
mirror-reflection BC at $$r = 0$$.

---

## Per-rank block partitioning

When launched on more than one MPI rank the global $$n_r \times n_z$$
mesh is split across a $$P_r \times P_z$$ Cartesian process grid whose
factorisation is chosen by `MPI.Dims_create`. Each rank owns a contiguous
rectangular block of interior cells; the `block_extents(n_global, np, ip)`
helper returns the block offset and local cell count for rank coordinate
`ip` along one direction. Remainder cells are spread over the
lowest-coordinate ranks so the load imbalance is at most one cell per rank
in either direction.

Single-rank runs reduce to a $$1 \times 1$$ process grid with
`MPI.PROC_NULL` neighbours, so the same code path serves serial and
parallel runs without branching.

The full discussion of the parallel topology, halo exchange, and
collective operations is in [MPI Parallelisation](parallelisation).

---

## State and flux arrays

For each rank the time stepper allocates the following arrays
(types `State` and `Grid` in `src/timeint.jl` and `src/grid.jl`):

| Array | Shape | Role |
|-------|-------|------|
| `state.U` | $$(4,\, n_r^\text{tot},\, n_z^\text{tot})$$ | Conservative state at cell centres |
| `state.U0` | same | Snapshot for SSP-RK |
| `state.L` | same | Spatial residual buffer |
| `state.Fr` | $$(4,\, n_r^\text{tot}+1,\, n_z^\text{tot})$$ | $$r$$-direction numerical fluxes at radial faces |
| `state.Fz` | $$(4,\, n_r^\text{tot},\, n_z^\text{tot}+1)$$ | $$z$$-direction numerical fluxes at axial faces |
| `state.lambda` | $$(n_r^\text{tot},\, n_z^\text{tot})$$ | Reaction-progress variable (programmed burn) |
| `state.rhoY` | $$(n_r^\text{tot},\, n_z^\text{tot})$$ | Conservative passive scalar (material tag) |
| `state.rhoY0` | same | RK snapshot of `rhoY` |
| `state.LrhoY` | same | Residual buffer for `rhoY` |

The conservative state is laid out with the variable index *first*
(`U[k, i, j]`) so that all four components of a single cell are contiguous
in memory; this is the layout that the HLLC and MUSCL kernels access in
inner loops. The fluxes are stored at faces, so `Fr[:, iface, j]` is the
flux through the face between cells `iface - 1` and `iface`.

Both the conservative state $$\mathbf{U}$$ and the passive scalar
$$\rho Y$$ have ghost layers and participate in halo exchange. The
reaction-progress variable $$\lambda$$ is also halo-exchanged after each
update (it is needed on neighbour ghost cells for HLLC on inter-rank
faces).
