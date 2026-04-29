---
title: MPI Parallelisation
nav_order: 10
description: "MPI Cartesian topology, geometric block decomposition, halo exchange, and global reductions."
---

# MPI Parallelisation
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Decomposition strategy

COSMO is parallelised by **structured geometric block decomposition** of
the global $$n_r \times n_z$$ mesh across a $$P_r \times P_z$$ Cartesian
process grid:

* the user supplies the total rank count via `mpiexec -n N`;
* `MPI.Dims_create(N, [0, 0])` chooses the factorisation $$P_r \times P_z$$
  closest to square;
* a Cartesian communicator is built with
  `MPI.Cart_create(comm, [P_r, P_z]; periodic=[false, false])` — there is
  no periodicity in either direction;
* `MPI.Cart_shift` returns the four nearest-neighbour ranks
  (`nbr_left`, `nbr_right`, `nbr_down`, `nbr_up`), with edge ranks
  receiving `MPI.PROC_NULL` for the missing direction(s).

The implementation lives in `build_topology` in
[`src/parallel.jl`](https://github.com/jman87/cosmo/blob/main/src/parallel.jl).
Single-rank runs reduce to a $$1 \times 1$$ process grid with all
neighbours `PROC_NULL`, so the same code path serves serial and parallel
runs without a separate code branch.

This is the simplest scalable decomposition for a structured uniform
grid:

* Communication is **only** with the four immediate neighbours, so the
  per-step communication volume scales as the local boundary length
  (square-root of the local cell count) rather than as the number of
  ranks.
* Global collectives are needed only for the CFL reduction (one
  `Allreduce` per step) and the I/O gather (one `Gatherv` per output
  field per snapshot, far less frequent than per-step).
* No graph partitioning library is required (in contrast to
  unstructured-mesh codes that depend on SCOTCH or ParMETIS).

---

## Block extents and load balance

For each direction the global cell count is split with the helper

```
block_extents(n_global, np, ip) -> (i_lo, n_local)
```

defined in `src/parallel.jl`. Block `ip` (0-indexed) gets cells

$$
\bigl[\,i_\text{lo} + 1,\;\; i_\text{lo} + n_\text{local}\,\bigr]
$$

of the global grid. With $$\text{base} = \lfloor n_\text{global} / n_p \rfloor$$
and $$\text{rem} = n_\text{global} - \text{base}\,n_p$$, the partition is

$$
n_\text{local} \;=\;
\begin{cases}
  \text{base} + 1, & i_p < \text{rem}, \\
  \text{base},     & i_p \;\geq\; \text{rem},
\end{cases}
$$

so the **load imbalance is at most one cell per rank in each direction**.
The lowest-coordinate ranks absorb the remainder, which makes the
distribution deterministic and reproducible across runs.

The same partitioner is used by the output writer when reassembling the
gathered field arrays into the global rectilinear grid
([Output](output)).

---

## Per-rank state

Each rank allocates **only its block plus 2 ghost layers on every side**
(see [Mesh and Grid](mesh-discretization)). The interior cell counts
are $$n_r^\text{loc},\,n_z^\text{loc}$$ and the total array size is
$$(n_r^\text{loc} + 2 n_g) \times (n_z^\text{loc} + 2 n_g)$$. The
ghost layers play three roles:

1. They supply the 4-cell MUSCL stencil at the first/last interior
   face (neighbouring interior cells, on-edge BC values, or off-rank
   neighbour values — the residual does not care about the source);
2. They store reflected values for the symmetry-axis and rigid-ground
   BCs;
3. They store off-rank neighbour state copied in by the halo exchange.

There is **no global state array**. The only places where the full
global mesh is materialised are the output writer (transiently, on rank 0)
and the per-rank face coordinates in `OutputWriter` (a tiny
$$\mathcal{O}(n_r + n_z)$$ vector, identical on every rank).

---

## Halo exchange

The halo exchange fills the four ghost slabs of the conservative state,
the burn fraction, and the passive scalar from the corresponding interior
slabs of the four neighbour ranks. The implementation lives in
[`src/parallel_halo.jl`](https://github.com/jman87/cosmo/blob/main/src/parallel_halo.jl).

### Conservative-state exchange

Each face of the rank-local block exchanges a slab of size
$$n_g$$ in the boundary-normal direction times the full perpendicular
extent. For the $$r$$ direction the slab is $$(4,\, n_g,\, n_z^\text{tot})$$
and for the $$z$$ direction it is $$(4,\, n_r^\text{tot},\, n_g)$$.

The exchange algorithm is:

1. **Pack** the relevant interior columns into a contiguous Float64
   buffer (one for each direction × side, four total in each spatial
   direction).
2. Call `MPI.Sendrecv!` to swap each pair of buffers between adjacent
   ranks. The `Sendrecv!` form serves two roles: it issues a matched
   send and receive in one call (avoiding deadlock without explicit
   non-blocking primitives), and it makes the operation a no-op when
   one side is `PROC_NULL` (an edge rank).
3. **Unpack** the received buffer into the appropriate ghost layer of
   the rank-local array.

The four exchanges per direction (left-send/right-recv,
right-send/left-recv, etc.) use distinct message tags
(`sendtag = 10, 11, 20, 21, …`) so that the directions do not race
in MPI implementations that use tag-based message matching. The same
pattern repeats for the $$z$$ direction with shape
$$(4,\, n_r^\text{tot},\, n_g)$$ and tags 20 / 21.

### Scalar-field exchange

A single Float64 field (the burn fraction $$\lambda$$ or the passive
scalar $$\rho Y$$) uses an analogous but smaller exchange in
`exchange_halo_scalar!` with shape $$(n_g,\, n_z^\text{tot})$$ /
$$(n_r^\text{tot},\, n_g)$$ and tags 30 / 31 / 40 / 41.

### Order within an RK substep

The residual call `rhs!` (in `src/timeint.jl`) drives the exchange
sequence:

```
exchange_halos!(state.U, grid, topo)              # conservative state
exchange_halo_scalar!(state.rhoY, grid, topo)     # passive scalar
apply_bcs!(state.U, grid, topo; ...)              # physical-edge ghosts only
apply_scalar_bcs!(state.rhoY, grid, topo)
compute_fluxes!(state, grid)                      # MUSCL + HLLC at every interior face
compute_residual!(state, grid)                    # FV divergence + geometric source
```

The "exchange first, BC second" ordering ensures that internal-seam
ghosts are filled by the neighbour exchange and physical-edge ghosts
are filled by the BC routine, with no double-write.

The burn fraction $$\lambda$$ is exchanged separately at the top of the
outer time step (after `update_burn!`) rather than in `rhs!`, because
$$\lambda$$ is held constant across the SSP-RK substages and a single
exchange per outer step suffices.

---

## Global reductions

Two MPI collectives are used at every step.

### CFL reduction

The CFL bound (see [Time Integration](time-integration)) requires the
global maximum of $$|\mathbf{u}| + c$$ across all interior cells. Each
rank computes its local maximum and then participates in

```
smax_global = MPI.Allreduce(smax_local, MPI.MAX, topo.comm)
```

A second `Allreduce` (also `MAX`) gathers the global maximum of the
burn fraction so that every rank can apply the detonation-velocity
safeguard with a consistent `lam_max_global > 0` test. Both reductions
are over a single Float64 — negligible cost compared to the per-step
compute.

### Output gather (less frequent)

When a snapshot is written, each rank packs its interior cells of each
output field into a flat buffer and calls `MPI.Gatherv!` with a
`VBuffer` containing per-rank receive sizes. Rank 0 receives the
concatenated stream and reassembles it into the global rectilinear
field array using the same `block_extents` partitioner. See
[Output](output) for the full description.

The gather happens once per output cadence (`output.frequency` time
steps), so its cost is amortised over many compute steps and does not
dominate the runtime.

---

## Single-rank fallback

A serial run uses

```
mpiexec -n 1 julia --project=cosmo run.jl input.jsonc
```

or simply

```
julia --project=cosmo run.jl input.jsonc
```

(MPI initialisation is unconditional; `MPI.Comm_size = 1` reduces to the
serial path automatically.) On a single rank:

* `MPI.Dims_create` returns $$(1, 1)$$;
* All four neighbour ranks are `MPI.PROC_NULL`, so each `Sendrecv!`
  is a no-op;
* `apply_bcs!` is invoked on every edge, so all ghost layers are
  filled by the physical-boundary BCs;
* `Allreduce` is the identity on a single value;
* `Gatherv!` simply copies the local buffer.

The serial code path is therefore exactly the same code as the parallel
code path; there is no separate non-MPI build target.

---

## Performance considerations

* The dominant cost per step is the residual assembly (MUSCL + HLLC on
  every interior face, plus the divergence-theorem update on every
  cell). Communication accounts for a small fraction of the runtime
  on any reasonable rank count for blast-class problems.
* The halo-exchange buffers are allocated freshly inside each call
  rather than reused. For very high rank counts on very small per-rank
  blocks this allocation can dominate; the simpler reallocation path
  is retained for code clarity, with the understanding that production
  blast runs use enough cells per rank that the per-call allocation is
  negligible.
* The `MPI.jl` JLL build provides a working MPI runtime out of the box.
  Switching to a system MPI (e.g. OpenMPI from Homebrew) is a one-time
  configuration:

```bash
julia --project=cosmo -e 'using MPIPreferences; MPIPreferences.use_system_binary()'
```

This is occasionally useful when the JLL build's MPI is incompatible
with the system batch scheduler or shared filesystem.
