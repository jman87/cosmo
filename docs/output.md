---
title: Output
nav_order: 11
description: "Cell-centred VTK rectilinear-grid snapshots, MPI gather to rank 0, master .pvd time series, and ParaView post-processing."
---

# Output
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Overview

At user-controlled intervals the solver writes the current cell-centred
state to disk in the **VTK XML rectilinear-grid** format
(`.vtr`). All snapshots are indexed by a single master `.pvd` (ParaView
Data) collection file that records the simulation time of each
snapshot and is the only file the user opens in ParaView.

The output stage is implemented in
[`src/output.jl`](https://github.com/jman87/cosmo/blob/main/src/output.jl)
and uses the [`WriteVTK.jl`](https://github.com/JuliaIO/WriteVTK.jl)
package for the on-disk encoding.

---

## Cell-centred fields

The solver stores conservative quantities as cell-volume averages, so
the VTK output uses **`VTKCellData()`** rather than point data. Each
snapshot writes the following fields per interior cell:

| Field name | Units | Definition |
|------------|-------|------------|
| `density`         | $$\text{lbf}\cdot\text{s}^2/\text{in}^4$$ | $$\rho$$ |
| `velocity_r`      | $$\text{in}/\text{s}$$ | $$u_r = (\rho u_r) / \rho$$ |
| `velocity_z`      | $$\text{in}/\text{s}$$ | $$u_z = (\rho u_z) / \rho$$ |
| `pressure`        | psi | EOS-blended $$p(\rho,\,e;\,\lambda)$$ |
| `specific_energy` | $$\text{in}^2/\text{s}^2$$ | $$e = (\rho E)/\rho - \tfrac{1}{2}(u_r^2 + u_z^2)$$ |
| `burn_fraction`   | — | $$\lambda \in [0, 1]$$ |
| `material_id`     | — | $$Y = \rho Y / \max(\rho,\,\rho_\text{floor})$$, clipped to $$[0, 1]$$ |

`material_id` is the explosive-origin mass fraction (passive scalar)
and tracks where products gas is, regardless of how far it has expanded
into the air; `burn_fraction` is the reaction-progress variable which
is fixed at 1 inside the original charge geometry once the burn front
has swept it ([Programmed Burn Model](burn-model)).

There is **no** Q2 → Q1 interpolation step (as in the original FEM
sibling solver) — the cell-centred FV state is inherently piecewise-
constant per cell, which maps directly onto the VTK cell-data layout.

---

## Rectilinear grid

A `.vtr` file represents the mesh as a **rectilinear grid**: a
1-D vector of $$r$$-coordinates ($$n_r + 1$$ values) and a 1-D vector of
$$z$$-coordinates ($$n_z + 1$$ values), with cells inferred by tensor
product. This is the most compact representation for COSMO's structured
uniform grid; it produces files perhaps an order of magnitude smaller
than an unstructured-grid VTK encoding of the same data.

The face-coordinate vectors are computed once at writer construction
(`build_writer` in `src/output.jl`):

$$
r_i^\text{face} \;=\; r_\text{min} + (i - 1)\,\Delta r, \quad i \in [1,\, n_r + 1],
$$

$$
z_j^\text{face} \;=\; z_\text{min} + (j - 1)\,\Delta z, \quad j \in [1,\, n_z + 1],
$$

with the leftmost coordinate clamped to zero when $$r_\text{min} = 0$$ to
avoid a tiny negative value from floating-point round-off.

---

## MPI gather

Each rank computes its **local** primitive fields into a flat
$$n_r^\text{loc} \cdot n_z^\text{loc}$$ buffer (one buffer per output
field). Rank 0 then collects the seven per-cell fields with one
`MPI.Gatherv!` per field:

* The per-rank receive sizes are pre-computed in `build_writer` from
  `block_extents`, indexed by the rank's Cartesian coordinates;
* `MPI.VBuffer(recv, sizes)` gathers the variable-size contributions
  into a single contiguous vector;
* Rank 0 reassembles the gathered vector into a global
  $$n_r \times n_z$$ array by walking the rank ordering on the
  Cartesian communicator and copying each rank's block into its
  $$(c_r, c_z)$$ slot.

The implementation lives in `_gather_field` in `src/output.jl`. Only
rank 0 holds the full global state; non-root ranks simply call
`Gatherv!` with `nothing` for the receive buffer and return.

This serial-on-rank-0 design is **not parallel I/O**, but for the
modest output sizes this solver targets (a 256 × 256 mesh produces a
few megabytes per snapshot) the gather cost is negligible compared to
the time-stepping cost between snapshots. The architecture trades
parallel-I/O complexity (and the dependency on a parallel HDF5 / NetCDF
build) for a simple serial-write path that produces files compatible
with every ParaView version.

---

## File names and master `.pvd`

Given `output.directory = "output"` and `output.pvd_name = "cosmo"`:

* Master file: `output/cosmo.pvd` — a small XML index that maps
  simulation times to snapshot filenames. ParaView opens *only* this
  file and traverses the index to load the time series.
* Per-snapshot file: `output/cosmo_NNNNNN.vtr` where `NNNNNN` is the
  zero-padded step number — one file per snapshot, written by rank 0.

The master `.pvd` is a `WriteVTK.CollectionFile` opened at the start of
the run with

```julia
pvd = paraview_collection(joinpath(out_dir, pvd_name))
```

and finalised with `vtk_save(pvd)` at the end of the run
(`close_writer!`). Each `write_snapshot!` call appends an entry to the
collection of the form `pvd[t] = vtk` after the `vtk_grid` block writes
the per-step `.vtr`.

If the run terminates early (CFL collapse, NaN, etc.) the writer still
finalises the `.pvd` in the `finally` block of the run loop, so the
partial time series is loadable in ParaView.

---

## Output frequency

A snapshot is written every `output.frequency` time steps. The
**initial state** (step 0) and the **final state** are always written,
regardless of the frequency setting. For a simulation with $$N$$ total
time steps and frequency $$f$$, the number of `.vtr` files written is
approximately

$$
N_\text{snap} \;\approx\; \left\lfloor\frac{N}{f}\right\rfloor + 2.
$$

Because $$\Delta t$$ is set adaptively by the CFL bound, the snapshot
*times* are not uniformly spaced — they fall on whichever simulation
times correspond to step multiples of $$f$$. ParaView's time slider
interpolates between these times when scrubbing the timeline.

For runs that need uniform output cadence in time (e.g. for Fourier
analysis of pressure histories), the user must set
`output.frequency = 1` and post-process by interpolating onto a uniform
time grid. This is rarely needed for blast applications.

---

## ParaView usage

1. **Open** `output/cosmo.pvd` in ParaView (`File → Open`). The time
   series loads automatically.
2. **Apply** in the Properties panel.
3. Use the **Play** button or the timeline slider to step through the
   snapshots.
4. Select a field from the **Coloring** drop-down (top toolbar) to
   visualise it. `pressure` and `material_id` are the most useful for
   blast diagnostics.
5. **Plot Over Line** or **Probe Location** to extract 1-D radial /
   axial profiles or point histories.

For visualising the full meridional plane and a 3-D representation of
the axisymmetric solution:

* Apply the **Reflect** filter (axis = $$z$$, plane normal = $$\hat{r}$$)
  to mirror the quarter-plane across the symmetry axis.
* Apply **Rotational Extrusion** (axis = $$z$$, sweep angle 360°) to
  generate the full 3-D blast cloud.

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---------|--------------|-----|
| Only one snapshot in `.pvd` | Run terminated before writing the final snapshot | Check the runner output for a non-finite-`dt` warning |
| Missing fields in ParaView | Older ParaView version | Upgrade to ParaView ≥ 5.10 |
| Coloring shows zero everywhere | Wrong field selected (e.g. `material_id` at $$t = 0$$ outside the charge) | Pick `density` or `pressure` |
| `.vtr` files exist but `.pvd` does not | Run terminated before `close_writer!` | The partial files are still loadable individually; combine them by hand into a manual `.pvd` if needed |
