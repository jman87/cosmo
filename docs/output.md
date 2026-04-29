---
title: Output
nav_order: 10
description: "VTK output format, field interpolation, and ParaView post-processing."
---

# Output
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Overview

At user-controlled intervals the solver writes the current field state to disk in the
**VTK XML** format used natively by ParaView.  Each snapshot is a separate `.vtu`
(unstructured grid) file; all snapshots are indexed by a master `.pvd` (ParaView Data)
file that records the simulation time associated with each snapshot.

---

## Q2 → Q1 interpolation

The conserved-variable fields live in the degree-2 Lagrange Q2 space (9-node quad elements).
Standard VTK unstructured-grid format can represent quadratic cells via the
`VTK_LAGRANGE_QUADRILATERAL` cell type, but support for this cell type is fragile across
ParaView versions.  The solver therefore **interpolates all output fields to the degree-1
Lagrange (Q1) space** (4-node bilinear quads, one DOF per vertex) before writing:

```
Q2 field (9 DOFs/cell)  →  fem.Expression  →  Q1 field (4 DOFs/cell)
```

Q1 DOFs coincide exactly with mesh vertices, producing standard VTK bilinear quad cells
that load cleanly in any ParaView version.  The sub-cell polynomial variation of the Q2
state is lost in the interpolation; the accuracy shown in ParaView is limited to
linear variation within each cell.

---

## Output fields

The following fields can be written; the default set is noted.

| Field name | Label in VTK | Default | Units |
|------------|-------------|---------|-------|
| `density` | `density_lbf_s2_in4` | ✓ | lbf·s²/in⁴ |
| `velocity_r` | `velocity_r_in_s` | ✓ | in/s |
| `velocity_z` | `velocity_z_in_s` | ✓ | in/s |
| `pressure` | `pressure_psi` | ✓ | psi |
| `material_id` | `material_id` | ✓ | — (0=air, 1=products) |
| `energy` | `total_energy_density_lbf_in2` | — | lbf/in² |

Velocity components are recovered from the conserved momentum by $$u_r = (\rho u_r)/\rho$$.
Pressure is computed from the blended JWL/ideal-gas EOS formula evaluated at the current
$$(\rho, e)$$ state.

---

## File naming

Given `output.pvd_name = "cosmo"` and `output.directory = "output"`:

* Master file: `output/cosmo.pvd`
* Scalar snapshot (serial): `output/cosmo_p0_NNNNNN.vtu`
* Parallel snapshot (N ranks): `output/cosmoNNNNNN.pvtu` (references per-rank `.vtu` files)

---

## Output frequency

A snapshot is written every `output.frequency` time steps.  The initial state (step 0)
and the final state are always written regardless of the frequency setting.  For a
simulation with $$N$$ total steps and frequency $$f$$, approximately $$\lfloor N/f \rfloor + 2$$
snapshots are written.

---

## ParaView usage

1. Open `output/cosmo.pvd` in ParaView.
2. Click **Apply** in the Properties panel to load the time series.
3. Use the **Play** button or timeline slider to step through snapshots.
4. Select fields from the **Coloring** drop-down (top toolbar).
5. Use **Plot Over Line** or **Probe Location** to extract 1-D profiles or point histories.

For axisymmetric results, apply the **Reflect** filter (reflect about the $$z$$-axis,
plane normal $$= \hat{r}$$) to reconstruct the full meridional plane, then apply
**Rotational Extrusion** to generate a 3-D representation.

---

## MPI parallel output

All I/O is performed in parallel through DOLFINx's native `VTKFile` writer.  Each MPI
rank writes its local portion of the domain to a rank-specific `.vtu` file; the writer
generates a `.pvtu` (parallel VTU) header file that references all rank files.  No
manual gather to rank 0 is required.

The output directory is created by rank 0; a `MPI.Barrier()` ensures the directory exists
before other ranks attempt to write.
