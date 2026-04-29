# Cell-centred VTK output on a 2-D rectilinear grid, with MPI gather.
#
# Each rank computes its local primitive fields, sends them to rank 0
# with a single Gatherv per field, and rank 0 writes the global rectilinear
# .vtr snapshot plus updates the master .pvd time series. Per-snapshot I/O
# is therefore serial on rank 0 — that is fine for the modest output sizes
# this solver targets and removes any need for parallel I/O libraries.

mutable struct OutputWriter
    out_dir::String
    pvd_path::String
    pvd::Union{Nothing, WriteVTK.CollectionFile}
    base_name::String
    rs_global::Vector{Float64}
    zs_global::Vector{Float64}
    counts_r::Vector{Cint}                    # number of interior r-cells per rank-r-coord
    counts_z::Vector{Cint}
    grid::Grid
    topo::Topology
end

function build_writer(cfg::CaseConfig, grid::Grid, topo::Topology)
    if topo.rank == 0
        isdir(cfg.out_dir) || mkpath(cfg.out_dir)
    end
    MPI.Barrier(topo.comm)

    pvd_path = joinpath(cfg.out_dir, cfg.pvd_name)
    pvd = topo.rank == 0 ? paraview_collection(pvd_path) : nothing

    # Global face coordinates (computed locally; identical on every rank).
    Pr, Pz = topo.dims
    nr_g, nz_g = grid.nr_global, grid.nz_global

    rs_global = Vector{Float64}(undef, nr_g + 1)
    zs_global = Vector{Float64}(undef, nz_g + 1)
    @inbounds for i in 1:(nr_g + 1)
        rs_global[i] = cfg.r_min + (i - 1) * grid.dr
    end
    @inbounds for j in 1:(nz_g + 1)
        zs_global[j] = cfg.z_min + (j - 1) * grid.dz
    end
    rs_global[1] = max(rs_global[1], 0.0)

    # Per-rank-coord interior cell counts, used by Gatherv for the field
    # arrays (one row per (cr, cz) with shape (nr_loc, nz_loc) on rank 0).
    counts_r = Vector{Cint}(undef, Pr)
    counts_z = Vector{Cint}(undef, Pz)
    for cr in 0:(Pr - 1)
        _, n = block_extents(nr_g, Pr, cr)
        counts_r[cr + 1] = Cint(n)
    end
    for cz in 0:(Pz - 1)
        _, n = block_extents(nz_g, Pz, cz)
        counts_z[cz + 1] = Cint(n)
    end

    return OutputWriter(cfg.out_dir, pvd_path, pvd, cfg.pvd_name,
                        rs_global, zs_global, counts_r, counts_z, grid, topo)
end

# Pull a flat Vector{Float64} of interior values from any 2-D field whose
# layout matches the local-cell dimensions (lambda, primitive arrays, etc.).
@inline function _pack_interior!(buf::Vector{Float64}, F::Array{Float64,2}, grid::Grid)
    ng = grid.ng
    nr = grid.nr
    nz = grid.nz
    @inbounds for j in 1:nz, i in 1:nr
        buf[(j - 1) * nr + i] = F[i + ng, j + ng]
    end
    return buf
end

# Gather an interior field onto rank 0 as a global (nr_g x nz_g) array.
# Each rank's contribution arrives as a contiguous block; rank 0 then
# scatters the blocks back into their (cr, cz) slot.
function _gather_field(local_buf::Vector{Float64}, w::OutputWriter)
    topo = w.topo
    Pr, Pz = topo.dims
    nr_g, nz_g = w.grid.nr_global, w.grid.nz_global

    # Block sizes per rank, in flattened-row-major order matching MPI's
    # rank ordering on the Cartesian communicator (row-major over (cr, cz)).
    sizes = Vector{Cint}(undef, topo.nproc)
    for r in 0:(topo.nproc - 1)
        cr, cz = MPI.Cart_coords(topo.comm, r)
        sizes[r + 1] = w.counts_r[cr + 1] * w.counts_z[cz + 1]
    end

    if topo.rank == 0
        recv = Vector{Float64}(undef, nr_g * nz_g)
        recv_buf = MPI.VBuffer(recv, sizes)
        MPI.Gatherv!(local_buf, recv_buf, topo.comm; root=0)

        # Reassemble per-rank blocks into the global (nr_g x nz_g) array.
        out = Array{Float64}(undef, nr_g, nz_g)
        offset = 0
        for r in 0:(topo.nproc - 1)
            cr, cz = MPI.Cart_coords(topo.comm, r)
            i_off, nr_loc = block_extents(nr_g, Pr, cr)
            j_off, nz_loc = block_extents(nz_g, Pz, cz)
            for j in 1:nz_loc, i in 1:nr_loc
                out[i_off + i, j_off + j] = recv[offset + (j - 1) * nr_loc + i]
            end
            offset += nr_loc * nz_loc
        end
        return out
    else
        MPI.Gatherv!(local_buf, nothing, topo.comm; root=0)
        return nothing
    end
end

"""
    write_snapshot!(writer, state, t, step)

Append one .vtr snapshot to the .pvd time series. Density, the two
velocity components, pressure, total internal energy, and the burn
fraction are written as cell data.
"""
function write_snapshot!(w::OutputWriter, state::State, t::Float64, step::Int)
    grid = w.grid
    nr   = grid.nr
    nz   = grid.nz
    ng   = grid.ng

    nl = nr * nz
    rho_l = Vector{Float64}(undef, nl)
    ur_l  = Vector{Float64}(undef, nl)
    uz_l  = Vector{Float64}(undef, nl)
    p_l   = Vector{Float64}(undef, nl)
    e_l   = Vector{Float64}(undef, nl)
    lam_l = Vector{Float64}(undef, nl)
    mat_l = Vector{Float64}(undef, nl)

    U    = state.U
    lam  = state.lambda
    rhoY = state.rhoY
    @inbounds for j in 1:nz, i in 1:nr
        ig = i + ng
        jg = j + ng
        r, u, v, p = primitives((U[1,ig,jg], U[2,ig,jg], U[3,ig,jg], U[4,ig,jg]),
                                lam[ig, jg])
        idx = (j - 1) * nr + i
        rho_l[idx] = r
        ur_l[idx]  = u
        uz_l[idx]  = v
        p_l[idx]   = p
        e_l[idx]   = U[4, ig, jg] / r - 0.5 * (u*u + v*v)
        lam_l[idx] = lam[ig, jg]
        # Material id is the explosive-origin mass fraction Y = rho*Y / rho.
        # It is passively advected with the flow, so the products' material
        # tag propagates outward as the charge expands into the air. Clipped
        # to [0, 1] to absorb tiny round-off excursions from the donor-cell
        # advection.
        Yv = rhoY[ig, jg] / max(r, RHO_FLOOR)
        mat_l[idx] = clamp(Yv, 0.0, 1.0)
    end

    rho_g = _gather_field(rho_l, w)
    ur_g  = _gather_field(ur_l,  w)
    uz_g  = _gather_field(uz_l,  w)
    p_g   = _gather_field(p_l,   w)
    e_g   = _gather_field(e_l,   w)
    lam_g = _gather_field(lam_l, w)
    mat_g = _gather_field(mat_l, w)

    if w.topo.rank == 0
        fname = @sprintf("%s_%06d", w.base_name, step)
        vtk_path = joinpath(w.out_dir, fname)

        vtk_grid(vtk_path, w.rs_global, w.zs_global) do vtk
            vtk["density",         VTKCellData()] = rho_g
            vtk["velocity_r",      VTKCellData()] = ur_g
            vtk["velocity_z",      VTKCellData()] = uz_g
            vtk["pressure",        VTKCellData()] = p_g
            vtk["specific_energy", VTKCellData()] = e_g
            vtk["burn_fraction",   VTKCellData()] = lam_g
            vtk["material_id",     VTKCellData()] = mat_g
            w.pvd[t] = vtk
        end
    end
    MPI.Barrier(w.topo.comm)
    return nothing
end

function close_writer!(w::OutputWriter)
    if w.topo.rank == 0 && w.pvd !== nothing
        vtk_save(w.pvd)
        w.pvd = nothing
    end
    return nothing
end
