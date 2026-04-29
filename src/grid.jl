# Local (per-rank) structured (r, z) grid with two layers of ghost cells.
#
# Each rank owns a Cartesian sub-block of the global mesh. The grid struct
# stores the local cell counts and the global indices of the rank's first
# interior cell. Cell-centre and face coordinates are stored at their
# physical locations so ghost cells naturally fall at the right place
# (negative-r when the rank borders the symmetry axis, etc.).
#
# Indexing convention
# -------------------
#   * Interior cells in radial direction: i in (ng+1)..(ng+nr)
#   * Interior cells in axial direction:  j in (ng+1)..(ng+nz)
#   * Total array width: nri = nr + 2*ng,   nzj = nz + 2*ng
#   * Face iface lies between cells (iface-1) and (iface).

const NGHOST = 2

struct Grid
    nr::Int                      # local interior cell count in r
    nz::Int                      # local interior cell count in z
    nr_global::Int
    nz_global::Int
    i_offset::Int                # global index of (local i = ng+1)  - 1
    j_offset::Int                # global index of (local j = ng+1)  - 1
    ng::Int
    nri::Int                     # nr + 2*ng
    nzj::Int                     # nz + 2*ng
    dr::Float64
    dz::Float64
    r_centers::Vector{Float64}   # length nri (local, in physical units)
    z_centers::Vector{Float64}   # length nzj
    r_faces::Vector{Float64}     # length nri+1
    z_faces::Vector{Float64}     # length nzj+1
    r_min_global::Float64
    r_max_global::Float64
    z_min_global::Float64
    z_max_global::Float64
end

"""
    build_grid(cfg, topo) -> Grid

Build a per-rank grid from the global config and the parallel topology.
The radial and axial cell partitions come from `block_extents`, so each
rank owns a contiguous block of interior cells.
"""
function build_grid(cfg::CaseConfig, topo::Topology)
    ng = NGHOST
    Pr, Pz = topo.dims
    cr, cz = topo.coords

    i_off, nr_loc = block_extents(cfg.nr, Pr, cr)
    j_off, nz_loc = block_extents(cfg.nz, Pz, cz)
    nr_loc >= 1 && nz_loc >= 1 ||
        error("Process grid $(Pr)x$(Pz) is too large for mesh $(cfg.nr)x$(cfg.nz); use fewer ranks")

    dr = (cfg.r_max - cfg.r_min) / cfg.nr
    dz = (cfg.z_max - cfg.z_min) / cfg.nz

    nri = nr_loc + 2*ng
    nzj = nz_loc + 2*ng

    r_centers = Vector{Float64}(undef, nri)
    z_centers = Vector{Float64}(undef, nzj)
    r_faces   = Vector{Float64}(undef, nri + 1)
    z_faces   = Vector{Float64}(undef, nzj + 1)

    # Local index i corresponds to global cell index (i_off + (i - ng))
    # whose centre is at r_min_global + (global_index - 0.5) * dr.
    @inbounds for i in 1:nri
        gi = i_off + (i - ng)         # 1-based global cell index
        r_centers[i] = cfg.r_min + (gi - 0.5) * dr
    end
    @inbounds for i in 1:(nri+1)
        gi_face = i_off + (i - ng) - 1     # 0-based global face index from r_min_global
        r_faces[i] = cfg.r_min + gi_face * dr
    end
    @inbounds for j in 1:nzj
        gj = j_off + (j - ng)
        z_centers[j] = cfg.z_min + (gj - 0.5) * dz
    end
    @inbounds for j in 1:(nzj+1)
        gj_face = j_off + (j - ng) - 1
        z_faces[j] = cfg.z_min + gj_face * dz
    end

    return Grid(nr_loc, nz_loc, cfg.nr, cfg.nz,
                i_off, j_off, ng, nri, nzj, dr, dz,
                r_centers, z_centers, r_faces, z_faces,
                cfg.r_min, cfg.r_max, cfg.z_min, cfg.z_max)
end

@inline interior_i(g::Grid) = (g.ng + 1):(g.ng + g.nr)
@inline interior_j(g::Grid) = (g.ng + 1):(g.ng + g.nz)

@inline rcell(g::Grid, i::Int) = g.r_centers[i]
@inline rface(g::Grid, iface::Int) = g.r_faces[iface]
