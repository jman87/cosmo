# MPI parallelisation: topology + geometric block partitioning.
#
# Decomposition
# -------------
# The (nr, nz) global grid is split into a Px-by-Pz Cartesian process grid
# whose dimensions MPI picks via `MPI.Dims_create`. Each rank owns a
# rectangular block of interior cells with `ng` ghost layers on every
# side. The same ghost layers serve both the physical-boundary BCs (filled
# locally on edge ranks) and inter-rank halo exchange (filled by
# `exchange_halos!` from neighbouring ranks; see `parallel_halo.jl`).

struct Topology
    comm::MPI.Comm
    rank::Int
    nproc::Int
    dims::NTuple{2,Int}            # (Pr, Pz)
    coords::NTuple{2,Int}          # this rank's (cr, cz) in the process grid
    nbr_left::Cint                 # rank to the -r side  (or MPI.PROC_NULL)
    nbr_right::Cint                # rank to the +r side
    nbr_down::Cint                 # rank to the -z side
    nbr_up::Cint                   # rank to the +z side
end

"""
    build_topology(comm) -> Topology

Build a 2-D Cartesian topology over `comm` whose dimensions are chosen by
`MPI.Dims_create`. Single-rank runs receive a 1x1 grid with PROC_NULL
neighbours, so the same code path covers serial and parallel.
"""
function build_topology(comm::MPI.Comm)
    nproc = MPI.Comm_size(comm)
    dims  = MPI.Dims_create(nproc, [0, 0])         # MPI picks the factorisation
    Pr    = Int(dims[1])
    Pz    = Int(dims[2])

    cart = MPI.Cart_create(comm, [Pr, Pz]; periodic=[false, false], reorder=true)
    cr, cz = MPI.Cart_coords(cart)

    nbr_left,  nbr_right = MPI.Cart_shift(cart, 0, 1)
    nbr_down,  nbr_up    = MPI.Cart_shift(cart, 1, 1)

    return Topology(cart, Int(MPI.Comm_rank(cart)), nproc,
                    (Pr, Pz), (Int(cr), Int(cz)),
                    Cint(nbr_left), Cint(nbr_right),
                    Cint(nbr_down), Cint(nbr_up))
end

"""
    block_extents(n_global, np, ip) -> (i_lo, n_local)

Geometric block partition of `n_global` cells into `np` blocks. Block
`ip` (0-indexed) gets cells `[i_lo+1 : i_lo+n_local]` of the global grid.
Remainder cells are spread over the lowest-coordinate ranks so the load
imbalance is at most one cell per rank.
"""
@inline function block_extents(n_global::Int, np::Int, ip::Int)
    base = div(n_global, np)
    rem  = n_global - base * np
    n_local = base + (ip < rem ? 1 : 0)
    i_lo    = ip * base + min(ip, rem)
    return i_lo, n_local
end

# Edge tests for applying physical BCs on the rank that owns the edge.
@inline at_r_min_edge(t::Topology) = t.coords[1] == 0
@inline at_r_max_edge(t::Topology) = t.coords[1] == t.dims[1] - 1
@inline at_z_min_edge(t::Topology) = t.coords[2] == 0
@inline at_z_max_edge(t::Topology) = t.coords[2] == t.dims[2] - 1
