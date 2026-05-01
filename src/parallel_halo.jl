# Halo exchange for the conservative state and the burn-fraction field.
#
# For a 4-component conservative state stored as U[k, i, j], the halo
# along the r direction is a slab of size (4, ng, nzj); along z it is
# (4, nri, ng). Each exchange packs the relevant interior columns into a
# contiguous Float64 buffer, calls `MPI.Sendrecv!`, and unpacks into the
# adjacent ghost layer. Ranks with no neighbour in a given direction
# receive PROC_NULL and the sendrecv is a no-op there.
#
# Loaded after grid.jl so the Grid struct is in scope for the type
# annotations on the function signatures.

# Pre-allocated communication buffers for halo exchange.
# Storing them here avoids heap allocation on the hot path (rhs! is called
# twice per SSP-RK2 step). Scalar buffers are shared between successive
# exchange_halo_scalar! calls, which is safe because they are sequential.
struct HaloBuffers
    # 4-component (U) — r direction, size 4*ng*nzj each
    u_send_left  ::Vector{Float64}
    u_recv_left  ::Vector{Float64}
    u_send_right ::Vector{Float64}
    u_recv_right ::Vector{Float64}
    # 4-component (U) — z direction, size 4*nri*ng each
    u_send_down  ::Vector{Float64}
    u_recv_down  ::Vector{Float64}
    u_send_up    ::Vector{Float64}
    u_recv_up    ::Vector{Float64}
    # scalar — r direction, size ng*nzj each
    s_send_left  ::Vector{Float64}
    s_recv_left  ::Vector{Float64}
    s_send_right ::Vector{Float64}
    s_recv_right ::Vector{Float64}
    # scalar — z direction, size nri*ng each
    s_send_down  ::Vector{Float64}
    s_recv_down  ::Vector{Float64}
    s_send_up    ::Vector{Float64}
    s_recv_up    ::Vector{Float64}
end

function build_halo_buffers(grid::Grid)
    ng  = grid.ng
    nri = grid.nri
    nzj = grid.nzj
    ur  = 4 * ng * nzj   # 4-component r-slab
    uz  = 4 * nri * ng   # 4-component z-slab
    sr  = ng * nzj        # scalar r-slab
    sz  = nri * ng        # scalar z-slab
    return HaloBuffers(
        Vector{Float64}(undef, ur), Vector{Float64}(undef, ur),
        Vector{Float64}(undef, ur), Vector{Float64}(undef, ur),
        Vector{Float64}(undef, uz), Vector{Float64}(undef, uz),
        Vector{Float64}(undef, uz), Vector{Float64}(undef, uz),
        Vector{Float64}(undef, sr), Vector{Float64}(undef, sr),
        Vector{Float64}(undef, sr), Vector{Float64}(undef, sr),
        Vector{Float64}(undef, sz), Vector{Float64}(undef, sz),
        Vector{Float64}(undef, sz), Vector{Float64}(undef, sz),
    )
end

"""
    exchange_halos!(U, grid, topo, hbufs)

Fill the four ghost-cell halos of `U` from neighbouring ranks (if any).
Ghost layers that face physical boundaries are not touched here — they
are filled by `apply_bcs!`.
"""
function exchange_halos!(U::Array{Float64,3}, grid::Grid, topo::Topology,
                         hbufs::HaloBuffers)
    ng  = grid.ng
    nri = grid.nri
    nzj = grid.nzj
    nr  = grid.nr
    nz  = grid.nz

    # ---- r direction --------------------------------------------------
    @inbounds for j in 1:nzj, i in 1:ng, k in 1:4
        hbufs.u_send_left[k + 4*((i-1) + ng*(j-1))]  = U[k, ng + i, j]
        hbufs.u_send_right[k + 4*((i-1) + ng*(j-1))] = U[k, nr + i, j]
    end

    MPI.Sendrecv!(hbufs.u_send_left,  hbufs.u_recv_right, topo.comm; dest=topo.nbr_left,  source=topo.nbr_right, sendtag=10, recvtag=10)
    MPI.Sendrecv!(hbufs.u_send_right, hbufs.u_recv_left,  topo.comm; dest=topo.nbr_right, source=topo.nbr_left,  sendtag=11, recvtag=11)

    if topo.nbr_left != MPI.PROC_NULL
        @inbounds for j in 1:nzj, i in 1:ng, k in 1:4
            U[k, i, j] = hbufs.u_recv_left[k + 4*((i-1) + ng*(j-1))]
        end
    end
    if topo.nbr_right != MPI.PROC_NULL
        @inbounds for j in 1:nzj, i in 1:ng, k in 1:4
            U[k, nr + ng + i, j] = hbufs.u_recv_right[k + 4*((i-1) + ng*(j-1))]
        end
    end

    # ---- z direction --------------------------------------------------
    @inbounds for j in 1:ng, i in 1:nri, k in 1:4
        hbufs.u_send_down[k + 4*((i-1) + nri*(j-1))] = U[k, i, ng + j]
        hbufs.u_send_up[k + 4*((i-1) + nri*(j-1))]   = U[k, i, nz + j]
    end

    MPI.Sendrecv!(hbufs.u_send_down, hbufs.u_recv_up,   topo.comm; dest=topo.nbr_down, source=topo.nbr_up,   sendtag=20, recvtag=20)
    MPI.Sendrecv!(hbufs.u_send_up,   hbufs.u_recv_down, topo.comm; dest=topo.nbr_up,   source=topo.nbr_down, sendtag=21, recvtag=21)

    if topo.nbr_down != MPI.PROC_NULL
        @inbounds for j in 1:ng, i in 1:nri, k in 1:4
            U[k, i, j] = hbufs.u_recv_down[k + 4*((i-1) + nri*(j-1))]
        end
    end
    if topo.nbr_up != MPI.PROC_NULL
        @inbounds for j in 1:ng, i in 1:nri, k in 1:4
            U[k, i, nz + ng + j] = hbufs.u_recv_up[k + 4*((i-1) + nri*(j-1))]
        end
    end
    return U
end

"""
    exchange_halo_scalar!(F, grid, topo, hbufs)

Halo exchange for a single scalar field stored as `F[i, j]`. Used for the
burn-fraction array, which is a per-cell value blended into the EOS.
The scalar buffers in `hbufs` are shared across successive calls and must
not be called concurrently.
"""
function exchange_halo_scalar!(F::Array{Float64,2}, grid::Grid, topo::Topology,
                                hbufs::HaloBuffers)
    ng  = grid.ng
    nri = grid.nri
    nzj = grid.nzj
    nr  = grid.nr
    nz  = grid.nz

    @inbounds for j in 1:nzj, i in 1:ng
        hbufs.s_send_left[i + ng*(j-1)]  = F[ng + i, j]
        hbufs.s_send_right[i + ng*(j-1)] = F[nr + i, j]
    end
    MPI.Sendrecv!(hbufs.s_send_left,  hbufs.s_recv_right, topo.comm; dest=topo.nbr_left,  source=topo.nbr_right, sendtag=30, recvtag=30)
    MPI.Sendrecv!(hbufs.s_send_right, hbufs.s_recv_left,  topo.comm; dest=topo.nbr_right, source=topo.nbr_left,  sendtag=31, recvtag=31)
    if topo.nbr_left != MPI.PROC_NULL
        @inbounds for j in 1:nzj, i in 1:ng
            F[i, j] = hbufs.s_recv_left[i + ng*(j-1)]
        end
    end
    if topo.nbr_right != MPI.PROC_NULL
        @inbounds for j in 1:nzj, i in 1:ng
            F[nr + ng + i, j] = hbufs.s_recv_right[i + ng*(j-1)]
        end
    end

    @inbounds for j in 1:ng, i in 1:nri
        hbufs.s_send_down[i + nri*(j-1)] = F[i, ng + j]
        hbufs.s_send_up[i + nri*(j-1)]   = F[i, nz + j]
    end
    MPI.Sendrecv!(hbufs.s_send_down, hbufs.s_recv_up,   topo.comm; dest=topo.nbr_down, source=topo.nbr_up,   sendtag=40, recvtag=40)
    MPI.Sendrecv!(hbufs.s_send_up,   hbufs.s_recv_down, topo.comm; dest=topo.nbr_up,   source=topo.nbr_down, sendtag=41, recvtag=41)
    if topo.nbr_down != MPI.PROC_NULL
        @inbounds for j in 1:ng, i in 1:nri
            F[i, j] = hbufs.s_recv_down[i + nri*(j-1)]
        end
    end
    if topo.nbr_up != MPI.PROC_NULL
        @inbounds for j in 1:ng, i in 1:nri
            F[i, nz + ng + j] = hbufs.s_recv_up[i + nri*(j-1)]
        end
    end
    return F
end
