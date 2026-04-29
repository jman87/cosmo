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

"""
    exchange_halos!(U, grid, topo)

Fill the four ghost-cell halos of `U` from neighbouring ranks (if any).
Ghost layers that face physical boundaries are not touched here — they
are filled by `apply_bcs!`.
"""
function exchange_halos!(U::Array{Float64,3}, grid::Grid, topo::Topology)
    ng  = grid.ng
    nri = grid.nri
    nzj = grid.nzj
    nr  = grid.nr
    nz  = grid.nz

    # ---- r direction --------------------------------------------------
    rsize = 4 * ng * nzj
    send_left  = Vector{Float64}(undef, rsize)
    recv_left  = Vector{Float64}(undef, rsize)
    send_right = Vector{Float64}(undef, rsize)
    recv_right = Vector{Float64}(undef, rsize)

    @inbounds for j in 1:nzj, i in 1:ng, k in 1:4
        send_left[k + 4*((i-1) + ng*(j-1))]  = U[k, ng + i, j]
        send_right[k + 4*((i-1) + ng*(j-1))] = U[k, nr + i, j]
    end

    MPI.Sendrecv!(send_left,  recv_right, topo.comm; dest=topo.nbr_left,  source=topo.nbr_right, sendtag=10, recvtag=10)
    MPI.Sendrecv!(send_right, recv_left,  topo.comm; dest=topo.nbr_right, source=topo.nbr_left,  sendtag=11, recvtag=11)

    if topo.nbr_left != MPI.PROC_NULL
        @inbounds for j in 1:nzj, i in 1:ng, k in 1:4
            U[k, i, j] = recv_left[k + 4*((i-1) + ng*(j-1))]
        end
    end
    if topo.nbr_right != MPI.PROC_NULL
        @inbounds for j in 1:nzj, i in 1:ng, k in 1:4
            U[k, nr + ng + i, j] = recv_right[k + 4*((i-1) + ng*(j-1))]
        end
    end

    # ---- z direction --------------------------------------------------
    zsize = 4 * nri * ng
    send_down = Vector{Float64}(undef, zsize)
    recv_down = Vector{Float64}(undef, zsize)
    send_up   = Vector{Float64}(undef, zsize)
    recv_up   = Vector{Float64}(undef, zsize)

    @inbounds for j in 1:ng, i in 1:nri, k in 1:4
        send_down[k + 4*((i-1) + nri*(j-1))] = U[k, i, ng + j]
        send_up[k + 4*((i-1) + nri*(j-1))]   = U[k, i, nz + j]
    end

    MPI.Sendrecv!(send_down, recv_up,   topo.comm; dest=topo.nbr_down, source=topo.nbr_up,   sendtag=20, recvtag=20)
    MPI.Sendrecv!(send_up,   recv_down, topo.comm; dest=topo.nbr_up,   source=topo.nbr_down, sendtag=21, recvtag=21)

    if topo.nbr_down != MPI.PROC_NULL
        @inbounds for j in 1:ng, i in 1:nri, k in 1:4
            U[k, i, j] = recv_down[k + 4*((i-1) + nri*(j-1))]
        end
    end
    if topo.nbr_up != MPI.PROC_NULL
        @inbounds for j in 1:ng, i in 1:nri, k in 1:4
            U[k, i, nz + ng + j] = recv_up[k + 4*((i-1) + nri*(j-1))]
        end
    end
    return U
end

"""
    exchange_halo_scalar!(F, grid, topo)

Halo exchange for a single scalar field stored as `F[i, j]`. Used for the
burn-fraction array, which is a per-cell value blended into the EOS.
"""
function exchange_halo_scalar!(F::Array{Float64,2}, grid::Grid, topo::Topology)
    ng  = grid.ng
    nri = grid.nri
    nzj = grid.nzj
    nr  = grid.nr
    nz  = grid.nz

    rsize = ng * nzj
    sl = Vector{Float64}(undef, rsize); rl = Vector{Float64}(undef, rsize)
    sr = Vector{Float64}(undef, rsize); rr = Vector{Float64}(undef, rsize)
    @inbounds for j in 1:nzj, i in 1:ng
        sl[i + ng*(j-1)] = F[ng + i, j]
        sr[i + ng*(j-1)] = F[nr + i, j]
    end
    MPI.Sendrecv!(sl, rr, topo.comm; dest=topo.nbr_left,  source=topo.nbr_right, sendtag=30, recvtag=30)
    MPI.Sendrecv!(sr, rl, topo.comm; dest=topo.nbr_right, source=topo.nbr_left,  sendtag=31, recvtag=31)
    if topo.nbr_left != MPI.PROC_NULL
        @inbounds for j in 1:nzj, i in 1:ng
            F[i, j] = rl[i + ng*(j-1)]
        end
    end
    if topo.nbr_right != MPI.PROC_NULL
        @inbounds for j in 1:nzj, i in 1:ng
            F[nr + ng + i, j] = rr[i + ng*(j-1)]
        end
    end

    zsize = nri * ng
    sd = Vector{Float64}(undef, zsize); rd = Vector{Float64}(undef, zsize)
    su = Vector{Float64}(undef, zsize); ru = Vector{Float64}(undef, zsize)
    @inbounds for j in 1:ng, i in 1:nri
        sd[i + nri*(j-1)] = F[i, ng + j]
        su[i + nri*(j-1)] = F[i, nz + j]
    end
    MPI.Sendrecv!(sd, ru, topo.comm; dest=topo.nbr_down, source=topo.nbr_up,   sendtag=40, recvtag=40)
    MPI.Sendrecv!(su, rd, topo.comm; dest=topo.nbr_up,   source=topo.nbr_down, sendtag=41, recvtag=41)
    if topo.nbr_down != MPI.PROC_NULL
        @inbounds for j in 1:ng, i in 1:nri
            F[i, j] = rd[i + nri*(j-1)]
        end
    end
    if topo.nbr_up != MPI.PROC_NULL
        @inbounds for j in 1:ng, i in 1:nri
            F[i, nz + ng + j] = ru[i + nri*(j-1)]
        end
    end
    return F
end
