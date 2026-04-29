# Ghost-cell boundary conditions.
#
# Physical-boundary BCs are applied only on the *physical* edges of the
# global domain, identified through the parallel topology. Ghost layers
# that face an internal partition seam are filled by `exchange_halos!`
# instead and must not be overwritten here.
#
# r = r_min : reflective if r_min == 0  (axisymmetric symmetry axis)
#             zero-gradient outflow otherwise
# r = r_max : zero-gradient outflow
# z = z_min : reflective if z_min == 0  (rigid ground / symmetry plane)
#             zero-gradient outflow otherwise
# z = z_max : zero-gradient outflow

"""
    apply_bcs!(U, grid, topo; reflect_r_min, reflect_z_min)

Fill the ghost cells of `U` on physical-boundary edges only. Halo cells
on internal seams are left untouched (the halo exchange fills those).
"""
function apply_bcs!(U::Array{Float64,3}, grid::Grid, topo::Topology;
                    reflect_r_min::Bool=true,
                    reflect_z_min::Bool=true)
    ng  = grid.ng
    nr  = grid.nr
    nz  = grid.nz
    nri = grid.nri
    nzj = grid.nzj

    # ---- r = r_min ------------------------------------------------------
    if at_r_min_edge(topo)
        @inbounds for j in 1:nzj
            for g in 1:ng
                src = reflect_r_min ? (2*ng + 1 - g) : (ng + 1)
                sgn = reflect_r_min ? -1.0 : 1.0
                U[1, g, j] = U[1, src, j]
                U[2, g, j] = sgn * U[2, src, j]   # u_r flipped on reflection
                U[3, g, j] = U[3, src, j]
                U[4, g, j] = U[4, src, j]
            end
        end
    end

    # ---- r = r_max  (always outflow) ------------------------------------
    if at_r_max_edge(topo)
        @inbounds for j in 1:nzj
            for g in 1:ng
                dst = nr + ng + g
                src = nr + ng                   # last interior cell
                U[1, dst, j] = U[1, src, j]
                U[2, dst, j] = U[2, src, j]
                U[3, dst, j] = U[3, src, j]
                U[4, dst, j] = U[4, src, j]
            end
        end
    end

    # ---- z = z_min ------------------------------------------------------
    if at_z_min_edge(topo)
        @inbounds for i in 1:nri
            for g in 1:ng
                src = reflect_z_min ? (2*ng + 1 - g) : (ng + 1)
                sgn = reflect_z_min ? -1.0 : 1.0
                U[1, i, g] = U[1, i, src]
                U[2, i, g] = U[2, i, src]
                U[3, i, g] = sgn * U[3, i, src]   # u_z flipped on reflection
                U[4, i, g] = U[4, i, src]
            end
        end
    end

    # ---- z = z_max  (always outflow) ------------------------------------
    if at_z_max_edge(topo)
        @inbounds for i in 1:nri
            for g in 1:ng
                dst = nz + ng + g
                src = nz + ng
                U[1, i, dst] = U[1, i, src]
                U[2, i, dst] = U[2, i, src]
                U[3, i, dst] = U[3, i, src]
                U[4, i, dst] = U[4, i, src]
            end
        end
    end
    return U
end

"""
    apply_scalar_bcs!(F, grid, topo)

Ghost-cell BC for a passive scalar (rho*Y, etc.). Reflective and outflow
conditions both reduce to "copy the interior value" — no sign flip is
needed for a scalar quantity. Internal seams are skipped (filled by halo
exchange).
"""
function apply_scalar_bcs!(F::Array{Float64,2}, grid::Grid, topo::Topology)
    ng  = grid.ng
    nr  = grid.nr
    nz  = grid.nz
    nri = grid.nri
    nzj = grid.nzj

    if at_r_min_edge(topo)
        @inbounds for j in 1:nzj, g in 1:ng
            F[g, j] = F[2*ng + 1 - g, j]            # mirror across axis
        end
    end
    if at_r_max_edge(topo)
        @inbounds for j in 1:nzj, g in 1:ng
            F[nr + ng + g, j] = F[nr + ng, j]
        end
    end
    if at_z_min_edge(topo)
        @inbounds for i in 1:nri, g in 1:ng
            F[i, g] = F[i, 2*ng + 1 - g]
        end
    end
    if at_z_max_edge(topo)
        @inbounds for i in 1:nri, g in 1:ng
            F[i, nz + ng + g] = F[i, nz + ng]
        end
    end
    return F
end
