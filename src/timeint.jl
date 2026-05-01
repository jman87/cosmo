# Explicit time integration: SSP-RK2 (Heun) and forward Euler.
#
# Spatial residual L(U) for cell (i, j) (axisymmetric finite volume):
#
#   L(U) = -1/V * [ (F_R*r_R - F_L*r_L)*dz + (G_T - G_B)*r_c*dr ]
#          + S_p
#
# where V = r_c * dr * dz, F is the r-flux, G is the z-flux, and the
# geometric pressure source is S_p = (0, p/r_c, 0, 0). The face areas are
# r-weighted (annular cell volume = 2*pi*r dr dz, divided by 2*pi).
#
# All EOS evaluations take a per-cell burn-fraction lambda so the solver
# transparently handles air (lambda = 0), TNT detonation products
# (lambda = 1), and the brief mixing region at the burn front.

struct State
    U  ::Array{Float64,3}      # (4, nri, nzj) current state
    U0 ::Array{Float64,3}      # snapshot for RK
    L  ::Array{Float64,3}      # residual buffer
    Fr ::Array{Float64,3}      # (4, nri+1, nzj) r-fluxes
    Fz ::Array{Float64,3}      # (4, nri,   nzj+1) z-fluxes
    lambda::Array{Float64,2}   # (nri, nzj) burn fraction (programmed burn,
                               #   position-based; controls EOS blend)
    rhoY  ::Array{Float64,2}   # (nri, nzj) rho * Y, where Y is the mass
                               #   fraction of explosive-origin material;
                               #   passively advected so the material tag
                               #   moves with the products as they expand
    rhoY0 ::Array{Float64,2}   # snapshot for RK
    LrhoY ::Array{Float64,2}   # residual for rhoY
    hbufs ::HaloBuffers        # pre-allocated MPI communication buffers
end

function build_state(grid::Grid)
    U  = zeros(Float64, 4, grid.nri, grid.nzj)
    U0 = similar(U)
    L  = similar(U)
    Fr = zeros(Float64, 4, grid.nri + 1, grid.nzj)
    Fz = zeros(Float64, 4, grid.nri,     grid.nzj + 1)
    lam   = zeros(Float64, grid.nri, grid.nzj)
    rhoY  = zeros(Float64, grid.nri, grid.nzj)
    rhoY0 = similar(rhoY)
    LrhoY = similar(rhoY)
    hbufs = build_halo_buffers(grid)
    return State(U, U0, L, Fr, Fz, lam, rhoY, rhoY0, LrhoY, hbufs)
end

# ---- residual assembly ---------------------------------------------------

function compute_fluxes!(state::State, grid::Grid)
    U   = state.U
    lam = state.lambda
    Fr  = state.Fr
    Fz  = state.Fz
    ng  = grid.ng

    # r-fluxes through interior + adjacent boundary faces.
    # Each j-strip writes to an independent column of Fr, so threads are safe.
    Threads.@threads for j in interior_j(grid)
        @inbounds for iface in (ng + 1):(grid.nr + ng + 1)
            i_L = iface - 1
            i_R = iface
            U_LL = (U[1,i_L-1,j], U[2,i_L-1,j], U[3,i_L-1,j], U[4,i_L-1,j])
            U_L  = (U[1,i_L,  j], U[2,i_L,  j], U[3,i_L,  j], U[4,i_L,  j])
            U_R  = (U[1,i_R,  j], U[2,i_R,  j], U[3,i_R,  j], U[4,i_R,  j])
            U_RR = (U[1,i_R+1,j], U[2,i_R+1,j], U[3,i_R+1,j], U[4,i_R+1,j])

            UL_face, UR_face = reconstruct_face(U_LL, U_L, U_R, U_RR)

            # Burn fraction is piecewise-constant (cell-centred); use the
            # owning cell's value on each side of the face. lambda only
            # changes during programmed-burn transitions, which span a
            # handful of cells over the duration of the run, so a simple
            # piecewise-constant choice is more than accurate enough.
            f = hllc_flux(UL_face, UR_face, lam[i_L, j], lam[i_R, j], 1)
            Fr[1,iface,j] = f[1]
            Fr[2,iface,j] = f[2]
            Fr[3,iface,j] = f[3]
            Fr[4,iface,j] = f[4]
        end
    end

    # z-fluxes. Each jface writes to an independent slice of Fz.
    Threads.@threads for jface in (ng + 1):(grid.nz + ng + 1)
        j_L = jface - 1
        j_R = jface
        @inbounds for i in interior_i(grid)
            U_LL = (U[1,i,j_L-1], U[2,i,j_L-1], U[3,i,j_L-1], U[4,i,j_L-1])
            U_L  = (U[1,i,j_L  ], U[2,i,j_L  ], U[3,i,j_L  ], U[4,i,j_L  ])
            U_R  = (U[1,i,j_R  ], U[2,i,j_R  ], U[3,i,j_R  ], U[4,i,j_R  ])
            U_RR = (U[1,i,j_R+1], U[2,i,j_R+1], U[3,i,j_R+1], U[4,i,j_R+1])

            UL_face, UR_face = reconstruct_face(U_LL, U_L, U_R, U_RR)
            g = hllc_flux(UL_face, UR_face, lam[i, j_L], lam[i, j_R], 2)
            Fz[1,i,jface] = g[1]
            Fz[2,i,jface] = g[2]
            Fz[3,i,jface] = g[3]
            Fz[4,i,jface] = g[4]
        end
    end
    return nothing
end

function compute_residual!(state::State, grid::Grid)
    fill!(state.L, 0.0)
    fill!(state.LrhoY, 0.0)
    Fr   = state.Fr
    Fz   = state.Fz
    L    = state.L
    LY   = state.LrhoY
    U    = state.U
    rhoY = state.rhoY
    lam  = state.lambda
    dr   = grid.dr
    dz   = grid.dz

    # Each j-strip writes to independent rows of L and LY — no race conditions.
    Threads.@threads for j in interior_j(grid)
        @inbounds for i in interior_i(grid)
            r_c = rcell(grid, i)
            r_L = rface(grid, i)        # left  face of cell i
            r_R = rface(grid, i + 1)    # right face of cell i
            invV = 1.0 / (r_c * dr * dz)

            _, _, _, p = primitives((U[1,i,j], U[2,i,j], U[3,i,j], U[4,i,j]),
                                    lam[i, j])

            for k in 1:4
                flux_r = (Fr[k, i+1, j] * r_R - Fr[k, i, j] * r_L) * dz
                flux_z = (Fz[k, i, j+1] - Fz[k, i, j]) * (r_c * dr)
                L[k, i, j] = -(flux_r + flux_z) * invV
            end

            # Axisymmetric geometric pressure source acts on r-momentum only
            L[2, i, j] += p / r_c

            # Passive-scalar advection of rho*Y. Donor-cell upwinding on
            # the HLLC mass flux (Fr[1], Fz[1]) gives a robust, monotonic
            # 1st-order update that carries the material tag with the flow.
            mfL = Fr[1, i,   j]
            mfR = Fr[1, i+1, j]
            mfB = Fz[1, i,   j]
            mfT = Fz[1, i, j+1]
            rho_im = max(U[1, i-1, j], RHO_FLOOR); Y_im = rhoY[i-1, j] / rho_im
            rho_ic = max(U[1, i,   j], RHO_FLOOR); Y_ic = rhoY[i,   j] / rho_ic
            rho_ip = max(U[1, i+1, j], RHO_FLOOR); Y_ip = rhoY[i+1, j] / rho_ip
            rho_jm = max(U[1, i, j-1], RHO_FLOOR); Y_jm = rhoY[i, j-1] / rho_jm
            rho_jp = max(U[1, i, j+1], RHO_FLOOR); Y_jp = rhoY[i, j+1] / rho_jp
            YfL = mfL >= 0.0 ? Y_im : Y_ic
            YfR = mfR >= 0.0 ? Y_ic : Y_ip
            YfB = mfB >= 0.0 ? Y_jm : Y_ic
            YfT = mfT >= 0.0 ? Y_ic : Y_jp
            fluxY_r = (mfR * r_R * YfR - mfL * r_L * YfL) * dz
            fluxY_z = (mfT * YfT - mfB * YfB) * (r_c * dr)
            LY[i, j] = -(fluxY_r + fluxY_z) * invV
        end
    end
    return nothing
end

function rhs!(state::State, grid::Grid, topo::Topology;
              reflect_r_min, reflect_z_min)
    exchange_halos!(state.U, grid, topo, state.hbufs)
    exchange_halo_scalar!(state.rhoY, grid, topo, state.hbufs)
    apply_bcs!(state.U, grid, topo;
               reflect_r_min=reflect_r_min, reflect_z_min=reflect_z_min)
    apply_scalar_bcs!(state.rhoY, grid, topo)
    compute_fluxes!(state, grid)
    compute_residual!(state, grid)
    return nothing
end

# ---- CFL time-step --------------------------------------------------------

"""
    compute_dt(state, grid, topo, charge, cfl, dt_max) -> Float64

Standard 2-D CFL bound:  dt = cfl * min(dr, dz) / max_cells(|u| + c).
The local maximum is reduced across MPI ranks. Whenever any cell on any
rank still holds unburned/products explosive, the detonation velocity is
also factored in — the gamma-equivalent JWL sound speed under-estimates
the real wave speed badly enough that an air-only CFL would crash the
explicit step before the burn finishes.
"""
function compute_dt(state::State, grid::Grid, topo::Topology, charge::Charge,
                    cfl::Float64, dt_max::Float64)
    U   = state.U
    lam = state.lambda
    h   = min(grid.dr, grid.dz)
    smax = 0.0
    lam_max_local = 0.0
    @inbounds for j in interior_j(grid), i in interior_i(grid)
        rho, u_r, u_z, p = primitives((U[1,i,j], U[2,i,j], U[3,i,j], U[4,i,j]),
                                      lam[i, j])
        c = mix_sound_speed(rho, p, lam[i, j])
        s = sqrt(u_r*u_r + u_z*u_z) + c
        s > smax && (smax = s)
        lam[i, j] > lam_max_local && (lam_max_local = lam[i, j])
    end
    smax_global = MPI.Allreduce(smax, MPI.MAX, topo.comm)
    lam_max_global = MPI.Allreduce(lam_max_local, MPI.MAX, topo.comm)
    if lam_max_global > 0.0
        smax_global = max(smax_global, charge.det_vel)
    end
    smax_global = max(smax_global, 1.0e-30)
    return min(cfl * h / smax_global, dt_max)
end

# ---- time-stepping kernels -----------------------------------------------

"""
    step_forward_euler!(state, grid, topo, dt; bcs)

Single-stage 1st-order update. Debugging only.
"""
function step_forward_euler!(state::State, grid::Grid, topo::Topology, dt::Float64;
                             reflect_r_min::Bool, reflect_z_min::Bool)
    rhs!(state, grid, topo;
         reflect_r_min=reflect_r_min, reflect_z_min=reflect_z_min)
    @turbo for I in eachindex(state.U)
        state.U[I] += dt * state.L[I]
    end
    @turbo for I in eachindex(state.rhoY)
        state.rhoY[I] += dt * state.LrhoY[I]
    end
    return nothing
end

"""
    step_ssp_rk2!(state, grid, topo, dt; bcs)

Two-stage SSP-RK2 (Heun). Time accuracy 2 matches the MUSCL spatial order.
"""
function step_ssp_rk2!(state::State, grid::Grid, topo::Topology, dt::Float64;
                       reflect_r_min::Bool, reflect_z_min::Bool)
    copyto!(state.U0,    state.U)
    copyto!(state.rhoY0, state.rhoY)

    rhs!(state, grid, topo;
         reflect_r_min=reflect_r_min, reflect_z_min=reflect_z_min)
    @turbo for I in eachindex(state.U)
        state.U[I] = state.U0[I] + dt * state.L[I]
    end
    @turbo for I in eachindex(state.rhoY)
        state.rhoY[I] = state.rhoY0[I] + dt * state.LrhoY[I]
    end

    rhs!(state, grid, topo;
         reflect_r_min=reflect_r_min, reflect_z_min=reflect_z_min)
    @turbo for I in eachindex(state.U)
        state.U[I] = 0.5 * state.U0[I] + 0.5 * (state.U[I] + dt * state.L[I])
    end
    @turbo for I in eachindex(state.rhoY)
        state.rhoY[I] = 0.5 * state.rhoY0[I] + 0.5 * (state.rhoY[I] + dt * state.LrhoY[I])
    end
    return nothing
end
