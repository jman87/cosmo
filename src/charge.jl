# Charge geometry, initial conditions, and programmed-burn rule.
#
# Two initial-condition modes are available:
#
# 1. PROGRAMMED_BURN  (default; matches the original Python solver)
#    * Inside the charge: unreacted explosive at rest, with rho = rho0_TNT
#      and specific internal energy e = e_CJ. Burn fraction lambda(r,z,t)
#      starts at 0 and reaches 1 when the spherical detonation front
#      (radius D*t from the initiation point) sweeps past the cell.
#    * Outside the charge: ambient air at STP.
#    * The mixed pressure mix_pressure(rho, e, lambda) gradually transitions
#      from ideal-gas behaviour (lambda = 0) to JWL behaviour (lambda = 1)
#      as the detonation front reaches each cell.
#
# 2. BRODE_BALLOON
#    * Inside the charge: high-internal-energy gas (e = Q_TNT) at rho0_TNT,
#      lambda = 1 from t = 0 (full JWL EOS).
#    * Outside the charge: ambient air at STP, lambda = 0.
#    * No subsequent burn update; lambda is frozen.
#    * Cheap and unconditionally robust, useful for very coarse meshes.

struct Charge
    shape::Symbol                # :sphere or :cylinder
    r_c::Float64                 # centre, in
    z_c::Float64
    radius::Float64              # in
    half_height::Float64         # cylinder only; equals radius for sphere
    rho0::Float64                # unreacted density, lbf*s^2/in^4
    ecj::Float64                 # CJ specific internal energy, in^2/s^2
    qtnt::Float64                # heat of detonation (used by BRODE), in^2/s^2
    det_vel::Float64             # detonation velocity, in/s
    ic_kind::InitialConditionKind
end

function build_charge(cfg::CaseConfig)
    shape  = Symbol(cfg.charge_shape)
    mass   = cfg.charge_mass                    # lbf*s^2/in (consistent units)
    volume = mass / TNT.rho0                    # in^3

    radius, half_h = if shape === :sphere
        r = (3.0 * volume / (4.0 * pi))^(1.0/3.0)
        (r, r)
    elseif shape === :cylinder
        h = cfg.charge_half_height
        h > 0 || error("cylinder_half_height must be positive")
        r = sqrt(volume / (2.0 * pi * h))
        (r, h)
    else
        error("Unknown charge shape: $(shape)")
    end

    return Charge(shape, cfg.charge_center_r, cfg.charge_center_z,
                  radius, half_h, TNT.rho0, TNT.ecj,
                  4.184e6 * J_KG_TO_IN2_S2,
                  cfg.detonation_velocity, cfg.ic_kind)
end

@inline function inside_charge(c::Charge, r::Real, z::Real)
    dr = r - c.r_c
    dz = z - c.z_c
    if c.shape === :sphere
        return (dr*dr + dz*dz) <= c.radius*c.radius
    else
        return (abs(dz) <= c.half_height && dr*dr <= c.radius*c.radius)
    end
end

# Distance from a cell point to the charge centre. The programmed-burn
# rule treats the detonation as a spherical wave from the centre, so the
# arrival time of the burn front at (r, z) is dist / D regardless of the
# charge shape (this is the same simplification the Python solver used).
@inline function dist_to_center(c::Charge, r::Real, z::Real)
    dr = r - c.r_c
    dz = z - c.z_c
    return sqrt(dr*dr + dz*dz)
end

"""
    set_initial_state!(U, lambda, grid, charge)

Fill the conservative-state array `U` and the burn-fraction array
`lambda` with the t = 0 condition selected by `charge.ic_kind`.
A small Cartesian sub-grid is sampled inside each cell to compute the
charge volume fraction, which is used to blend the in-charge and ambient
states across cells that the geometric boundary cuts through.
"""
function set_initial_state!(U::Array{Float64,3}, lambda::Array{Float64,2},
                            rhoY::Array{Float64,2},
                            grid::Grid, charge::Charge)
    NSUB = 5
    inv_n = 1.0 / NSUB

    # Air at STP (u = 0 so total energy = internal energy)
    rhoE_air = RHO_AIR_STP * E_AIR_STP

    # In-charge specific internal energy depends on the IC mode.
    e_in = charge.ic_kind === BRODE_BALLOON ? charge.qtnt : charge.ecj
    rhoE_in = charge.rho0 * e_in

    @inbounds for j in 1:grid.nzj
        zc = grid.z_centers[j]
        for i in 1:grid.nri
            rc = grid.r_centers[i]

            # Sub-cell volume fraction inside the charge.
            phi = 0.0
            for sj in 0:(NSUB-1), si in 0:(NSUB-1)
                rs = rc + ((si + 0.5) * inv_n - 0.5) * grid.dr
                zs = zc + ((sj + 0.5) * inv_n - 0.5) * grid.dz
                phi += inside_charge(charge, rs, zs) ? 1.0 : 0.0
            end
            phi *= inv_n * inv_n

            rho_cell  = phi * charge.rho0 + (1.0 - phi) * RHO_AIR_STP
            rhoE_cell = phi * rhoE_in     + (1.0 - phi) * rhoE_air

            U[1, i, j] = rho_cell
            U[2, i, j] = 0.0
            U[3, i, j] = 0.0
            U[4, i, j] = rhoE_cell

            # Burn fraction at t = 0:
            #   PROGRAMMED_BURN: lambda = 0 until the detonation front arrives
            #   BRODE_BALLOON  : lambda = phi (immediate JWL state inside)
            lambda[i, j] = charge.ic_kind === BRODE_BALLOON ? phi : 0.0

            # Material mass fraction Y: 1 inside the charge geometry, 0 in
            # ambient air. rho*Y is the conservative form that the
            # passive-scalar advection update evolves with the flow.
            rhoY[i, j] = rho_cell * phi
        end
    end
    return U, lambda, rhoY
end

"""
    update_burn!(lambda, grid, charge, t)

Refresh the programmed-burn fraction at the current time. A cell at
distance d from the charge centre is fully burned (lambda = 1) when
D * t >= d, otherwise lambda = 0. A binary rule is the simplest robust
choice and matches the Python solver. Only used when the IC mode is
PROGRAMMED_BURN; for BRODE_BALLOON the lambda field is frozen at t = 0.
"""
function update_burn!(lambda::Array{Float64,2}, grid::Grid, charge::Charge, t::Float64)
    if charge.ic_kind !== PROGRAMMED_BURN
        return lambda
    end
    Dt = charge.det_vel * t
    @inbounds for j in 1:grid.nzj, i in 1:grid.nri
        rc = grid.r_centers[i]
        zc = grid.z_centers[j]
        if inside_charge(charge, rc, zc) && dist_to_center(charge, rc, zc) <= Dt
            lambda[i, j] = 1.0
        end
    end
    return lambda
end
