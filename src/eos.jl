# Equation of state and physical constants.
#
# Two EOS are supported:
#
#   1. IdealGas (gamma-law) for ambient air.
#         p = (gamma - 1) * rho * e
#
#   2. JWL (Jones-Wilkins-Lee) for detonation products.
#         p = A (1 - omega/(R1 eta)) exp(-R1/eta)
#           + B (1 - omega/(R2 eta)) exp(-R2/eta)
#           + omega rho e
#         eta = rho / rho0
#
# A burn fraction lambda in [0, 1] blends the two:
#
#         p_mix = (1 - lambda) * p_air + lambda * p_jwl
#
# lambda is set by the programmed-burn rule (see charge.jl): lambda = 1
# everywhere the detonation front has reached, lambda = 0 elsewhere.
#
# Imperial units throughout:
#       length   in
#       time     s
#       mass     lbf*s^2/in
#       pressure psi          (lbf/in^2)
#       density  lbf*s^2/in^4
#       energy   in^2/s^2

# --- Universal constants and unit conversions ----------------------------
const PA_TO_PSI           = 1.450377e-4
const KG_M3_TO_LBF_S2_IN4 = 0.001940320 / 20736          # ≈ 9.357e-8
const J_KG_TO_IN2_S2      = 1550.003                      # 1 J/kg = 1550.003 in^2/s^2
const G_C                 = 386.088                       # in/s^2 (weight in lbf -> mass in lbf*s^2/in via m = W / g_c)

# --- Ideal-gas EOS --------------------------------------------------------
struct IdealGas
    gamma::Float64
end

const AIR = IdealGas(1.4)

# --- JWL EOS --------------------------------------------------------------
struct JWL
    A::Float64        # psi
    B::Float64        # psi
    R1::Float64
    R2::Float64
    omega::Float64
    rho0::Float64     # lbf*s^2/in^4
    pcj::Float64      # psi  (informational)
    ecj::Float64      # in^2/s^2  (CJ specific internal energy)
    dcj::Float64      # in/s    (CJ detonation velocity)
end

# TNT JWL parameters (Dobratz & Crawford 1985 / LLNL EOS handbook),
# converted from SI to imperial.
const TNT = let
    A_si    = 373.77e9    # Pa
    B_si    =   3.747e9   # Pa
    rho0_si = 1630.0      # kg/m^3
    pcj_si  =  21.0e9     # Pa
    ecj_si  =   6.0e6     # J/kg  (specific energy at the CJ state)
    dcj_si  =   6930.0    # m/s   (TNT detonation velocity)
    JWL(
        A_si    * PA_TO_PSI,
        B_si    * PA_TO_PSI,
        4.15,
        0.90,
        0.35,
        rho0_si * KG_M3_TO_LBF_S2_IN4,
        pcj_si  * PA_TO_PSI,
        ecj_si  * J_KG_TO_IN2_S2,
        dcj_si  * 39.3700787,                # m/s -> in/s
    )
end

# --- Ambient air (standard atmosphere) -----------------------------------
const RHO_AIR_STP = 1.225 * KG_M3_TO_LBF_S2_IN4               # lbf*s^2/in^4
const P_ATM       = 14.696                                    # psi
const E_AIR_STP   = P_ATM / ((AIR.gamma - 1.0) * RHO_AIR_STP) # in^2/s^2

# --- Numerical floors ----------------------------------------------------
const RHO_FLOOR = 1.0e-12
const P_FLOOR   = 1.0e-6
const E_FLOOR   = 1.0e-6

# --- EOS evaluations ------------------------------------------------------
@inline pressure_ideal(eos::IdealGas, rho, e) = (eos.gamma - 1.0) * rho * e

@inline function pressure_jwl(eos::JWL, rho, e)
    eta = rho / eos.rho0
    inv_eta = 1.0 / eta
    t1 = eos.A * (1.0 - eos.omega / (eos.R1 * eta)) * exp(-eos.R1 * inv_eta)
    t2 = eos.B * (1.0 - eos.omega / (eos.R2 * eta)) * exp(-eos.R2 * inv_eta)
    t3 = eos.omega * rho * e
    return t1 + t2 + t3
end

# Approximate squared sound speed for each EOS, used only by the CFL
# estimator. For JWL we use the gamma-equivalent c^2 = (omega + 1) * p / rho
# which bounds the true value from above for the parameter ranges used here
# and is the standard CFL estimate in production blast codes.
@inline sound_speed_ideal(eos::IdealGas, rho, p) =
    sqrt(eos.gamma * max(p, P_FLOOR) / max(rho, RHO_FLOOR))
@inline sound_speed_jwl(eos::JWL, rho, p) =
    sqrt((eos.omega + 1.0) * max(p, P_FLOOR) / max(rho, RHO_FLOOR))

"""
    mix_pressure(rho, e, lambda) -> p

Reaction-progress weighted pressure. lambda in [0, 1] blends air (gamma
law) with TNT detonation products (JWL):

    p = (1 - lambda) * p_air + lambda * p_jwl
"""
@inline function mix_pressure(rho, e, lambda)
    rho_s = max(rho, RHO_FLOOR)
    e_s   = max(e,   E_FLOOR)
    p_air = pressure_ideal(AIR, rho_s, e_s)
    if lambda <= 0.0
        return max(p_air, P_FLOOR)
    end
    p_jwl = pressure_jwl(TNT, rho_s, e_s)
    return max((1.0 - lambda) * p_air + lambda * p_jwl, P_FLOOR)
end

"""
    mix_sound_speed(rho, p, lambda) -> c

Same blend as `mix_pressure`, applied to the per-EOS sound-speed estimates.
Used only by the CFL bound, so an inexpensive over-estimate is fine.
"""
@inline function mix_sound_speed(rho, p, lambda)
    c_air = sound_speed_ideal(AIR, rho, p)
    if lambda <= 0.0
        return c_air
    end
    c_jwl = sound_speed_jwl(TNT, rho, p)
    return (1.0 - lambda) * c_air + lambda * c_jwl
end

"""
    primitives(U, lambda) -> (rho, u_r, u_z, p)

Recover primitive variables from `U = (rho, rho*u_r, rho*u_z, rho*E)` using
the EOS blend selected by `lambda` (0 = air, 1 = TNT products).
"""
@inline function primitives(U::NTuple{4,Float64}, lambda::Float64)
    rho = max(U[1], RHO_FLOOR)
    u_r = U[2] / rho
    u_z = U[3] / rho
    e   = U[4] / rho - 0.5 * (u_r*u_r + u_z*u_z)
    p   = mix_pressure(rho, e, lambda)
    return rho, u_r, u_z, p
end

"""
    conservative(rho, u_r, u_z, p, lambda) -> NTuple{4,Float64}

Build the conservative state vector from primitives. Inverts the same
EOS blend as `primitives`. For lambda == 0 this is exact; for lambda > 0
we invert the blended pressure by solving for e via Newton from the
ideal-gas guess, which converges in 2-4 iterations.
"""
@inline function conservative(rho, u_r, u_z, p, lambda)
    e = if lambda <= 0.0
        p / ((AIR.gamma - 1.0) * rho)
    else
        e_guess = p / ((AIR.gamma - 1.0) * rho)
        _invert_mix(rho, p, lambda, e_guess)
    end
    rhoE = rho * (e + 0.5 * (u_r*u_r + u_z*u_z))
    return (rho, rho*u_r, rho*u_z, rhoE)
end

@inline function _invert_mix(rho, p_target, lambda, e0)
    e = e0
    for _ in 1:8
        p_now  = mix_pressure(rho, e, lambda)
        # d p / d e at fixed rho: ideal contributes (gamma-1)*rho;
        # JWL contributes omega*rho. Both are independent of e.
        dp_de  = (1.0 - lambda) * (AIR.gamma - 1.0) * rho +
                 lambda * TNT.omega * rho
        de     = (p_target - p_now) / max(dp_de, 1.0e-30)
        e += de
        abs(de) < 1.0e-10 * max(abs(e), 1.0e-30) && break
    end
    return e
end
