# HLLC approximate Riemann solver for the compressible Euler equations.
#
# Reference: Toro, Spruce & Speares, "Restoration of the contact surface
# in the HLL Riemann solver", Shock Waves 4 (1994).
#
# The two-EOS extension (air on one side, JWL detonation products on the
# other) is handled by passing a separate burn fraction `lambdaL`,
# `lambdaR` in to each side; primitives and sound speeds use the local
# blend, which is the simplest physically-consistent way to apply HLLC
# at a contact-discontinuity-like burn-front interface.

@inline function _native_flux(U::NTuple{4,Float64}, p::Float64, dim::Int)
    rho = max(U[1], RHO_FLOOR)
    u_r = U[2] / rho
    u_z = U[3] / rho
    rhoE = U[4]
    if dim == 1
        return (rho*u_r,
                rho*u_r*u_r + p,
                rho*u_r*u_z,
                (rhoE + p) * u_r)
    else
        return (rho*u_z,
                rho*u_r*u_z,
                rho*u_z*u_z + p,
                (rhoE + p) * u_z)
    end
end

"""
    hllc_flux(UL, UR, lamL, lamR, dim) -> NTuple{4,Float64}

Numerical flux through a face whose outward normal is the `dim` direction
(`1` = r, `2` = z). UL, UR are the conservative states immediately to
the left and right of the face, and lamL, lamR are the local burn
fractions used to evaluate the EOS.
"""
@inline function hllc_flux(UL::NTuple{4,Float64}, UR::NTuple{4,Float64},
                   lamL::Float64, lamR::Float64, dim::Int)
    rhoL, urL, uzL, pL = primitives(UL, lamL)
    rhoR, urR, uzR, pR = primitives(UR, lamR)
    cL = mix_sound_speed(rhoL, pL, lamL)
    cR = mix_sound_speed(rhoR, pR, lamR)

    # Normal/tangential velocity components for this face.
    if dim == 1
        unL = urL; utL = uzL
        unR = urR; utR = uzR
    else
        unL = uzL; utL = urL
        unR = uzR; utR = urR
    end

    # Davis wave-speed estimates (Eq. 10.48, Toro 2009 3rd ed.).
    SL = min(unL - cL, unR - cR)
    SR = max(unL + cL, unR + cR)

    FL = _native_flux(UL, pL, dim)
    FR = _native_flux(UR, pR, dim)

    if SL >= 0.0
        return FL
    elseif SR <= 0.0
        return FR
    end

    # Contact-wave speed (Eq. 10.37, Toro).
    denom = rhoL * (SL - unL) - rhoR * (SR - unR)
    Sstar = (pR - pL + rhoL*unL*(SL - unL) - rhoR*unR*(SR - unR)) / denom

    EL = UL[4] / max(rhoL, RHO_FLOOR)
    ER = UR[4] / max(rhoR, RHO_FLOOR)

    if Sstar >= 0.0
        coef     = rhoL * (SL - unL) / (SL - Sstar)
        E_starL  = EL + (Sstar - unL) * (Sstar + pL / (rhoL * (SL - unL)))
        if dim == 1
            U_star = (coef, coef*Sstar, coef*utL, coef*E_starL)
        else
            U_star = (coef, coef*utL, coef*Sstar, coef*E_starL)
        end
        return ntuple(k -> FL[k] + SL * (U_star[k] - UL[k]), 4)
    else
        coef     = rhoR * (SR - unR) / (SR - Sstar)
        E_starR  = ER + (Sstar - unR) * (Sstar + pR / (rhoR * (SR - unR)))
        if dim == 1
            U_star = (coef, coef*Sstar, coef*utR, coef*E_starR)
        else
            U_star = (coef, coef*utR, coef*Sstar, coef*E_starR)
        end
        return ntuple(k -> FR[k] + SR * (U_star[k] - UR[k]), 4)
    end
end
