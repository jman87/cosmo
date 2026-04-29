# MUSCL reconstruction with a minmod slope limiter.
#
# For each cell we form a limited slope from neighbouring cell averages and
# extrapolate a left/right state at every face. minmod is the most
# diffusive of the standard TVD limiters but it is also the most robust:
# it never produces overshoots, never amplifies negative-density excursions
# near very strong shocks, and is the limiter of choice when correctness
# is more important than crispness. (Other choices: van Leer, MC, superbee.)

@inline function minmod(a::Float64, b::Float64)
    if a*b <= 0.0
        return 0.0
    elseif abs(a) < abs(b)
        return a
    else
        return b
    end
end

"""
    limited_slope(uL, uC, uR) -> Float64

Minmod-limited centred slope between the left, centre and right cell
averages. The slope is per *cell*, so the half-cell extrapolation to a
face uses `±0.5 * slope`.
"""
@inline limited_slope(uL::Float64, uC::Float64, uR::Float64) =
    minmod(uC - uL, uR - uC)

"""
    reconstruct_face(U_LL, U_L, U_R, U_RR) -> (UL_face, UR_face)

Build the left- and right-of-face states between cell L and cell R using
a 4-point MUSCL stencil. Each component is reconstructed independently;
limiting on conservative variables is the simplest robust choice for a
fixed-grid Eulerian solver.
"""
@inline function reconstruct_face(
    U_LL::NTuple{4,Float64},
    U_L ::NTuple{4,Float64},
    U_R ::NTuple{4,Float64},
    U_RR::NTuple{4,Float64},
)
    UL_face = ntuple(k -> U_L[k] + 0.5 * limited_slope(U_LL[k], U_L[k], U_R[k]),  4)
    UR_face = ntuple(k -> U_R[k] - 0.5 * limited_slope(U_L[k],  U_R[k], U_RR[k]), 4)
    return UL_face, UR_face
end
