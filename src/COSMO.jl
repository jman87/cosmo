"""
    COSMO

Two-dimensional axisymmetric Eulerian finite-volume blast solver.

The conservative state U = (rho, rho*u_r, rho*u_z, rho*E) is advanced on a
structured (r, z) grid using:

* HLLC approximate Riemann solver  (Toro, Spruce & Speares, 1994)
* MUSCL spatial reconstruction with the minmod slope limiter (TVD)
* SSP-RK2 (Heun) explicit time integration  (Shu & Osher, 1988)
* Two built-in EOS: ideal gas (gamma = 1.4) for air, JWL for TNT
* Programmed-burn IC by default; Brode compressed-gas balloon optional
* MPI parallelism with geometric block decomposition + halo exchange

Imperial units throughout (matching the original Python solver):
    length   in
    time     s
    mass     lbf*s^2/in
    pressure psi          (lbf/in^2)
    density  lbf*s^2/in^4
    energy   in^2/s^2     (specific internal energy)
"""
module COSMO

using JSON
using Printf
using WriteVTK
using MPI

include("eos.jl")
include("config.jl")
include("parallel.jl")
include("grid.jl")
include("parallel_halo.jl")
include("charge.jl")
include("reconstruct.jl")
include("riemann.jl")
include("boundary.jl")
include("timeint.jl")
include("output.jl")
include("solver.jl")

export run_case

end # module COSMO
