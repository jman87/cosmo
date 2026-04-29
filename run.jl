#!/usr/bin/env julia
# Entry-point script for cosmo.
#
# Usage:
#   serial:    julia --project=cosmo cosmo/run.jl path/to/input.json
#   parallel:  mpiexec -n 4 julia --project=cosmo cosmo/run.jl path/to/input.json
#
# `MPI.Init` is called on entry; serial runs work transparently because
# MPI sees a single-rank world and the geometric partitioner picks a
# trivial 1x1 process grid.

using Pkg
Pkg.activate(@__DIR__)

include(joinpath(@__DIR__, "src", "COSMO.jl"))
using .COSMO
using MPI

if length(ARGS) != 1
    if MPI.Initialized() ? MPI.Comm_rank(MPI.COMM_WORLD) == 0 : true
        println(stderr, "Usage: julia --project=cosmo cosmo/run.jl <input.json>")
    end
    exit(1)
end

MPI.Initialized() || MPI.Init()
try
    COSMO.run_case(ARGS[1])
finally
    MPI.Finalized() || MPI.Finalize()
end
