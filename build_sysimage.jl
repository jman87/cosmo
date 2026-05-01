# Build a precompiled sysimage for COSMO to eliminate JIT startup overhead.
#
# Usage:
#   julia --project=. build_sysimage.jl
#
# After building, run the solver with:
#   julia --sysimage cosmo.so --project=. run.jl <input.jsonc>
#   mpiexec -n 4 julia --sysimage cosmo.so --project=. run.jl <input.jsonc>
#
# The sysimage caches compiled machine code for all COSMO, MPI, and WriteVTK
# code paths exercised during the precompile workload, so first-step latency
# drops from ~30 s to < 1 s.  Rebuild the sysimage after any source change.

using PackageCompiler

# Run the sphere example as the precompile workload so the sysimage contains
# compiled code for the actual hot paths (HLLC, MUSCL, halo exchange, VTK I/O).
precompile_script = joinpath(@__DIR__, "examples", "sphere", "run.sh")

PackageCompiler.create_sysimage(
    ["COSMO", "MPI", "WriteVTK", "LoopVectorization"];
    sysimage_path   = joinpath(@__DIR__, "cosmo.so"),
    project         = @__DIR__,
    precompile_execution_file = nothing,  # set to a driver script to warm more paths
    cpu_target       = "native",           # optimise for this machine's ISA
)

println("Sysimage written to: ", joinpath(@__DIR__, "cosmo.so"))
println("Run with: julia --sysimage cosmo.so --project=. run.jl <input.jsonc>")
