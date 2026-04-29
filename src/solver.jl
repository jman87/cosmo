# Top-level driver: load config, build grid, set ICs, run the time loop.

"""
    run_case(input_path::AbstractString)

Run a single blast simulation defined by the JSON file at `input_path`.
The function expects MPI to be initialised already (see `cosmo/run.jl`).
On a single-rank invocation it falls back to a 1x1 process grid, so the
serial and parallel code paths are identical.
"""
function run_case(input_path::AbstractString)
    if !MPI.Initialized()
        MPI.Init()
    end
    topo = build_topology(MPI.COMM_WORLD)
    rank0 = topo.rank == 0

    cfg = load_config(input_path)
    grid = build_grid(cfg, topo)
    charge = build_charge(cfg)
    state = build_state(grid)

    set_initial_state!(state.U, state.lambda, state.rhoY, grid, charge)

    reflect_r_min = abs(cfg.r_min) < 1.0e-12
    reflect_z_min = abs(cfg.z_min) < 1.0e-12

    writer = build_writer(cfg, grid, topo)

    if rank0
        @printf "==========================================================\n"
        @printf "  cosmo  --  axisymmetric Eulerian FV blast solver\n"
        @printf "==========================================================\n"
        @printf "  ranks          : %d   process grid %dx%d\n" topo.nproc topo.dims[1] topo.dims[2]
        @printf "  global mesh    : %d x %d cells\n" grid.nr_global grid.nz_global
        @printf "  domain         : r in [%.3f, %.3f]  z in [%.3f, %.3f] in\n" cfg.r_min cfg.r_max cfg.z_min cfg.z_max
        @printf "  dr, dz         : %.4e, %.4e in\n" grid.dr grid.dz
        @printf "  charge         : %s   mass %.4e lbf*s^2/in (weight %.3f lbf)   radius %.4f in\n" charge.shape cfg.charge_mass (cfg.charge_mass * G_C) charge.radius
        @printf "  IC             : %s   D = %.3e in/s\n" charge.ic_kind charge.det_vel
        @printf "  EOS            : air (ideal, gamma=1.4) + JWL (TNT)\n"
        @printf "  ambient air    : rho=%.3e  p=%.3f psi\n" RHO_AIR_STP P_ATM
        @printf "  scheme         : %s   CFL=%.2f   t_end=%.3e s\n" cfg.scheme cfg.cfl cfg.t_end
        @printf "  BC: r_min %s,   z_min %s\n" (reflect_r_min ? "reflect" : "outflow") (reflect_z_min ? "reflect" : "outflow")
        @printf "----------------------------------------------------------\n"
    end

    update_burn!(state.lambda, grid, charge, 0.0)
    exchange_halo_scalar!(state.lambda, grid, topo)
    write_snapshot!(writer, state, 0.0, 0)

    stepper! = cfg.scheme === :ssp_rk2 ? step_ssp_rk2! : step_forward_euler!

    t = 0.0
    step = 0
    wall_t0 = time()

    while t < cfg.t_end
        # Refresh programmed-burn front (no-op for BRODE).
        update_burn!(state.lambda, grid, charge, t)
        exchange_halo_scalar!(state.lambda, grid, topo)

        dt = compute_dt(state, grid, topo, charge, cfg.cfl, cfg.dt_max)
        dt = min(dt, cfg.t_end - t)
        if !(dt > 0.0) || !isfinite(dt)
            rank0 && @warn "Non-positive or non-finite dt; stopping" t dt
            break
        end

        stepper!(state, grid, topo, dt;
            reflect_r_min=reflect_r_min, reflect_z_min=reflect_z_min)
        t += dt
        step += 1

        if step % cfg.out_freq == 0
            write_snapshot!(writer, state, t, step)
            if rank0
                wall = time() - wall_t0
                @printf "  step %6d  t = %.4e s  dt = %.4e s  wall = %.1f s\n" step t dt wall
            end
        end
    end

    write_snapshot!(writer, state, t, step)
    close_writer!(writer)

    if rank0
        wall = time() - wall_t0
        @printf "----------------------------------------------------------\n"
        @printf "  finished: %d steps, t = %.4e s, wall = %.1f s\n" step t wall
        @printf "  results : %s.pvd\n" writer.pvd_path
        @printf "==========================================================\n"
    end

    MPI.Barrier(topo.comm)
    return nothing
end
