# JSON input-file loader.
#
# Layout matches the original Python solver. New optional knobs:
#   * charge.detonation_velocity (in/s)  -- defaults to TNT CJ value
#   * charge.explosive_material  ("tnt") -- only TNT for now
#   * charge.initial_condition   ("programmed_burn" | "brode") -- default
#                                   programmed_burn now that we have JWL
#
# Unknown keys are ignored (e.g. legacy shock_capturing / solver_method).

@enum InitialConditionKind PROGRAMMED_BURN BRODE_BALLOON

struct CaseConfig
    nr::Int
    nz::Int
    r_min::Float64
    r_max::Float64
    z_min::Float64
    z_max::Float64

    charge_shape::String                 # "sphere" or "cylinder"
    charge_mass::Float64                 # mass of charge in lbf*s^2/in (consistent units throughout)
    charge_center_r::Float64
    charge_center_z::Float64
    charge_half_height::Float64          # cylinder only
    charge_material::String              # "tnt"
    detonation_velocity::Float64         # in/s
    ic_kind::InitialConditionKind

    t_end::Float64
    cfl::Float64
    dt_max::Float64
    scheme::Symbol                        # :forward_euler | :ssp_rk2

    out_dir::String
    out_freq::Int
    pvd_name::String
end

"""
    load_config(path) -> CaseConfig

Parse a JSON input file. Optional keys fall back to sensible defaults.
"""
# Read a JSON file that may contain `//` line comments and `/* ... */`
# block comments (neither of which the strict JSON spec allows). String
# literals are preserved untouched; only comments outside string literals
# are stripped before the result is handed to JSON.parse.
function read_json_with_comments(path::AbstractString)
    raw = read(path, String)
    # Strip /* ... */ block comments first (non-greedy, multi-line).
    cleaned = replace(raw, r"/\*[\s\S]*?\*/" => "")
    # Strip // line comments while preserving strings. A simple regex that
    # matches // outside of double-quoted strings:
    cleaned = replace(cleaned, r"(\"(?:[^\"\\]|\\.)*\")|//[^\n]*" => s"\1")
    return JSON.parse(cleaned)
end

function load_config(path::AbstractString)
    isfile(path) || error("Config file not found: $path")
    raw = read_json_with_comments(path)

    mesh = _required(raw, "mesh")
    domain = _required(raw, "domain")
    charge = _required(raw, "charge")
    time_ = _required(raw, "time")
    out = _required(raw, "output")

    shape = lowercase(string(_required(charge, "shape")))
    shape in ("sphere", "cylinder") || error("Unknown charge shape: $shape")

    half_h = if shape == "cylinder"
        Float64(_required(charge, "cylinder_half_height"))
    else
        0.0
    end

    material = lowercase(string(get(charge, "explosive_material", "tnt")))
    material == "tnt" || error("Unsupported explosive_material: $material (only 'tnt' is built in)")

    det_vel = Float64(get(charge, "detonation_velocity", TNT.dcj))

    ic_str = lowercase(string(get(charge, "initial_condition", "programmed_burn")))
    ic_kind = if ic_str in ("programmed_burn", "burn")
        PROGRAMMED_BURN
    elseif ic_str in ("brode", "balloon")
        BRODE_BALLOON
    else
        error("Unknown initial_condition: $ic_str (use 'programmed_burn' or 'brode')")
    end

    center = _required(charge, "center")

    scheme_str = lowercase(string(get(time_, "scheme", "ssp_rk2")))
    scheme = if scheme_str in ("ssp_rk2", "rk2", "heun")
        :ssp_rk2
    elseif scheme_str in ("forward_euler", "fe", "euler")
        :forward_euler
    else
        error("Unsupported time scheme: $scheme_str (use 'ssp_rk2' or 'forward_euler')")
    end

    # Charge size: accept either "mass" (lbf*s^2/in, the consistent-unit
    # value used throughout the solver) or "weight" (lbf, converted via
    # m = W / g_c). Exactly one must be supplied. The legacy "weight_lbm"
    # key is rejected with a clear message.
    has_mass   = haskey(charge, "mass")
    has_weight = haskey(charge, "weight")
    if haskey(charge, "weight_lbm")
        error("'weight_lbm' is no longer accepted; use 'mass' (lbf*s^2/in) or 'weight' (lbf)")
    end
    (has_mass || has_weight) ||
        error("Charge size must be specified as either 'mass' (lbf*s^2/in) or 'weight' (lbf)")
    !(has_mass && has_weight) ||
        error("Charge size: specify either 'mass' or 'weight', not both")
    charge_mass = has_mass ?
        Float64(charge["mass"]) :
        Float64(charge["weight"]) / G_C

    return CaseConfig(
        Int(_required(mesh, "nr")),
        Int(_required(mesh, "nz")),
        Float64(get(domain, "r_min", 0.0)),
        Float64(_required(domain, "r_max")),
        Float64(get(domain, "z_min", 0.0)),
        Float64(_required(domain, "z_max")), shape,
        charge_mass,
        Float64(_required(center, "r")),
        Float64(_required(center, "z")),
        half_h,
        material,
        det_vel,
        ic_kind, Float64(_required(time_, "t_end")),
        Float64(get(time_, "cfl", 0.4)),
        Float64(get(time_, "dt_max", Inf)),
        scheme, String(get(out, "directory", "output")),
        Int(get(out, "frequency", 10)),
        String(get(out, "pvd_name", "cosmo")),
    )
end

@inline function _required(d::AbstractDict, key::AbstractString)
    haskey(d, key) || error("Missing required config field: '$key'")
    return d[key]
end
