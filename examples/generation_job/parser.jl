using ArgParse



function parser(args; efield=false, kw...)
    s = ArgParseSettings(; kw...)
    if efield
        @add_arg_table! s begin
            "incident_angle"
                help = "incident angle of light from normal in degrees"
                arg_type = Float64
                required = true
            "beam_sigma"
                help = "standard deviation of the gaussian beam; 0 means uniform"
                arg_type = Float64
                required = true
        end
    end
    @add_arg_table! s begin
        "radius"
            help = "mirror radius"
            arg_type = Float64
            required = true
        "nrings"
            help = "number of rings"
            arg_type = UInt32
            required = true
        "roughness_rms"
            help = "RMS height of surface roughness"
            arg_type = Float64
            default = 0.0
        "roughness_sigma"
            help = "standard deviation of frequency of mirror roughness"
            arg_type = Float64
            default = 0.0
    end
    return parse_args(args, s)
end



function filename(prefix, argdict, efield)
    r     = argdict["radius"]
    n     = argdict["nrings"]
    rms   = argdict["roughness_rms"]
    sigma = argdict["roughness_sigma"]
    name = "$prefix-r$r-n$n-rms$rms-s$sigma"
    if efield
        angle = argdict["incident_angle"]
        beamsigma = argdict["beam_sigma"]
        name *= "-a$angle-bs$beamsigma"
    end
    return name
end
