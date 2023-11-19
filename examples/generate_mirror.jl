#!/usr/bin/env julia

#SBATCH -t 1 --mem 1G
#SBATCH -q standby --requeue

# Generate a mirror and write it to a file

using Mirrors, ArgParse, Printf



"""
    parse(args)

Parse command line arguments and return `r`, `n`, `rms`, and `sigma`.
"""
function parse(args)
    s = ArgParseSettings(description="Generate a mirror and write it to a file named
                                      `mirror-\$r-\$n-\$rms-\$sigma`.")
    @add_arg_table s begin
        "r"
            help = "mirror radius"
            arg_type = Float64
            required = true
        "n"
            help = "number of rings"
            arg_type = UInt32
            required = true
        "rms"
            help = "RMS height of surface roughness"
            arg_type = Float64
            default = 0.0
        "sigma"
            help = "standard deviation of frequency of mirror roughness"
            arg_type = Float64
            default = 0.0
    end
    parsed = parse_args(args, s)
    return parsed["r"], parsed["n"], parsed["rms"], parsed["sigma"]
end



function main(args=ARGS)
    # Parse
    r, n, rms, sigma = parse(args)
    # Alert if this job has been requeued
    # TODO
    # Create and write mirror
    mirror = Mirror(r, n, rms, sigma)
    filename = @sprintf "mirror-%09.4f-%04i-%07.4f-%09.4f" r n rms sigma
    write(mirror, filename)
end



if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
