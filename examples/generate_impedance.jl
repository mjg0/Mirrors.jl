#!/usr/bin/env julia

#SBATCH -t 60 --mem 250G --n 100 -N 1
#SBATCH -q standby --requeue

# Generate and write an impedance matrix file given a mirror file.

using Mirrors, ArgParse



"""
    parse(args)

Parse command line arguments and return `infile` and `outfile`.
"""
function parse(args)
    s = ArgParseSettings(description="Read a `Mirror` from `infile`, generate the
                                      corresponding impedance matrix, and write it to
                                      `outfile`.",
                         epilog="Checkpointing occurs automatically--given the same output
                                 file name, it can be interrupted and restarted an arbitrary
                                 amount of times and still successfully generate `outfile`.
                                 Checkpoint files are cleaned up automatically")
    @add_arg_table s begin
        "infile"
            help = "an existing mirror file"
            required = true
        "outfile"
            help = "the name of the impedance file to be written"
            required = true
    end
    parsed = parse_args(args, s)
    return parsed["infile"], parsed["outfile"]
end



"""
    writeimpedance(infile, outfile)

Generate and write an impedance to `outfile` given a `Mirror` stored in `infile`.

Checkpointing and resumption happens automatically.
"""
function writeimpedance(infile, outfile)
    # Read in mirror
    mirror = Mirror(infile)
    n = 4 * length(mirror)
    # Short circuit if the output file has already been written
    if filesize(outfile) == n^2*sizeof(ComplexF64) return end
    # Name of checkpoint file given column number
    checkfilename(j) = "outfile.checkpoint-$(lpad(j, ndigits(length(mirror)), '0'))-of-$n"
    # Fill Z, checkpointing (and skipping columns that have been checkpointed) as we go
    Z = Matrix{ComplexF64}(undef, n, n)
    @inbounds @fastmath Threads.@threads :dynamic for (j, p2) in collect(enumerate(mirror))
        # Short-circuit if the checkpoint file has already been written
        if filesize(checkfilename(j)) != n*4*sizeof(ComplexF64)
            for (i, p1) in enumerate(mirror)
                Zblock = view(Z, 4i-3:4i, 4j-3:4j)
                if singular && i==j # main diagonal
                    singularblockfill!(f, Zblock, mirror, p1)
                else                # off-diagonal
                    nonsingularblockfill!(f, Zblock, mirror, p1, p2)
                end
            end
            # Write checkpoint file
            write(checkfilename(j), Z[:,4j-3:4j])
            println("Wrote checkpoint $(checkfilename(j))")
        end
    end
    # Combine checkpoint files into a single output file
    open(outfile, "w") do f
        for j in 1:length(mirror)
            write(f, read(checkfilename(j)))
        end
    end
    println("Wrote output file $outfile")
    # Clean up checkpoint files
    rm.(checkfilename.(1:length(mirror)))
end



function main(args=ARGS)
    # Parse
    infile, outfile = parse(args)
    # Alert if this job has been requeued
    # TODO
    # Fill and write impedance
    writeimpedance(infile, outfile)
end



if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
