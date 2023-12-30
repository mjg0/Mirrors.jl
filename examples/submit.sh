#!/bin/bash

dir="$(dirname "$0")"

module purge
module load julia parallel slurm-auto-array
mkdir -p job-logs

#                  radius      nrings  roughness rms      roughness sigma
mirror_params='::: 3 10 30 ::: 100 ::: 0 0.01 0.1 ::: 1 3 10'
#                  angle             beam sigma
efield_params='::: 0 15 30 45 75 ::: 0 1 3 10'



submit_array_job() {
    local script="$1"
    local params="$2"
    local t="$3"
    local n="$4"
    local mem="$5"
    shift 5

    local threads=auto
    grep -qE 'current$' <<< "$script" && threads=1

    bash -c "parallel echo $params" | slurm-auto-array --time "$t" \
                                                       --ntasks "$n" --nodes 1 \
                                                       --mem "$mem" \
                                                       -J "${script#generate_}" \
                                                       -l "job-logs/$script-%A.log" \
                                                       -o "job-logs/$script-%A-%a.out" \
                                                       --parsable \
                                                       "$@" -- \
                                                 julia --project="$dir/.." \
                                                       --threads="$threads" \
                                                       --heap-size-hint="$mem" \
                                                       "$dir/$script"
}



mirror_job_id="$(submit_array_job generate_mirror "$mirror_params" 5 1 2G)"

impedance_job_id="$(submit_array_job generate_impedance "$mirror_params" 90 40 350G -d afterok:$mirror_job_id)"

efield_job_id="$(submit_array_job generate_efield "$efield_params $mirror_params" 5 1 1G -d afterok:$mirror_job_id)"

current_job_id="$(submit_array_job generate_current "$efield_params $mirror_params" 90 128 500G \
                  -d afterok:$impedance_job_id:$efield_job_id)"

reflectance_job_id="$(submit_array_job generate_reflectance "$efield_params $mirror_params" 30 1 1G \
                      -d afterok:$current_job_id)"

plot_job_id="$(submit_array_job generate_plots "$efield_params $mirror_params" 5 1 4G -d afterok:$reflectance_job_id)"



echo "Mirror job ID:         $efield_job_id"
echo "Impedance job ID:      $impedance_job_id"
echo "Electric field job ID: $efield_job_id"
echo "Current job ID:        $current_job_id"
echo "Reflectance job ID:    $reflectance_job_id"
echo "Plotting job ID:       $plot_job_id"