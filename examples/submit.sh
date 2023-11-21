#!/bin/bash

dir="$(dirname "$0")"

module purge
module load julia parallel slurm-auto-array

mirror_params='::: 3 10 30 ::: 100 ::: 0.01 0.1 0.5 ::: 1 3 10'
efield_params='::: 0 15 30 45 75 ::: 0 1 3 10'



submit_array_job() {
    script="$1"
    params="$2"
    shift; shift
    bash -c "parallel echo $params" | slurm-auto-array --qos standby \
                                                       --requeue \
                                                       --parsable \
                                                       -J "$script" \
                                                       -o "$script-%A-%a.out" \
                                                       "$@" -- \
                                                       julia --project "$dir/.." --threads auto "$dir/$script"
}



mirror_job_id="$(     submit_array_job generate_mirror      "$mirror_params" \
                                       -t 1  --mem 1G)"

impedance_job_id="$(  submit_array_job generate_impedance   "$mirror_params" \
                                       -t 60 --mem 250G -N 1 -n 100 -d afterok:$mirror_job_id)"

efield_job_id="$(     submit_array_job generate_efield      "$mirror_params $efield_params" \
                                       -t 1  --mem 1G               -d afterok:$mirror_job_id)"

current_job_id="$(    submit_array_job generate_current     "$mirror_params $efield_params" \
                                       -t 60 --mem 500G -N 1 -n 120 -d afterok:$impedance_job_id:$efield_job_id)"

reflectance_job_id="$(submit_array_job generate_reflectance "$mirror_params $efield_params" \
                                       -t 10 --mem 10G -d afterok:$current_job_id)"



echo "Mirror job ID:         $efield_job_id"
echo "Impedance job ID:      $impedance_job_id"
echo "Electric field job ID: $efield_job_id"
echo "Current job ID:        $current_job_id"
echo "Reflectance job ID:    $reflectance_job_id"
