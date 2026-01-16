#!/bin/bash

# parameter values
depvals=(1 4)
typevals=(2 3 1)
theta_vals=(1 2)
scaled_vals=(1 4)

mkdir -p logstab2tabE2

for type in "${typevals[@]}"; do
  for dep in "${depvals[@]}"; do
    for varz in "${theta_vals[@]}"; do
      for scaled in "${scaled_vals[@]}"; do

        output_file="$HOME/PWWA/simulations/logstab2tabE2/output_type${type}_dep${dep}_varz${varz}_scaled${scaled}.txt"

        qsub -o "$output_file" \
             -v type=$type,dep=$dep,varz=$varz,scaled=$scaled <<'EOF' &

#!/bin/bash
#PBS -S /bin/bash
#PBS -l select=1:ncpus=96
#PBS -l walltime=100:00:00
#PBS -N SimTab2TabE2
#PBS -j oe

cd "$PBS_O_WORKDIR"

source $HOME/spack-1.0/share/spack/setup-env.sh
spack load r

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

Rscript "$HOME/PWWA/simulations/4_SimTab2TabE2.R" \
        "$type" "$dep" "$varz" "$scaled"

EOF

      done
    done
  done
done
