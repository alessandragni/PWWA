#!/bin/bash

# parameter values
theta_vals=(0.5 1 2)
cens_vals=(0.25 0.5) 

mkdir -p logstab1tabE1

for varz in "${theta_vals[@]}"; do
  for cens in "${cens_vals[@]}"; do
    output_file="$HOME/PWWA/simulations/logstab1tabE1/output_varz${varz}_cens${cens}.txt"

    qsub -o "$output_file" \
         -v varz=$varz,cens=$cens <<'EOF' &
         

#!/bin/bash
#PBS -S /bin/bash
#PBS -l select=1:ncpus=96
#PBS -l walltime=100:00:00
#PBS -N SimTab1TabE1
#PBS -j oe

cd "$PBS_O_WORKDIR"

source $HOME/spack-1.0/share/spack/setup-env.sh
spack load r

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

Rscript "$HOME/PWWA/simulations/3_SimTab1TabE1.R" \
        "$varz" "$cens"

EOF

  done
done
