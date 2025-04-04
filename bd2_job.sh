#!/bin/bash
#SBATCH --job-name=np_test_job        # Job name
#SBATCH --output=output.txt           # Standard output file
#SBATCH --error=error.txt             # Standard error file
#SBATCH --partition=batch             # Partition or queue name
#SBATCH --nodes=4                     # Number of nodes
#SBATCH --ntasks 256                  # Number of tasks for distributed parallelism
#SBATCH --ntasks-per-node=64          # Number of tasks per node for MPI apps
#SBATCH --cpus-per-task=1             # Number of CPU cores per task for non-MPI / multithreaded apps
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=jewebste@butler.edu

cd /home/jewebste/CarmichaelAlgorithms

mpirun ./parallel
