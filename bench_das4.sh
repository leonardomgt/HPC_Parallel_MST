#!/bin/bash
sbatch --nodes=1 main.job /var/scratch/hpcl1918/mouse_gene.mtx
sbatch --nodes=2 main.job /var/scratch/hpcl1918/mouse_gene.mtx
sbatch --nodes=4 main.job /var/scratch/hpcl1918/mouse_gene.mtx
sbatch --nodes=8 main.job /var/scratch/hpcl1918/mouse_gene.mtx
sbatch --nodes=16 main.job /var/scratch/hpcl1918/mouse_gene.mtx

sbatch --nodes=1 main.job /var/scratch/hpcl1918/ldoor.mtx
sbatch --nodes=2 main.job /var/scratch/hpcl1918/ldoor.mtx
sbatch --nodes=4 main.job /var/scratch/hpcl1918/ldoor.mtx
sbatch --nodes=8 main.job /var/scratch/hpcl1918/ldoor.mtx
sbatch --nodes=16 main.job /var/scratch/hpcl1918/ldoor.mtx

sbatch --nodes=1 main.job /var/scratch/hpcl1918/nlpkkt240.mtx
sbatch --nodes=2 main.job /var/scratch/hpcl1918/nlpkkt240.mtx
sbatch --nodes=4 main.job /var/scratch/hpcl1918/nlpkkt240.mtx
sbatch --nodes=8 main.job /var/scratch/hpcl1918/nlpkkt240.mtx
sbatch --nodes=16 main.job /var/scratch/hpcl1918/nlpkkt240.mtx

squeue