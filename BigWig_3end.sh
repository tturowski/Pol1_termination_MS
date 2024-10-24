#!/bin/bash
#SBATCH --job-name=BigWig_3end_nc
#SBATCH --output=BigWig_3end_nc_%j.out
#SBATCH --error=BigWig_3end_nc_%j.err
#SBATCH --nodes=12
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=72:00:00

# Load any required modules
# module load your_module_name

# Activate the Conda environment
source /home/${USER}/.bashrc
eval "$(conda shell.bash hook)"
conda activate processing

# Navigate to the directory containing your SAM files
cd /home/tomasz.turowski/99_manuscripts/Pol1_termination_MS/04_BigWig

# Loop through each SAM file and run SAM2profilesGenomic.py
for f in *sam; do
    SAM2profilesGenomic.py -f $f -u 3end -n -s polyA &
<<<<<<< HEAD
#    SAM2profilesGenomic.py -f $f -u 5end &
=======
    SAM2profilesGenomic.py -f $f -u 5end &
>>>>>>> 6a552995de66d68ed29acc0378f0d70b8a2c35fe
done

# Wait for all background jobs to finish
wait

# Deactivate the Conda environment
conda deactivate
