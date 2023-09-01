#!/usr/bin/env bash
#SBATCH --job-name=emapper
#SBATCH -p skylake
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 20
#SBATCH --gres=tmpfs:10G
#SBATCH --ntasks-per-core=1
#SBATCH --time=02:00:00
#SBATCH --mem=50GB
#SBATCH -o /home/a1645424/hpcfs/analysis/popgen/scripts/joblogs/%x_%a_%A_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alastair.ludington@adelaide.edu.au

DIR="/home/a1645424/hpcfs/analysis/popgen"
OUT="${DIR}/results/annotate-outlier-genes"
DB="/home/a1645424/hpcfs/database/funannotate_db/eggnog_db"

mkdir -p "${OUT}"

export EGGNOG_DATA_DIR="${DB}"

source '/home/a1645424/hpcfs/micromamba/etc/profile.d/micromamba.sh'
micromamba activate emapper

emapper.py \
    --cpu "${SLURM_CPUS_PER_TASK}" \
    -i "${OUT}/genes-needing-annotation.fa" \
    --itype 'proteins' \
    --pident 50 \
    --query_cover 70 \
    --output 'annotated-outlier-genes' \
    --output_dir "${OUT}" \
    --scratch_dir "${TMPDIR}" \
    --temp_dir "${TMPDIR}" \
    --override


micromamba deactivate

