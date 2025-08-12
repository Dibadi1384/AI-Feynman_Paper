#!/bin/bash
#SBATCH --job-name=Aifeynman-Test
#SBATCH --mail-user=d2daroon@uwaterloo.ca
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --partition=gpu_a100
#SBATCH --gres=gpu:a100:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --time=1-00:00:00
#SBATCH --output=out-%j.log
#SBATCH --chdir=/work/d2daroon/Work1

echo "[$(date)] Job $SLURM_JOB_ID starting on $(hostname)"
echo "Working dir: $(pwd)"
echo "CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES"

# 1) Create & activate virtualenv in cwd (fresh each run)
if [ ! -d feyn ]; then
    python3 -m venv feyn
fi
source feyn/bin/activate

# 2) Install packages with NumPy<2.0 before aifeynman
pip install --upgrade pip
pip install "numpy<2.0"
pip install aifeynman --no-deps
pip install tensorflow sortedcontainers scikit-learn torch seaborn torchvision openpyxl

# 3) Run the Python driver unbuffered
echo "[$(date)] Running run_feynman_all.py"
python3 -u feynex.py
ret=$?
echo "[$(date)] Script exited with code $ret"

deactivate
echo "[$(date)] Job finished"