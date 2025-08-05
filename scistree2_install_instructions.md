Instead of using `conda` environments, it is best to use `containers` in the HPC. This is because `conda` environments create several small files which are not optimal for HPC structure.

To do so, we first create a configuration file in `yaml`:
`scistree_env.yml`
```yml
name: scistree2
channels:
    - bioconda
    - conda-forge
dependencies:
    - python=3.11
    - pandas
    - numpy
    - tree
    - ete3
    - pyqt=5
```

Then we load `hpc-container-wrapper` and run `conda-containerize` to create the container. Here we store it in the `scratch` directory:

```sh
module load hpc-container-wrapper
conda-containerize new --prefix "$VSC_SCRATCH/containers/scistree2" scistree2_env.yml
```

Since [ScisTree2](https://github.com/yufengwudcs/ScisTree2/tree/python) is not available via conda, we have to manually install it. To do so, we can use the `conda-containerize update` command and provide a `post-install` script to it:

`scistree2_install.sh`
```sh
git clone https://github.com/yufengwudcs/ScisTree2.git
cd ScisTree2
pip install .
```

Then run:
```sh
conda-containerize update --post-install scistree2_install.sh "$VSC_SCRATCH/containers/scistree2"
```
After it is complete we can remove the `ScisTree2` directory that was created when the repo was cloned:
```sh
rm -rf ScisTree2
```

To use it, we export the bin path to `$PATH` so we can use it
```sh
 export PATH="$VSC_SCRATCH/containers/scistree2/bin:$PATH"
 unset PYTHONPATH
```

# Running on HPC
since we have many files, we can use slurm arrays for that. To do so, we write a script that will be submitted to slurms and in that script we specify all files:

scistree_run_HPC.sh
```sh
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=50:00:00
#SBATCH --job-name=scistree2_run
#SBATCH --account=ap_itg_mpu
#SBATCH --mail-user=gnegreira@itg.be
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err
#SBATCH --array=1-NUM_JOBS  # will be replaced dynamically

# Load modules
module load calcua/all
module load hpc-container-wrapper

# Activate container
export PATH="$VSC_SCRATCH/containers/scistree2/bin:$PATH"
unset PYTHONPATH

# Create a sorted list of files
FILE_LIST=($(ls "$VSC_SCRATCH/atrandi_paper/inputs/nucleotide_variants"/*for_scistree2.tsv))

# Select the file corresponding to this array task
FILE=${FILE_LIST[$SLURM_ARRAY_TASK_ID-1]}

echo "Processing file: $FILE"
python ScisTree2.py "$FILE"
```

To submit it, we just need to provide the array to `sbatch`:
```sh
chmod 777 scistree2_run_HPC.sh
NUM_FILES=$(ls "$VSC_SCRATCH/atrandi_paper/inputs/nucleotide_variants"/*for_scistree2.tsv | wc -l)
sbatch --array=1-$NUM_FILES scistree2_run_HPC.sh
```
Slurm should take care of the rest