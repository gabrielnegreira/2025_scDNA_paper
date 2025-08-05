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

Now we export the bin path to `$PATH` so we can use it
```sh
 export PATH="$VSC_SCRATCH/containers/scistree2/bin:$PATH"
 unset PYTHONPATH
```
