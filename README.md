# wfdiff

Disclaimer: The current of the code is modified from a previous version by [Lion Krischer](https://github.com/krischer).

The code is functional with Python version 3.6 and up, but the test suite is currently broken.

One of the plot capacities is currently provided by the Basemap module, which reached its EoL in 2020.  Although this does not prevent the code from being used with the python 3.8 version, compatibility may not be ensured with future versions.

I would like to re-emphasize the original author's warning:
>***:warning: This package is work in progress and NOT YET READY FOR PRODUCTIVE USE.***:

## Installation

### Installing Conda

`wfdiff` has a couple of dependencies (please check [env_wfdiff.yml](https://github.com/uafgeotools/wfdiff/blob/master/setup.py) if you prefer to manually install dependencies). We recommand that you install wfdiff using `conda`. If Anaconda or (Miniconda) is not available on you system, please download and install [Anaconda](https://www.anaconda.com/products/individual) for your system.

If you have Anaconda already install and you need to update it, you can do so with

```bash
conda update conda
```

### Create the conda environment and install wfdiff

The following will download the latest version of `wfdiff`, create a conda environment (named wfdiff) and install `wfdiff` and all of its dependencies.

```bash
git clone https://github.com/krischer/wfdiff.git
cd wfdiff
conda env create -f env_wfdiff.yml
```

Activate the newly create environment and install wfdiff

 ```bash
 conda activate wfdiff
 ```

 (On windows just with `$ activate wfdiff`). Remember to activate it everytime you want to use `wfdiff`. You can quit the `wfdiff` environment with `conda deactivate`.

You can also update an existing wfdiff environment with

```bash
conda env update -n wfdiff --file env_wfdiff.yml
```

 ---

**Note:** Depending on your cluster, the `mpi4py` shipping with Anaconda might not work with the MPI on your machine. It is best to uninstall the `mpi4py` shipping with Anaconda:

```bash
conda remove mpi4py
```

Now make sure the correct mpi is active (e.g. `mpicc` and consorts point to the correct executables) and install `mpi4py` with

```bash
conda install pip
pip install mpi4py
```

This will cause `mpi4py` to be compiled with the MPI compiler on your system which should resolve any issues.

## Running wfdiff

To run the code, run the `run_wfdiff_test.py` with

```bash
python run_wfdiff_test.py
```

As these calculations can potentially take a long time, you can also run it with MPI:

```bash
mpirun -n 2 python run_wfdiff.py
```

Note that the number of available ressources may differ on your machine.
