wfdiff
======

Disclaimer: The current version of the code by [Lion Krischer](https://github.com/krischer/wfdiff) does not work under the originally planned python version (2.7, 3.3, 3.4). A quick workaround is to install it on python 3.8 (the Obspy module is not yet available for python 3.9).
This workaround also seems to work in python 3.6 and 3.7.

Although the code seems to be fully functional when installed with version 3.8, the test suite currently appears to be broken.

One of the plot capacities is currently provided by the Basemap module, which reached its EoL in 2020.  Although this does not prevent the code from being used with the python 3.8 version, compatibility may not be ensured with future versions.

I would like to re-emphasize the original author's warning:
>***:warning: This package is work in progress and NOT YET READY FOR PRODUCTIVE USE.***:

## Installation
### Installing Python and Dependencies

Installation instructions are slightly modified from original doc [here](http://krischer.github.io/wfdiff/). 

`wfdiff` has a couple of dependencies. If you are well versed in Python, then the installation should not prove an issue, otherwise, please follow the advice here and download the [Anaconda Python distribution](https://www.anaconda.com/products/individual) for your system. I encourage chosing Python 3.8 at the moment, but 3.6 and 3.7 should also work with these instructions.

After downloading and installing Anaconda, update it with

```
$ conda update conda
```
Then create a new environment. You don’t have to, but using a new environment grants a clean seperation between Python installations.

```
 $ conda create -n wfdiff python=3.8
 ```
 
 This will create a new conda environment based on Python 3.8 named `wfdiff`. Activate it with
 
 ```
 $ conda activate wfdiff
 ```
 
 (On windows just with `$ activate wfdiff`). Remember to activate it everytime you want to use `wfdiff`. You can quit the `wfdiff` environment with `conda deactivate`.
 
 Now install ObsPy with
 
 ```
 $ conda install -c obspy obspy
 ```
 and the remaining dependencies with
 
 ```
 $ conda install basemap pandas flake8 pytest nose mpi4py tqdm pyasdf
 $ conda install -c conda-forge basemap-data-hires
 ```
---
**Note:** Depending on your cluster, the `mpi4py` shipping with Anaconda might not work with the MPI on your machine. It is best to uninstall the `mpi4py` shipping with Anaconda:

```
$ conda remove mpi4py
```

Now make sure the correct mpi is active (e.g. `mpicc` and consorts point to the correct executables) and install `mpi4py` with

```
$ conda install pip
$ pip install mpi4py
```

This will cause `mpi4py` to be compiled with the MPI compiler on your system which should resolve any issues.

### Install wfdiff
You should work with the latest version of the code and install it with

```
$ git clone https://github.com/krischer/wfdiff.git
$ cd wfdiff
$ pip install -v -e .
```
## Running wfdiff
To run the code, run the `run_wfdiff_test.py` with
```
$ python run_wfdiff_test.py
```

As these calculations can potentially take a long time, you can also run it with MPI:
```
$ mpirun -n 8 python run_wfdiff.py
```
Note that the number of available ressources may differ on your machine.

The public interfaces of `wfdiff` can correctly deal with MPI and will parallelize the code if possible. Please make sure to not run anything else in the Python file or be aware of what is actually happening.
