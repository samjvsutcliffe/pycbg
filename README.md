[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5179973.svg)](https://doi.org/10.5281/zenodo.5179973)

PyCBG is a Python tool that should be helpful in running [CB-Geo MPM](https://github.com/cb-geo/mpm) simulations, either for generating expected input files at the preprocessing stage, or for postprocessing results.

Install PyCBG
=============

`pycbg` can be installed using `pip` (the latter being itself installed on Debian-based systems with `sudo apt install python3-pip`) which will also install its depedencies. Type the following command from the root of `pycbg` (the directory where this present file is): 

```console
$ python3 -m pip install pycbg
```

Installation automatically includes the following dependencies: 
 - `numpy`
 - `gmsh`
 - `pandas`
 - `matplotlib`
 - `versioneer`
 - `pyreadline`
 - `sphinx` (at least version `3.3.1`)
 - `sphinx_rtd_theme`

## Directly from the package

If you want to install a specific version of PyCBG (e.g. the most recent one), you can download the package and install it with `pip`.

**Downloading PyCBG**

You can download the package with `git`:
```console
$ git clone git@forgemia.inra.fr:mpm-at-recover/pycbg.git
```

Alternativly, you can download it manually from [ForgeMia](https://forgemia.inra.fr/mpm-at-recover/pycbg).

**Installing PyCBG**

From the root of the `pycbg` directory (the one you just downloaded), type the following command: 

```console
$ python3 -m pip install -e .
```

Command line usage
==================

While PyCBG is essentially a Python module, installation also provides a new Python executable `pycbg` (`pip` should automatically install it inside a directory in your `$PATH`) with all PyCBG features being already imported. The executable may serve to create a PyCBG interactive session, build the documentation or get PyCBG's version.

## Complete description
```console
$ pycbg -h
usage: pycbg [-h] [-v] [-p] [-i] [-n] [-d [BUILD_DIR]] [PYCBG_SCRIPT]

Manage CB-Geo MPM simulations using PyCBG Python module

positional arguments:
  PYCBG_SCRIPT          pycbg script to be run. By default, the following import lines are added at the top of the file: `from pycbg.preprocessing import *`, `from pycbg.postprocessing import *` and
                        `from pycbg.MPMxDEM import *`. To deactivate this behaviour, use the -n (or --no-import) option

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         print pycbg version
  -p, --pip-show        alias for `pip show pycbg`
  -i, --interactive     run in an interactive IPython session. Using both the -i and -n options simply creates a IPython interactive session
  -n, --no-import       deactivates automatic import of pycbg
  -d [BUILD_DIR], --build-doc [BUILD_DIR]
                        build pycbg's documentation in BUILD_DIR, its path being relative to the current working directory. If BUILD_DIR isn't specified, it will be set to `${PWD}/pycbg_doc`. If
                        BUILD_DIR is `..`, it is set to `../pycbg_doc`. If -d and PYCBG_SCRIPT are specified, the documentation is build before running the script

$ pycbg-gif -h
usage: pycbg-gif [-h] [-c [COLOR_VAR]] [-l [COLOR_LABEL]] [-o [OUTPUT_FILE]] [-j [N_JOBS]] [SIMULATION_DIR] [PROJECTION_PLANE]

Make a gif of material points' positions for the given simulation

positional arguments:
  SIMULATION_DIR        Path to the simulation's directory to be plotted
  PROJECTION_PLANE      Plane on which to project material points' positions. For instance, `1,2` will project positions on the (y,z) plane.

optional arguments:
  -h, --help            show this help message and exit
  -c [COLOR_VAR], --colored-by [COLOR_VAR]
                        Variable to use for material points' coloring. All variables in the csv results files are available, their name being their column's header. Can be a comma-separated list of
                        several variables. If not provided, all material points are black.
  -l [COLOR_LABEL], --color-label [COLOR_LABEL]
                        Label for colorbar, can be a latex expression. If a list of variables was specified with the c option, then it should be a comma-separated list of labels. Default to variables'
                        names.
  -o [OUTPUT_FILE], --output-file [OUTPUT_FILE]
                        Path of the output file, without the '.gif' extension. If a list of variables was specified with the -c option, all file are suffixed with the name of the coloring variable.
                        Default is 'video'.
  -j [N_JOBS], --parallel-jobs [N_JOBS]
                        Number of cores to use. Default to 1.
```

## Usage
One can easily run a PyCBG script:
```console
$ pycbg my_script.py
```

Or experiment with PyCBG's functions and classes:
```console
$ pycbg -i
Python 3.9.7 (default, Sep  3 2021, 12:45:31) 
Type 'copyright', 'credits' or 'license' for more information
IPython 8.0.1 -- An enhanced Interactive Python. Type '?' for help.

In [1]: mesh = Mesh(...)
```

## Get version
To get the installed version of PyCBG, simply run:
```console
$ pycbg -v
v1.0.2+107.g3997bd9
```

## Make GIFs of material point's positions
A GIF of the material points' positions can be made using the `pycbg-gif` command. Each material point can be colored by any variable in the csv files. For instance:
```console
$ pycbg-gif sims_dir 0,2 -c stress_zz,volume -l '$\sigma_{zz}$',V
```
Note that latex strings should be escaped.

Documentation
=============

## On ReadTheDocs

The latest build of the documentation at [ReadTheDocs](https://readthedocs.org/) is available online [here](https://pycbg.readthedocs.io/en/latest/). If nothing appears under `Classes overview` left of your screen, please reload the page (this is probably a small bug on ReadTheDocs).

## Local build using the command line

The documentation can also be built locally using `sphinx`:
```
pycbg -d
```

This will create a folder in the current working directory named `pycbg_doc` containing the documentation's built.
It can then be accesed by opening `pycbg_doc/_build/html/index.html` in your browser.

The `pycbg_doc` folder name can be modified if necessary when executing the bash script:
```
pycbg -d my_folder_name
```

