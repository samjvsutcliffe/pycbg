[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5179973.svg)](https://doi.org/10.5281/zenodo.5179973)

PyCBG is a python module able to generate [CB-Geo MPM](https://github.com/cb-geo/mpm)'s input files.

Install PyCBG
=============

`pycbg` can be installed using `pip` (the latter being itself installed on Debian-based systems with `sudo apt install python3-pip`) which will also install its depedencies. Type the following command from the root of `pycbg` (the directory where this present file is): 

```
python3 -m pip install pycbg
```

Installation automatically includes the following dependencies: 
 - `numpy`
 - `gmsh`
 - `pandas`
 - `matplotlib`
 - `versioneer`
 - `sphinx` (at least version `3.3.1`)
 - `sphinx_rtd_theme`

## Directly from the package

If you want to install a specific version of PyCBG (e.g. the most recent one), you can download the package and install it with `pip`.

**Downloading PyCBG**

You can download the package with `git`:
```
git clone git@forgemia.inra.fr:mpm-at-recover/pycbg.git
```

Alternativly, you can download it manually from [ForgeMia](https://forgemia.inra.fr/mpm-at-recover/pycbg).

**Installing PyCBG**

From the root of the `pycbg` directory (the one you just downloaded), type the following command: 

```
python3 -m pip install -e .
```

Command line usage
==================

After installing PyCBG, a new python executable `pycbg` is available (`pip` should automatically install it inside a directory in your `$PATH`). It allows to easily create a PyCBG interactive session, build the documentation or getting PyCBG's version. Here is its description:
```console
$ pycbg -h
usage: pycbg [-h] [-v] [-i] [-n] [-d [BUILD_DIR]] [PYCBG_SCRIPT]

Manage CG-Geo MPM simulations using PyCBG Python module

positional arguments:
  PYCBG_SCRIPT          pycbg script to be run. By default, the following import lines are added at the top of the file: `from pycbg.preprocessing import *`, `from pycbg.postprocessing import *`. To
                        deactivate this behaviour, use the -n (or --no-import) option

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         print pycbg version
  -i, --interactive     run in an interactive IPython session. Using both the -i and -n options simply creates a IPython interactive session
  -n, --no-import       deactivates automatic import of pycbg when running PYCBG_SCRIPT
  -d [BUILD_DIR], --build-doc [BUILD_DIR]
                        build pycbg's documentation in BUILD_DIR, its path being relative to the current working directory. If the directory already exists, it is removed without prompt before building
                        the doc. If BUILD_DIR isn't specified, it will be set to `${PWD}/pycbg_doc`. If -d and PYCBG_SCRIPT are specified, the documentation is build before running the script
```

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
sh build_doc.sh my_folder_name
```

