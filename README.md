[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5179973.svg)](https://doi.org/10.5281/zenodo.5179973)

PyCBG is a python module able to generate [CB-Geo MPM](https://github.com/cb-geo/mpm)'s input files.

Install PyCBG
=============

## From PyPI

The last upload of PyCBG on [PyPI](https://pypi.org/) can be simply installed using `pip`:

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

Documentation
=============

## On ReadTheDocs

The latest build of the documentation at [ReadTheDocs](https://readthedocs.org/) is available online [here](https://pycbg.readthedocs.io/en/latest/). If nothing appears under `Classes overview` left of your screen, please reload the page (this is probably a small bug on ReadTheDocs).

## Local build

The documentation can also be built locally using `sphinx`. From the root of PyCBG, run the following command:
```
sh build_doc.sh
```

This will execute a bash script that builds the documentation in a folder named `pycbg_doc`, located in the same directory as pycbg. 
The documentation can then be accesed by opening `pycbg_doc/_build/html/index.html` in your browser.

The pycbg_doc folder name can be modified if necessary when executing the bash script:
```
sh build_doc.sh my_folder_name
```
