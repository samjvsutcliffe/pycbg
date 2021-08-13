PyCBG is a python module able to generate [CB-Geo MPM](https://github.com/cb-geo/mpm)'s input files.

Install PyCBG
=============

## From PyPI

The last upload of PyCBG on PyPI can be simply installed using pip:

```
python3 -m pip install pycbg
```

Installation should automatically include the following dependencies: 
 - `numpy`
 - `gmsh`
 - `pandas`
 - `matplotlib`
 - `versioneer`
 - `sphinx` (at least version `3.3.1`)
 - `sphinx_rtd_theme`

## Directly from the package

If you want to install a specific version of PyCBG (e.g. the most recent one), you can download the package and install it with pip.

**Downloading PyCBG**

You can download the package with git:
```
git clone git@forgemia.inra.fr:mpm-at-recover/pycbg.git
```

Altenativly, you can download it manually from ForgeMia.

**Installing PyCBG**

From the root of the `pycbg` directory (the one you just downloaded), type the following command: 

```
python3 -m pip install -e .
```

Documentation
=============

## On ReadTheDocs

The latest build of the documentation is available online [here](https://pycbg.readthedocs.io/en/latest/). If nothing appears under `Classes overview`, please reload the page (this is probably a small bug on ReadTheDocs).

## Local built

The documentation can also be build locally using `sphinx`. From the root of PyCBG, run the following command:
```
sh build_doc.sh
```

This will execute a bash script that builds the documentation in a folder named `pycbg_doc`, located in the same directory as pycbg. 
The documentation can then be accesed by opening `pycbg_doc/_build/html/index.html` in your browser.

The pycbg_doc folder name can be modified when executing the bash script:
```
sh build_doc.sh my_folder_name
```
