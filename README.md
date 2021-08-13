PyCBG is a python module able to generate [CB-Geo MPM](https://github.com/cb-geo/mpm)'s input files.

Install PyCBG
=============

`pycbg` can be installed using pip which will also install its depedencies. Type the following command from the root of pycbg (the directory where this present file is): 

```
python3 -m pip install pycbg
```

The following modules should be automatically installed: 
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

You can then open `pycbg/doc/_built/index.html` to open the documentation in your browser.