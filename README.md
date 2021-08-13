PyCBG is a python module able to generate [CB-Geo MPM](https://github.com/cb-geo/mpm)'s input files.

Install PyCBG
=============

`pycbg` can be installed using pip which will also install its depedencies. Type the following command from the root of pycbg (the directory where this present file is): 

```
python3 -m pip install -e .
```

The following modules should be automatically installed: 
 - `numpy`
 - `gmsh`
 - `pandas`
 - `matplotlib`
 - `versioneer`
 - `sphinx` (at least version `3.3.1`)
 - `sphinx_rtd_theme`

Documentation
=============

## On ReadTheDocs

The latest build of the documentation is available online [here](https://pycbg.readthedocs.io/en/latest/). If nothing appears under `Classes overview`, please reload the page (this is probably a small bug on ReadTheDocs).

## Build documentation locally

The documentation can also be build locally using `sphinx`. From the root of PyCBG, run the following command:
```
sh build_doc.sh
```

You can then open `pycbg/doc/_built/index.html` to open the documentation in your browser.