Install pycbg
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
 - `sphinx` (at least version `3.3.1`)
 - `sphinx_rtd_theme`

Build the documentation
=======================

The documentation can be build using `sphinx`. From the root of pycbg, run the following command :

```
sh build_doc.sh
```

You can then open `pycbg/doc/_built/index.html` to open the documentation in your browser.