Install pycbg
=============

`pycbg` can be installed as a Python3 module using pip. Type the following command from the root of pycbg (the directory where this present file is): 

```
python3 -m pip install -e .
```

Installation should automatically include the following dependencies: 
 - `numpy`
 - `gmsh`
 - `pandas`
 - `matplotlib`
 - `sphinx` (at least version `3.3.1`)
 - `sphinx_rtd_theme`

Build the documentation
=======================

The documentation can be build using `sphinx`. From the root of pycbg, run the following command:
```
sh build_doc.sh
```

This will execute a bash script that builds the documentation in a folder named `pycbg_doc`, located in the same directory as pycbg. 
The documentation can then be accesed by opening `pycbg_doc/_build/html/index.html` in your browser.

The name of the folder can be specified when executing the bash script:
```
sh build_doc.sh my_folder_name
```