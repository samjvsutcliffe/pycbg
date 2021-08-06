#!/bin/bash

cp -r doc ../${1:-pycbg_doc}
cd ../${1:-pycbg_doc}
make html