#!/bin/bash

rm -fr ${1:-pycbg_doc}
cp -r "${0%/*}"/doc ${1:-pycbg_doc}
cd ${1:-pycbg_doc}
make html