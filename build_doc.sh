#!/bin/bash

if [[ ${1:-pycbg_doc} -ef .. ]]
then
    BUILD_DIR="${1:-pycbg_doc}/pycbg_doc"
else
    BUILD_DIR=${1:-pycbg_doc}
fi

build_doc () {
    echo "build dir: ${1}"
    rm -fr ${1}
    cp -r "${0%/*}"/doc ${1}
    cd ${1}
    make html
}

if [ -d "$BUILD_DIR" ]
then
    read -r -p "The directory $BUILD_DIR already exists. Overwrite it? [y/N] " response
    case "$response" in
        [yY][eE][sS]|[yY]) 
            build_doc $BUILD_DIR
            ;;
        *)
            echo "Aborting documention's build"
            ;;
    esac
else
    build_doc $BUILD_DIR
fi

