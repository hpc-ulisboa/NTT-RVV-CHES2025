#!/bin/bash

if [ "$0" = "$BASH_SOURCE" ]; then
    echo "Error: Don't run this script, source it"
    exit 1
fi

TOOLCHAIN="llvm-EPI-development-toolchain-cross"

# Set CC and CXX
export CC="$(pwd)/$TOOLCHAIN/bin/clang"
export CXX="$(pwd)/$TOOLCHAIN/bin/clang++"

echo "CC and CXX set to:"
echo "CC: $CC"
echo "CXX: $CXX"