#!/bin/bash

if [ "$0" = "$BASH_SOURCE" ]; then
    echo "Error: Don't run this script, source it"
    exit 1
fi

TOOLCHAIN="llvm-EPI-release-toolchain-cross"

# Get the directory of the script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Set CC and CXX relative to the script location
export CC="$SCRIPT_DIR/$TOOLCHAIN/bin/clang"
export CXX="$SCRIPT_DIR/$TOOLCHAIN/bin/clang++"

echo "CC and CXX set to:"
echo "CC: $CC"
echo "CXX: $CXX"