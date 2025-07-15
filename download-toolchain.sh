#!/bin/bash

TOOLCHAIN_ARCHIVE=llvm-EPI-development-toolchain-cross-latest.tar.bz2
URL="https://ssh.hca.bsc.es/epi/ftp/$TOOLCHAIN_ARCHIVE"

wget -N "$URL"

# Download successful?
if [ $? -ne 0 ]; then
    echo "Failed to download the compiler toolchain"
    exit 1
fi

# Extract
echo "Extracting"

tar -xjf "$TOOLCHAIN_ARCHIVE" --verbose

# Extraction successful?
if [ $? -ne 0 ]; then
    echo "Failed to extract the compiler toolchain"
    exit 1
fi

echo "Compiler toolchain has been downloaded and extracted successfully"