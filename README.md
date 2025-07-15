# NTT-RVV-CHES2025
NTT vectorization for RISC-V Vector (RVV) - CHES2025

# RISC-V Vector Compiler
The provided code makes use of the [EPI Vector Intrinsics](https://admin.hca.bsc.es/epi/ftp/doc/intrinsics/EPI/epi-intrinsics.html), and therefore requires an RVV 1.0 EPI LLVM compiler toolchain. The easiest way to obtain this is to run the `download-toolchain.sh` from the repository's root, which will download the latest version of the EPI cross-compilation toolchain.

### Other compilers
Other compiler toolchains can be obtained at (BSC's FTP server)[https://ssh.hca.bsc.es/epi/ftp/]. For compiling natively (on a RISC-V machine), look for a `llvm-EPI-development-toolchain-native` toolchain.

# Compilation
To compile the vectorized NTT/INTT implementations, `cd` into `ntt/`, and run `make`. This will create `ntt/bin/ntt-rvv.o`, which will need to be linked with any binaries that make use of OpenFHE

To compile OpenFHE, first set the environment variables to use the EPI toolchain. This can be done with the provided script, after the toolchain has been downloaded: `source set-compilers.sh`. You may then compile OpenFHE normally.