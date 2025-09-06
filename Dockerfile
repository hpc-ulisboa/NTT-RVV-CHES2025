# image for building
FROM ghcr.io/gem5/ubuntu-24.04_all-dependencies:latest AS builder

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        git wget unzip \
        build-essential cmake \
        python3-dev libssl-dev \
        libgmp-dev libmpfr-dev libmpc-dev \
        && rm -rf /var/lib/apt/lists/*

WORKDIR /workspace
RUN git clone --depth 1 https://github.com/hpc-ulisboa/NTT-RVV-CHES2025.git && \
    cd NTT-RVV-CHES2025 && \
    # Replace the SSH URL of the OpenFHE submodule with HTTPS
    git config submodule.openfhe-development.url https://github.com/openfheorg/openfhe-development.git && \
    # Now initialise & update submodules using the corrected URL
    git submodule update --init --recursive

WORKDIR /workspace/NTT-RVV-CHES2025
RUN ./patch-openfhe.sh
RUN ./download-toolchain.sh

WORKDIR /workspace/NTT-RVV-CHES2025/gem5
RUN scons build/RISCV/gem5.opt -j$(nproc)

#   The script only exports variables; we source it and persist the
#   variables for later RUN steps by writing them to a file that we later
#   `source` in each layer.
ENV COMPILER_ENV=/opt/rvv-toolchain/env.sh
RUN mkdir -p /opt/rvv-toolchain && \
    { \
        echo "export TOOLCHAIN_ROOT=/workspace/NTT-RVV-CHES2025/llvm-EPI-release-toolchain-cross"; \
        echo "export CC=\$TOOLCHAIN_ROOT/bin/clang"; \
        echo "export CXX=\$TOOLCHAIN_ROOT/bin/clang++"; \
    } > $COMPILER_ENV

# Use Bash for the rest of the build (so `source` works)
SHELL ["/bin/bash", "-c"]

WORKDIR /workspace/NTT-RVV-CHES2025/openfhe-development
RUN mkdir build && cd build && \
    . $COMPILER_ENV && \
    cmake -DCMAKE_INSTALL_PREFIX=/opt/openfhe -DBUILD_STATIC=ON .. && \
    make -j$(nproc) && \
    make install

WORKDIR /workspace/NTT-RVV-CHES2025/ntt
RUN . $COMPILER_ENV && make

WORKDIR /workspace/NTT-RVV-CHES2025/example-ofhe-app
RUN mkdir -p build && cd build && \
    . $COMPILER_ENV && \
    cmake -DBUILD_STATIC=ON .. && \
    make -j$(nproc)

# lightweight runtime image
FROM ghcr.io/gem5/ubuntu-24.04_all-dependencies:latest AS runtime

# Copy everything we built plus gem5 model
COPY --from=builder /workspace/NTT-RVV-CHES2025/ntt          /opt/NTT-RVV-CHES2025/ntt
COPY --from=builder /workspace/NTT-RVV-CHES2025/example-ofhe-app/build /opt/NTT-RVV-CHES2025/example-ofhe-app/build
COPY --from=builder /workspace/NTT-RVV-CHES2025/gem5/build/RISCV/gem5.opt /opt/gem5/gem5.opt
COPY --from=builder /workspace/NTT-RVV-CHES2025/gem5-model /opt/NTT-RVV-CHES2025/gem5-model
COPY --from=builder /opt/openfhe /opt/openfhe

# Add the RVV tool‑chain to PATH (the tool‑chain script already added it to
# $PATH when sourced; we replicate that here)
ENV PATH=/workspace/NTT-RVV-CHES2025/toolchain/bin:$PATH \
    LD_LIBRARY_PATH=/workspace/NTT-RVV-CHES2025/toolchain/lib:$LD_LIBRARY_PATH

WORKDIR /opt/NTT-RVV-CHES2025

# Users can override this with `docker run … <command>`
CMD ["/bin/bash", "-c", "\
    echo '--- Running RVV NTT speed‑up test ---'; \
    /opt/gem5/gem5.opt /opt/NTT-RVV-CHES2025/gem5-model/main.py /opt/NTT-RVV-CHES2025/ntt/bin/test; \
    echo; \
    echo '--- Running OpenFHE neural‑net example (will take a long time) ---'; \
    /opt/gem5/gem5.opt /opt/NTT-RVV-CHES2025/gem5-model/main.py /opt/NTT-RVV-CHES2025/example-ofhe-app/build/neural-net \
"]