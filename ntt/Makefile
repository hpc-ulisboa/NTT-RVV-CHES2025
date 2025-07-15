CC = ../llvm-EPI-development-toolchain-cross/bin/clang++

C_FLAGS = -O2 -mepi
SRC = ntt-rvv.cpp

.PHONY: all
all: bin/ntt-rvv.o

bin/ntt-rvv.o: $(SRC) | $(BIN_DIR)
	mkdir -p bin
	$(CC) $(C_FLAGS) $^ -o $@ -c

.PHONY: clean
clean:
	rm -rf bin