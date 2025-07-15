CC = ../llvm-EPI-release-toolchain-cross/bin/clang++

C_FLAGS = -O2 -mepi
SRC = ntt-rvv.cpp

.PHONY: all
all: bin/ntt-rvv.o bin/test

bin/ntt-rvv.o: $(SRC)
	mkdir -p bin
	$(CC) $(C_FLAGS) $^ -o $@ -c

bin/test: test.cpp bin/ntt-rvv.o
	mkdir -p bin
	$(CC) $(C_FLAGS) -static $^ -o $@ -lm

.PHONY: clean
clean:
	rm -rf bin