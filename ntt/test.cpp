#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "ntt-rvv.h"

#define RDCYCLE(var) __asm__ __volatile__("rdcycle %0" : "=r"(var))

typedef unsigned __int128 uint128_t;

uint64_t *input_scalar;
uint64_t *input_test;
uint64_t *rootOfUnityTable;
uint64_t *rootOfUnityInverseTable;
uint64_t *preconRootOfUnityTable;
uint64_t *preconRootOfUnityInverseTable;

uint32_t msb(uint32_t t)
{
    uint32_t r = 1;

    while (t >>= 1)
    {
        r++;
    }
    return r;
}

uint64_t xorshift_state = 1;
uint64_t generate_random_number_64_bits()
{
    uint64_t x = xorshift_state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    return xorshift_state = x;
}

bool arrays_match(uint64_t *input_scalar, uint64_t *input_test, uint32_t p)
{
    for (uint32_t i = 0; i < p; i++)
        if (input_scalar[i] != input_test[i])
            return false;

    return true;
}

// Cooley-Tuckey
void ntt_scalar(
    const uint32_t p, uint32_t t, uint32_t logt, const uint64_t modulus, uint64_t *element,
    const uint64_t *rootOfUnityTable, const uint64_t *preconRootOfUnityTable)
{
    for (uint32_t m{1}; m < p; m <<= 1, t >>= 1, --logt)
    {
        for (uint32_t i{0}; i < m; ++i)
        {
            uint64_t omega{rootOfUnityTable[i + m]};
            uint64_t preconOmega{preconRootOfUnityTable[i + m]};
            for (uint32_t j1{i << logt}, j2{j1 + t}; j1 < j2; ++j1)
            {
                uint64_t omegaFactor{element[j1 + t]};

                uint64_t q = (((uint128_t)omegaFactor * preconOmega) >> 64);

                uint64_t mult1 = static_cast<uint64_t>(omegaFactor * omega);
                uint64_t mult2 = static_cast<uint64_t>(q * modulus);
                omegaFactor = static_cast<uint64_t>(mult1 - mult2);
                if (omegaFactor >= modulus)
                    omegaFactor -= modulus;

                auto loVal{element[j1 + 0]};

                uint64_t aux2 = loVal + omegaFactor;
                if (aux2 >= modulus)
                    aux2 -= modulus;
                element[j1 + 0] = aux2;

                aux2 = loVal - omegaFactor;
                if (aux2 >= modulus)
                    aux2 += modulus;
                element[j1 + t] = aux2;
            }
        }
    }
}

// Gentelman-Sande
void intt_scalar(
    const uint32_t p, const uint64_t modulus, uint64_t *element, const uint64_t *rootOfUnityInverseTable,
    const uint64_t *preconRootOfUnityInverseTable, const uint64_t cycloOrderInv, const uint64_t preconCycloOrderInv)
{
    for (uint32_t i{0}; i < p; i += 2)
    {
        uint64_t omega{rootOfUnityInverseTable[(i + p) >> 1]};
        uint64_t preconOmega{preconRootOfUnityInverseTable[(i + p) >> 1]};
        uint64_t omegaFactor{element[i + 1]};
        auto loVal{element[i + 0]};

        uint64_t aux2 = loVal + omegaFactor;
        if (aux2 >= modulus)
            aux2 -= modulus;

        uint64_t q = (((uint128_t)aux2 * preconCycloOrderInv) >> 64);

        uint64_t mult1 = static_cast<uint64_t>(aux2 * cycloOrderInv);
        uint64_t mult2 = static_cast<uint64_t>(q * modulus);
        aux2 = static_cast<uint64_t>(mult1 - mult2);
        if (aux2 >= modulus)
            aux2 -= modulus;

        element[i + 0] = aux2;

        aux2 = loVal - omegaFactor;
        if (aux2 >= modulus)
            aux2 += modulus;

        q = (((uint128_t)aux2 * preconOmega) >> 64);

        mult1 = static_cast<uint64_t>(aux2 * omega);
        mult2 = static_cast<uint64_t>(q * modulus);
        aux2 = static_cast<uint64_t>(mult1 - mult2);
        if (aux2 >= modulus)
            aux2 -= modulus;

        q = (((uint128_t)aux2 * preconCycloOrderInv) >> 64);

        mult1 = static_cast<uint64_t>(aux2 * cycloOrderInv);
        mult2 = static_cast<uint64_t>(q * modulus);
        aux2 = static_cast<uint64_t>(mult1 - mult2);
        if (aux2 >= modulus)
            aux2 -= modulus;

        element[i + 1] = aux2;
    }
    for (uint32_t m{p >> 2}, t{2}, logt{2}; m >= 1; m >>= 1, t <<= 1, ++logt)
    {
        for (uint32_t i{0}; i < m; ++i)
        {
            uint64_t omega{rootOfUnityInverseTable[i + m]};
            uint64_t preconOmega{preconRootOfUnityInverseTable[i + m]};

            for (uint32_t j1{i << logt}, j2{j1 + t}; j1 < j2; ++j1)
            {

                auto hiVal{element[j1 + t]};
                auto loVal{element[j1 + 0]};

                uint64_t aux2 = loVal + hiVal;
                if (aux2 >= modulus)
                    aux2 -= modulus;

                element[j1 + 0] = aux2;

                aux2 = loVal - hiVal;
                if (aux2 >= modulus)
                    aux2 += modulus;

                uint64_t q = (((uint128_t)aux2 * preconOmega) >> 64);

                uint64_t mult1 = static_cast<uint64_t>(aux2 * omega);
                uint64_t mult2 = static_cast<uint64_t>(q * modulus);
                aux2 = static_cast<uint64_t>(mult1 - mult2);
                if (aux2 >= modulus)
                    aux2 -= modulus;

                element[j1 + t] = aux2;
            }
        }
    }
}

void test_ntt_wrapper(const uint32_t p, const uint64_t modulus)
{
    const uint32_t t = p >> 1;
    const uint32_t logt = msb(t);

    for (uint32_t i = 0; i < p; i++)
    {
        // Input values in the correct ring
        input_scalar[i] = generate_random_number_64_bits() % modulus;
        input_test[i] = input_scalar[i];

        // random values for the auxiliary tables
        rootOfUnityTable[i] = generate_random_number_64_bits() % modulus;
        preconRootOfUnityTable[i] = generate_random_number_64_bits() % modulus;
    }

    uint64_t start_cycles, end_cycles;

    ntt_scalar(p, t, logt, modulus, input_scalar, rootOfUnityTable, preconRootOfUnityTable);
    RDCYCLE(start_cycles);
    ntt_scalar(p, t, logt, modulus, input_scalar, rootOfUnityTable, preconRootOfUnityTable);
    RDCYCLE(end_cycles);
    uint64_t cycles_scalar = end_cycles - start_cycles;

    ntt_korn_lambiote_vector(p, t, logt, modulus, input_test, rootOfUnityTable, preconRootOfUnityTable);
    RDCYCLE(start_cycles);
    ntt_korn_lambiote_vector(p, t, logt, modulus, input_test, rootOfUnityTable, preconRootOfUnityTable);
    RDCYCLE(end_cycles);
    uint64_t cycles_vector = end_cycles - start_cycles;

    printf("scalar: %lu, vector: %lu cycles\n", cycles_scalar, cycles_vector);
}

void test_intt_wrapper(const uint32_t p, const uint64_t modulus)
{
    uint64_t cycloOrderInv = generate_random_number_64_bits() % modulus;
    uint64_t preconCycloOrderInv = generate_random_number_64_bits() % modulus;

    for (uint32_t i = 0; i < p; i++)
    {
        // Input values in the correct ring
        input_scalar[i] = generate_random_number_64_bits() % modulus;
        input_test[i] = input_scalar[i];

        // random values for the auxiliary tables
        rootOfUnityInverseTable[i] = generate_random_number_64_bits() % modulus;
        preconRootOfUnityInverseTable[i] = generate_random_number_64_bits() % modulus;
    }

    uint64_t start_cycles, end_cycles;

    intt_scalar(p, modulus, input_scalar, rootOfUnityInverseTable, preconRootOfUnityInverseTable, cycloOrderInv, preconCycloOrderInv);
    RDCYCLE(start_cycles);
    intt_scalar(p, modulus, input_scalar, rootOfUnityInverseTable, preconRootOfUnityInverseTable, cycloOrderInv, preconCycloOrderInv);
    RDCYCLE(end_cycles);
    uint64_t cycles_scalar = end_cycles - start_cycles;

    intt_pease_vector_mulh(p, modulus, input_test, rootOfUnityInverseTable, preconRootOfUnityInverseTable, cycloOrderInv, preconCycloOrderInv);
    RDCYCLE(start_cycles);
    intt_pease_vector_mulh(p, modulus, input_test, rootOfUnityInverseTable, preconRootOfUnityInverseTable, cycloOrderInv, preconCycloOrderInv);
    RDCYCLE(end_cycles);
    uint64_t cycles_vector = end_cycles - start_cycles;

    printf("scalar: %lu, vector: %lu cycles\n", cycles_scalar, cycles_vector);
}

int main()
{
    constexpr uint64_t modulus = 10023405405;
    constexpr uint32_t p = 1024 * 1;

    const uint32_t t = p >> 1;
    const uint32_t logt = msb(t);
    const uint32_t stages = logt + 1;

    uint32_t max_gvl = __builtin_epi_vsetvl(p, __epi_e64, __epi_m1);

    // Allocate memory for the arrays
    input_scalar = (uint64_t *)malloc(sizeof(uint64_t) * p);
    input_test = (uint64_t *)malloc(sizeof(uint64_t) * p);
    rootOfUnityTable = (uint64_t *)malloc(sizeof(uint64_t) * p);
    rootOfUnityInverseTable = (uint64_t *)malloc(sizeof(uint64_t) * p);
    preconRootOfUnityTable = (uint64_t *)malloc(sizeof(uint64_t) * p);
    preconRootOfUnityInverseTable = (uint64_t *)malloc(sizeof(uint64_t) * p);

    printf("p: %u\n", p);

    test_ntt_wrapper(p, modulus);
    test_intt_wrapper(p, modulus);

    // Free allocated memory
    free(input_scalar);
    free(input_test);
    free(rootOfUnityTable);
    free(rootOfUnityInverseTable);
    free(preconRootOfUnityTable);
    free(preconRootOfUnityInverseTable);

    return 0;
}