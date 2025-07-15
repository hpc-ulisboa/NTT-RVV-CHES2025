
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ntt-rvv.h"

uint8_t initialized_mem = 1;
uint64_t *aux_increment = NULL;
uint64_t *auxiliary_array = NULL;

typedef unsigned __int128 uint128_t;

void free_ntts_mem()
{
    if (!initialized_mem)
    {
        free(aux_increment);
        free(auxiliary_array);
    }
    initialized_mem = 1;
}

void ntt_korn_lambiote_vector(
    const uint32_t p, uint32_t t, uint32_t logt, const uint64_t modulus, uint64_t *element,
    const uint64_t *rootOfUnityTable, const uint64_t *preconRootOfUnityTable)
{
    uint8_t element_or_aux = 0;
    uint32_t stopvalue_sft1 = p >> 1;
    uint32_t gvl = __builtin_epi_vsetvl(p, __epi_e64, __epi_m1);
    uint32_t max_gvl = __builtin_epi_vsetvl(p, __epi_e64, __epi_m1);
    uint64_t *input_stage_array = element;
    uint64_t *output_stage_array = NULL;
    uint64_t *swap_pointer_aux = NULL;

    uint32_t log_of_gvl = log2(gvl);

    __epi_1xi64 v_coef_mod = __builtin_epi_vmv_v_x_1xi64(modulus, gvl);
    __epi_1xi64 v_shift_32 = __builtin_epi_vmv_v_x_1xi64(32, gvl);
    __epi_1xi64 v_32_dimension = __builtin_epi_vmv_v_x_1xi64(4294967295, gvl);
    __epi_1xi64 v_U;
    __epi_1xi64 v_V;

    if (initialized_mem)
    {
        initialized_mem = 0;
        aux_increment = (uint64_t *)malloc(gvl * sizeof(uint64_t));
        auxiliary_array = (uint64_t *)malloc(p * sizeof(uint64_t));
        for (uint32_t i = 0; i < gvl; i++)
        {
            aux_increment[i] = i;
        }
    }

    output_stage_array = auxiliary_array;
    __epi_1xi64 v_index_1 = __builtin_epi_vmv_v_x_1xi64(1, gvl);
    __epi_1xi64 v_index_original = __builtin_epi_vload_unsigned_1xi64(&aux_increment[0], gvl);
    __epi_1xi64 v_index_stage_aux = __builtin_epi_vmv_v_x_1xi64(1, gvl);
    __epi_1xi64 v_m = __builtin_epi_vmv_v_x_1xi64(2, gvl);

    uint32_t m = 1;

    __epi_1xi64 v_S = __builtin_epi_vmv_v_x_1xi64(rootOfUnityTable[1], gvl);
    __epi_1xi64 v_precomp_aux = __builtin_epi_vmv_v_x_1xi64(preconRootOfUnityTable[1], gvl);

    for (uint32_t i = 0; i < p; i += (max_gvl << 1)) // First Loop (Broadcast)
    {
        v_U = __builtin_epi_vload_unsigned_1xi64(&(input_stage_array[(i >> 1)]), max_gvl);
        v_V = __builtin_epi_vload_unsigned_1xi64(&(input_stage_array[(i >> 1) + stopvalue_sft1]), max_gvl);

        __epi_1xi64 v_q = __builtin_epi_vmulh_1xi64(v_V, v_precomp_aux, max_gvl);

        __epi_1xi64 v_mult = __builtin_epi_vmul_1xi64(v_V, v_S, max_gvl);
        __epi_1xi64 v_aux = __builtin_epi_vmul_1xi64(v_q, v_coef_mod, max_gvl);
        v_V = __builtin_epi_vsub_1xi64(v_mult, v_aux, max_gvl);
        __epi_1xi1 mask = __builtin_epi_vmsleu_1xi64(v_coef_mod, v_V, max_gvl);
        v_V = __builtin_epi_vsub_1xi64_mask(v_V, v_V, v_coef_mod, mask, max_gvl);

        __epi_1xi64 v_result_0 = __builtin_epi_vadd_1xi64(v_U, v_V, max_gvl);

        // Reduction v_result_0
        mask = __builtin_epi_vmsleu_1xi64(v_coef_mod, v_result_0, max_gvl);
        v_result_0 = __builtin_epi_vsub_1xi64_mask(v_result_0, v_result_0, v_coef_mod, mask, max_gvl);

        // store results
        __builtin_epi_vstore_strided_unsigned_1xi64(&(output_stage_array[i]), v_result_0, 16, max_gvl);

        __epi_1xi64 v_result_1 = __builtin_epi_vsub_1xi64(v_U, v_V, max_gvl); // masked as well? for speedup?
        // Reduction v_result_1
        mask = __builtin_epi_vmsleu_1xi64(v_coef_mod, v_result_1, max_gvl);
        v_result_1 = __builtin_epi_vadd_1xi64_mask(v_result_1, v_result_1, v_coef_mod, mask, max_gvl);

        // store results
        __builtin_epi_vstore_strided_unsigned_1xi64(&(output_stage_array[i + 1]), v_result_1, 16, max_gvl);
    }

    m <<= 1;
    element_or_aux ^= 1;
    __epi_1xi64 v_omega_og = __builtin_epi_vload_unsigned_1xi64(&rootOfUnityTable[0], max_gvl);
    __epi_1xi64 v_precomp_og = __builtin_epi_vload_unsigned_1xi64(&preconRootOfUnityTable[0], max_gvl);

    swap_pointer_aux = output_stage_array;
    output_stage_array = input_stage_array;
    input_stage_array = swap_pointer_aux;

    for (; m < gvl; m <<= 1, element_or_aux ^= 1) // Intermediate Loops vrgather
    {
        __epi_1xi64 v_index = __builtin_epi_vand_1xi64(v_index_original, v_index_stage_aux, gvl);
        v_index = __builtin_epi_vadd_1xi64(v_index, v_m, gvl);
        __epi_1xi64 v_S = __builtin_epi_vrgather_1xi64(v_omega_og, v_index, gvl);
        __epi_1xi64 v_precomp_aux = __builtin_epi_vrgather_1xi64(v_precomp_og, v_index, gvl);

        v_index_stage_aux = __builtin_epi_vsll_1xi64(v_index_stage_aux, v_index_1, gvl);
        v_index_stage_aux = __builtin_epi_vor_1xi64(v_index_stage_aux, v_index_1, gvl);
        v_m = __builtin_epi_vsll_1xi64(v_m, v_index_1, gvl);

        for (uint32_t i = 0; i < p; i += (max_gvl << 1))
        {
            v_U = __builtin_epi_vload_unsigned_1xi64(&(input_stage_array[(i >> 1)]), max_gvl);
            v_V = __builtin_epi_vload_unsigned_1xi64(&(input_stage_array[(i >> 1) + stopvalue_sft1]), max_gvl);

            __epi_1xi64 v_q = __builtin_epi_vmulh_1xi64(v_V, v_precomp_aux, max_gvl);

            __epi_1xi64 v_mult = __builtin_epi_vmul_1xi64(v_V, v_S, max_gvl);
            __epi_1xi64 v_aux = __builtin_epi_vmul_1xi64(v_q, v_coef_mod, max_gvl);
            v_V = __builtin_epi_vsub_1xi64(v_mult, v_aux, max_gvl);
            __epi_1xi1 mask = __builtin_epi_vmsleu_1xi64(v_coef_mod, v_V, max_gvl);
            v_V = __builtin_epi_vsub_1xi64_mask(v_V, v_V, v_coef_mod, mask, max_gvl);

            __epi_1xi64 v_result_0 = __builtin_epi_vadd_1xi64(v_U, v_V, max_gvl);

            // Reduction v_result_0
            mask = __builtin_epi_vmsleu_1xi64(v_coef_mod, v_result_0, max_gvl);
            v_result_0 = __builtin_epi_vsub_1xi64_mask(v_result_0, v_result_0, v_coef_mod, mask, max_gvl);

            __epi_1xi64 v_result_1 = __builtin_epi_vsub_1xi64(v_U, v_V, max_gvl); // masked as well? for speedup?
            // Reduction v_result_1
            mask = __builtin_epi_vmsleu_1xi64(v_coef_mod, v_result_1, max_gvl);
            v_result_1 = __builtin_epi_vadd_1xi64_mask(v_result_1, v_result_1, v_coef_mod, mask, max_gvl);

            // store results
            __builtin_epi_vstore_strided_unsigned_1xi64(&(output_stage_array[i]), v_result_0, 16, max_gvl);
            __builtin_epi_vstore_strided_unsigned_1xi64(&(output_stage_array[i + 1]), v_result_1, 16, max_gvl);
        }

        swap_pointer_aux = output_stage_array;
        output_stage_array = input_stage_array;
        input_stage_array = swap_pointer_aux;
    }

    for (; m < p; m <<= 1, element_or_aux ^= 1) // Final Loops vload(with skips)
    {
        for (uint32_t i = 0; i < (m >> log_of_gvl); i++) // m/256       m>>8
        {
            uint64_t index = i << log_of_gvl;

            __epi_1xi64 v_S = __builtin_epi_vload_unsigned_1xi64(&rootOfUnityTable[index + m], max_gvl);
            __epi_1xi64 v_precomp_aux = __builtin_epi_vload_unsigned_1xi64(&preconRootOfUnityTable[index + m], max_gvl);

            for (uint32_t j = 0; j < (stopvalue_sft1 / m); j++)
            {
                v_U = __builtin_epi_vload_unsigned_1xi64(&(input_stage_array[index + j * m]), max_gvl);
                v_V = __builtin_epi_vload_unsigned_1xi64(&(input_stage_array[index + j * m + stopvalue_sft1]), max_gvl);

                __epi_1xi64 v_q = __builtin_epi_vmulh_1xi64(v_V, v_precomp_aux, max_gvl);

                __epi_1xi64 v_mult = __builtin_epi_vmul_1xi64(v_V, v_S, max_gvl);
                __epi_1xi64 v_aux = __builtin_epi_vmul_1xi64(v_q, v_coef_mod, max_gvl);
                v_V = __builtin_epi_vsub_1xi64(v_mult, v_aux, max_gvl);
                __epi_1xi1 mask = __builtin_epi_vmsleu_1xi64(v_coef_mod, v_V, max_gvl);
                v_V = __builtin_epi_vsub_1xi64_mask(v_V, v_V, v_coef_mod, mask, max_gvl);

                __epi_1xi64 v_result_0 = __builtin_epi_vadd_1xi64(v_U, v_V, max_gvl);

                // Reduction v_result_0
                mask = __builtin_epi_vmsleu_1xi64(v_coef_mod, v_result_0, max_gvl);
                v_result_0 = __builtin_epi_vsub_1xi64_mask(v_result_0, v_result_0, v_coef_mod, mask, max_gvl);

                __epi_1xi64 v_result_1 = __builtin_epi_vsub_1xi64(v_U, v_V, max_gvl); // masked as well? for speedup?
                // Reduction v_result_1
                mask = __builtin_epi_vmsleu_1xi64(v_coef_mod, v_result_1, max_gvl);
                v_result_1 = __builtin_epi_vadd_1xi64_mask(v_result_1, v_result_1, v_coef_mod, mask, max_gvl);

                // store results
                __builtin_epi_vstore_strided_unsigned_1xi64(&(output_stage_array[((index + j * m) << 1)]), v_result_0, 16, max_gvl);
                __builtin_epi_vstore_strided_unsigned_1xi64(&(output_stage_array[((index + j * m) << 1) + 1]), v_result_1, 16, max_gvl);
            }
        }

        swap_pointer_aux = output_stage_array;
        output_stage_array = input_stage_array;
        input_stage_array = swap_pointer_aux;
    }

    // for half of the m values, the results will end up in the auxiliary array rather than the element array
    if (element_or_aux)
    {
        for (uint32_t i = 0; i < p; i += max_gvl)
        {
            __epi_1xi64 v_auxarr = __builtin_epi_vload_unsigned_1xi64(&(input_stage_array[i]), max_gvl);
            __builtin_epi_vstore_unsigned_1xi64(&element[i], v_auxarr, max_gvl);
        }
    }
}

void intt_pease_vector_mulh(
    const uint32_t p, const uint64_t modulus, uint64_t *element, const uint64_t *rootOfUnityInverseTable,
    const uint64_t *preconRootOfUnityInverseTable, const uint64_t cycloOrderInv, const uint64_t preconCycloOrderInv)
{

    uint8_t element_or_aux = 0;
    uint32_t stopvalue_sft1 = p >> 1;
    uint32_t gvl = __builtin_epi_vsetvl(p, __epi_e64, __epi_m1);
    uint32_t max_gvl = __builtin_epi_vsetvl(p, __epi_e64, __epi_m1);
    __epi_1xi64 v_coef_mod = __builtin_epi_vmv_v_x_1xi64(modulus, gvl);
    __epi_1xi64 v_shift_32 = __builtin_epi_vmv_v_x_1xi64(32, gvl);
    __epi_1xi64 v_32_dimension = __builtin_epi_vmv_v_x_1xi64(4294967295, gvl);
    __epi_1xi64 v_preconCycloOrderInv = __builtin_epi_vmv_v_x_1xi64(preconCycloOrderInv, gvl);
    __epi_1xi64 v_cycloOrderInv = __builtin_epi_vmv_v_x_1xi64(cycloOrderInv, gvl);
    uint64_t *input_stage_array = element;
    uint64_t *output_stage_array = NULL;
    uint64_t *swap_pointer_aux = NULL;

    uint32_t log_of_gvl = log2(gvl);

    if (initialized_mem)
    {
        initialized_mem = 0;
        aux_increment = (uint64_t *)malloc(gvl * sizeof(uint64_t));
        auxiliary_array = (uint64_t *)malloc(p * sizeof(uint64_t));
        for (uint32_t i = 0; i < gvl; i++)
        {
            aux_increment[i] = i;
        }
    }

    output_stage_array = auxiliary_array;

    __epi_1xi64 v_index_1 = __builtin_epi_vmv_v_x_1xi64(1, gvl);
    __epi_1xi64 v_index_original = __builtin_epi_vload_unsigned_1xi64(&aux_increment[0], gvl);
    __epi_1xi64 v_index_stage_aux = __builtin_epi_vmv_v_x_1xi64((gvl >> 1) - 1, gvl);
    __epi_1xi64 v_m = __builtin_epi_vmv_v_x_1xi64(gvl >> 1, gvl);

    uint32_t m = p >> 1;

    for (uint32_t i = 0; i < (m >> log_of_gvl); i++) // m/256
    {
        uint64_t index = i << log_of_gvl;

        __epi_1xi64 v_S = __builtin_epi_vload_unsigned_1xi64(&rootOfUnityInverseTable[index + m], max_gvl);
        __epi_1xi64 v_precomp_aux = __builtin_epi_vload_unsigned_1xi64(&preconRootOfUnityInverseTable[index + m], max_gvl);

        __epi_1xi64 v_U = __builtin_epi_vload_strided_unsigned_1xi64(&(input_stage_array[((index) << 1)]), 16, max_gvl);
        __epi_1xi64 v_V = __builtin_epi_vload_strided_unsigned_1xi64(&(input_stage_array[((index) << 1) + 1]), 16, max_gvl);

        __epi_1xi64 v_result_0 = __builtin_epi_vadd_1xi64(v_U, v_V, max_gvl);

        // Reduction v_result_0
        __epi_1xi1 mask = __builtin_epi_vmsleu_1xi64(v_coef_mod, v_result_0, max_gvl);
        v_result_0 = __builtin_epi_vsub_1xi64_mask(v_result_0, v_result_0, v_coef_mod, mask, max_gvl);

        __epi_1xi64 v_q = __builtin_epi_vmulh_1xi64(v_V, v_preconCycloOrderInv, max_gvl);

        __epi_1xi64 v_mult = __builtin_epi_vmul_1xi64(v_result_0, v_cycloOrderInv, max_gvl);
        __epi_1xi64 v_aux = __builtin_epi_vmul_1xi64(v_q, v_coef_mod, max_gvl);
        v_result_0 = __builtin_epi_vsub_1xi64(v_mult, v_aux, max_gvl);
        mask = __builtin_epi_vmsleu_1xi64(v_coef_mod, v_result_0, max_gvl);
        v_result_0 = __builtin_epi_vsub_1xi64_mask(v_result_0, v_result_0, v_coef_mod, mask, max_gvl);

        __epi_1xi64 v_result_1 = __builtin_epi_vsub_1xi64(v_U, v_V, max_gvl); // masked as well? for speedup?
        // Reduction v_result_1
        mask = __builtin_epi_vmsleu_1xi64(v_coef_mod, v_result_1, max_gvl);
        v_result_1 = __builtin_epi_vadd_1xi64_mask(v_result_1, v_result_1, v_coef_mod, mask, max_gvl);

        v_q = __builtin_epi_vmulh_1xi64(v_result_1, v_precomp_aux, max_gvl);

        v_mult = __builtin_epi_vmul_1xi64(v_result_1, v_S, max_gvl);
        v_aux = __builtin_epi_vmul_1xi64(v_q, v_coef_mod, max_gvl);
        v_result_1 = __builtin_epi_vsub_1xi64(v_mult, v_aux, max_gvl);
        mask = __builtin_epi_vmsleu_1xi64(v_coef_mod, v_result_1, max_gvl);
        v_result_1 = __builtin_epi_vsub_1xi64_mask(v_result_1, v_result_1, v_coef_mod, mask, max_gvl);

        v_q = __builtin_epi_vmulh_1xi64(v_result_1, v_preconCycloOrderInv, max_gvl);

        v_mult = __builtin_epi_vmul_1xi64(v_result_1, v_cycloOrderInv, max_gvl);
        v_aux = __builtin_epi_vmul_1xi64(v_q, v_coef_mod, max_gvl);
        v_result_1 = __builtin_epi_vsub_1xi64(v_mult, v_aux, max_gvl);
        mask = __builtin_epi_vmsleu_1xi64(v_coef_mod, v_result_1, max_gvl);
        v_result_1 = __builtin_epi_vsub_1xi64_mask(v_result_1, v_result_1, v_coef_mod, mask, max_gvl);

        __builtin_epi_vstore_unsigned_1xi64(&(output_stage_array[index]), v_result_0, max_gvl);
        __builtin_epi_vstore_unsigned_1xi64(&(output_stage_array[index + stopvalue_sft1]), v_result_1, max_gvl);
    }

    swap_pointer_aux = output_stage_array;
    output_stage_array = input_stage_array;
    input_stage_array = swap_pointer_aux;
    m >>= 1;
    element_or_aux ^= 1;

    for (; m >= gvl; m >>= 1, element_or_aux ^= 1)
    {
        for (uint32_t i = 0; i < (m >> log_of_gvl); i++) // m/256
        {
            uint64_t index = i << log_of_gvl;

            __epi_1xi64 v_S = __builtin_epi_vload_unsigned_1xi64(&rootOfUnityInverseTable[index + m], max_gvl);
            __epi_1xi64 v_precomp_aux = __builtin_epi_vload_unsigned_1xi64(&preconRootOfUnityInverseTable[index + m], max_gvl);

            for (uint32_t j = 0; j < (stopvalue_sft1 / m); j += 2)
            {
                __epi_1xi64 v_U = __builtin_epi_vload_strided_unsigned_1xi64(&(input_stage_array[((index + j * m) << 1)]), 16, max_gvl);
                __epi_1xi64 v_V = __builtin_epi_vload_strided_unsigned_1xi64(&(input_stage_array[((index + j * m) << 1) + 1]), 16, max_gvl);
                __epi_1xi64 _v_U = __builtin_epi_vload_strided_unsigned_1xi64(&(input_stage_array[((index + (j + 1) * m) << 1)]), 16, max_gvl);
                __epi_1xi64 _v_V = __builtin_epi_vload_strided_unsigned_1xi64(&(input_stage_array[((index + (j + 1) * m) << 1) + 1]), 16, max_gvl);

                __epi_1xi64 v_result_0 = __builtin_epi_vadd_1xi64(v_U, v_V, max_gvl);
                __epi_1xi64 _v_result_0 = __builtin_epi_vadd_1xi64(_v_U, _v_V, max_gvl);

                // Reduction v_result_0
                __epi_1xi1 mask = __builtin_epi_vmsleu_1xi64(v_coef_mod, v_result_0, max_gvl);
                v_result_0 = __builtin_epi_vsub_1xi64_mask(v_result_0, v_result_0, v_coef_mod, mask, max_gvl);
                __epi_1xi1 _mask = __builtin_epi_vmsleu_1xi64(v_coef_mod, _v_result_0, max_gvl);
                _v_result_0 = __builtin_epi_vsub_1xi64_mask(_v_result_0, _v_result_0, v_coef_mod, mask, max_gvl);

                // Reduction v_result_1
                __epi_1xi64 v_result_1 = __builtin_epi_vsub_1xi64(v_U, v_V, max_gvl); // masked as well? for speedup?
                mask = __builtin_epi_vmsleu_1xi64(v_coef_mod, v_result_1, max_gvl);
                v_result_1 = __builtin_epi_vadd_1xi64_mask(v_result_1, v_result_1, v_coef_mod, mask, max_gvl);
                __epi_1xi64 _v_result_1 = __builtin_epi_vsub_1xi64(_v_U, _v_V, max_gvl); // masked as well? for speedup?
                _mask = __builtin_epi_vmsleu_1xi64(v_coef_mod, _v_result_1, max_gvl);
                _v_result_1 = __builtin_epi_vadd_1xi64_mask(_v_result_1, _v_result_1, v_coef_mod, _mask, max_gvl);

                __epi_1xi64 v_q = __builtin_epi_vmulh_1xi64(v_result_1, v_precomp_aux, max_gvl);
                __epi_1xi64 _v_q = __builtin_epi_vmulh_1xi64(_v_result_1, v_precomp_aux, max_gvl);

                __epi_1xi64 v_mult = __builtin_epi_vmul_1xi64(v_result_1, v_S, max_gvl);
                __epi_1xi64 v_aux = __builtin_epi_vmul_1xi64(v_q, v_coef_mod, max_gvl);
                v_result_1 = __builtin_epi_vsub_1xi64(v_mult, v_aux, max_gvl);
                mask = __builtin_epi_vmsleu_1xi64(v_coef_mod, v_result_1, max_gvl);
                v_result_1 = __builtin_epi_vsub_1xi64_mask(v_result_1, v_result_1, v_coef_mod, mask, max_gvl);
                __epi_1xi64 _v_mult = __builtin_epi_vmul_1xi64(_v_result_1, v_S, max_gvl);
                __epi_1xi64 _v_aux = __builtin_epi_vmul_1xi64(_v_q, v_coef_mod, max_gvl);
                _v_result_1 = __builtin_epi_vsub_1xi64(_v_mult, _v_aux, max_gvl);
                _mask = __builtin_epi_vmsleu_1xi64(v_coef_mod, _v_result_1, max_gvl);
                _v_result_1 = __builtin_epi_vsub_1xi64_mask(_v_result_1, _v_result_1, v_coef_mod, _mask, max_gvl);

                __builtin_epi_vstore_unsigned_1xi64(&(output_stage_array[index + j * m]), v_result_0, max_gvl);
                __builtin_epi_vstore_unsigned_1xi64(&(output_stage_array[index + j * m + stopvalue_sft1]), v_result_1, max_gvl);
                __builtin_epi_vstore_unsigned_1xi64(&(output_stage_array[index + (j + 1) * m]), _v_result_0, max_gvl);
                __builtin_epi_vstore_unsigned_1xi64(&(output_stage_array[index + (j + 1) * m + stopvalue_sft1]), _v_result_1, max_gvl);
            }
        }

        swap_pointer_aux = output_stage_array;
        output_stage_array = input_stage_array;
        input_stage_array = swap_pointer_aux;
    }

    __epi_1xi64 v_omega_og = __builtin_epi_vload_unsigned_1xi64(&rootOfUnityInverseTable[0], max_gvl);
    __epi_1xi64 v_precomp_og = __builtin_epi_vload_unsigned_1xi64(&preconRootOfUnityInverseTable[0], max_gvl);

    for (; m >= 2; m >>= 1, element_or_aux ^= 1)
    {
        __epi_1xi64 v_index = __builtin_epi_vand_1xi64(v_index_original, v_index_stage_aux, gvl);
        v_index = __builtin_epi_vadd_1xi64(v_index, v_m, gvl);
        __epi_1xi64 v_S = __builtin_epi_vrgather_1xi64(v_omega_og, v_index, gvl);
        __epi_1xi64 v_precomp_aux = __builtin_epi_vrgather_1xi64(v_precomp_og, v_index, gvl);

        v_index_stage_aux = __builtin_epi_vsrl_1xi64(v_index_stage_aux, v_index_1, gvl);
        v_m = __builtin_epi_vsrl_1xi64(v_m, v_index_1, gvl);

        for (uint32_t i = 0; i < p; i += (max_gvl << 2))
        {
            uint32_t ii = i + (max_gvl << 1);

            __epi_1xi64 v_U = __builtin_epi_vload_strided_unsigned_1xi64(&(input_stage_array[i]), 16, max_gvl);
            __epi_1xi64 v_V = __builtin_epi_vload_strided_unsigned_1xi64(&(input_stage_array[i + 1]), 16, max_gvl);
            __epi_1xi64 _v_U = __builtin_epi_vload_strided_unsigned_1xi64(&(input_stage_array[ii]), 16, max_gvl);
            __epi_1xi64 _v_V = __builtin_epi_vload_strided_unsigned_1xi64(&(input_stage_array[ii + 1]), 16, max_gvl);

            __epi_1xi64 v_result_0 = __builtin_epi_vadd_1xi64(v_U, v_V, max_gvl);
            __epi_1xi64 _v_result_0 = __builtin_epi_vadd_1xi64(_v_U, _v_V, max_gvl);

            // Reduction v_result_0
            __epi_1xi1 mask = __builtin_epi_vmsleu_1xi64(v_coef_mod, v_result_0, max_gvl);
            v_result_0 = __builtin_epi_vsub_1xi64_mask(v_result_0, v_result_0, v_coef_mod, mask, max_gvl);
            __epi_1xi1 _mask = __builtin_epi_vmsleu_1xi64(v_coef_mod, _v_result_0, max_gvl);
            _v_result_0 = __builtin_epi_vsub_1xi64_mask(_v_result_0, _v_result_0, v_coef_mod, _mask, max_gvl);

            // Reduction v_result_1
            __epi_1xi64 v_result_1 = __builtin_epi_vsub_1xi64(v_U, v_V, max_gvl); // masked as well? for speedup?
            mask = __builtin_epi_vmsleu_1xi64(v_coef_mod, v_result_1, max_gvl);
            v_result_1 = __builtin_epi_vadd_1xi64_mask(v_result_1, v_result_1, v_coef_mod, mask, max_gvl);
            __epi_1xi64 _v_result_1 = __builtin_epi_vsub_1xi64(_v_U, _v_V, max_gvl); // masked as well? for speedup?
            _mask = __builtin_epi_vmsleu_1xi64(v_coef_mod, _v_result_1, max_gvl);
            _v_result_1 = __builtin_epi_vadd_1xi64_mask(_v_result_1, _v_result_1, v_coef_mod, _mask, max_gvl);

            __epi_1xi64 v_q = __builtin_epi_vmulh_1xi64(v_result_1, v_precomp_aux, max_gvl);
            __epi_1xi64 _v_q = __builtin_epi_vmulh_1xi64(_v_result_1, v_precomp_aux, max_gvl);

            __epi_1xi64 v_mult = __builtin_epi_vmul_1xi64(v_result_1, v_S, max_gvl);
            __epi_1xi64 v_aux = __builtin_epi_vmul_1xi64(v_q, v_coef_mod, max_gvl);
            v_result_1 = __builtin_epi_vsub_1xi64(v_mult, v_aux, max_gvl);
            mask = __builtin_epi_vmsleu_1xi64(v_coef_mod, v_result_1, max_gvl);
            v_result_1 = __builtin_epi_vsub_1xi64_mask(v_result_1, v_result_1, v_coef_mod, mask, max_gvl);

            __epi_1xi64 _v_mult = __builtin_epi_vmul_1xi64(_v_result_1, v_S, max_gvl);
            __epi_1xi64 _v_aux = __builtin_epi_vmul_1xi64(_v_q, v_coef_mod, max_gvl);
            _v_result_1 = __builtin_epi_vsub_1xi64(_v_mult, _v_aux, max_gvl);
            _mask = __builtin_epi_vmsleu_1xi64(v_coef_mod, _v_result_1, max_gvl);
            _v_result_1 = __builtin_epi_vsub_1xi64_mask(_v_result_1, _v_result_1, v_coef_mod, _mask, max_gvl);

            __builtin_epi_vstore_unsigned_1xi64(&(output_stage_array[i >> 1]), v_result_0, max_gvl);
            __builtin_epi_vstore_unsigned_1xi64(&(output_stage_array[(i >> 1) + stopvalue_sft1]), v_result_1, max_gvl);
            __builtin_epi_vstore_unsigned_1xi64(&(output_stage_array[ii >> 1]), _v_result_0, max_gvl);
            __builtin_epi_vstore_unsigned_1xi64(&(output_stage_array[(ii >> 1) + stopvalue_sft1]), _v_result_1, max_gvl);
        }

        swap_pointer_aux = output_stage_array;
        output_stage_array = input_stage_array;
        input_stage_array = swap_pointer_aux;
    }

    if (m != 0)
    {
        v_omega_og = __builtin_epi_vmv_v_x_1xi64(rootOfUnityInverseTable[1], max_gvl);
        v_precomp_og = __builtin_epi_vmv_v_x_1xi64(preconRootOfUnityInverseTable[1], max_gvl);

        for (uint32_t i = 0; i < p; i += (max_gvl << 2))
        {
            uint32_t ii = i + (max_gvl << 1);

            __epi_1xi64 v_U = __builtin_epi_vload_strided_unsigned_1xi64(&(input_stage_array[i]), 16, max_gvl);
            __epi_1xi64 v_V = __builtin_epi_vload_strided_unsigned_1xi64(&(input_stage_array[i + 1]), 16, max_gvl);
            __epi_1xi64 _v_U = __builtin_epi_vload_strided_unsigned_1xi64(&(input_stage_array[ii]), 16, max_gvl);
            __epi_1xi64 _v_V = __builtin_epi_vload_strided_unsigned_1xi64(&(input_stage_array[ii + 1]), 16, max_gvl);

            __epi_1xi64 v_result_0 = __builtin_epi_vadd_1xi64(v_U, v_V, max_gvl);
            __epi_1xi64 _v_result_0 = __builtin_epi_vadd_1xi64(_v_U, _v_V, max_gvl);

            // Reduction v_result_0
            __epi_1xi1 mask = __builtin_epi_vmsleu_1xi64(v_coef_mod, v_result_0, max_gvl);
            v_result_0 = __builtin_epi_vsub_1xi64_mask(v_result_0, v_result_0, v_coef_mod, mask, max_gvl);
            __epi_1xi1 _mask = __builtin_epi_vmsleu_1xi64(v_coef_mod, _v_result_0, max_gvl);
            _v_result_0 = __builtin_epi_vsub_1xi64_mask(_v_result_0, _v_result_0, v_coef_mod, _mask, max_gvl);

            // Reduction v_result_1
            __epi_1xi64 v_result_1 = __builtin_epi_vsub_1xi64(v_U, v_V, max_gvl); // masked as well? for speedup?
            mask = __builtin_epi_vmsleu_1xi64(v_coef_mod, v_result_1, max_gvl);
            v_result_1 = __builtin_epi_vadd_1xi64_mask(v_result_1, v_result_1, v_coef_mod, mask, max_gvl);
            __epi_1xi64 _v_result_1 = __builtin_epi_vsub_1xi64(_v_U, _v_V, max_gvl); // masked as well? for speedup?
            _mask = __builtin_epi_vmsleu_1xi64(v_coef_mod, _v_result_1, max_gvl);
            _v_result_1 = __builtin_epi_vadd_1xi64_mask(_v_result_1, _v_result_1, v_coef_mod, _mask, max_gvl);

            __epi_1xi64 v_q = __builtin_epi_vmulh_1xi64(v_result_1, v_precomp_og, max_gvl);
            __epi_1xi64 _v_q = __builtin_epi_vmulh_1xi64(_v_result_1, v_precomp_og, max_gvl);

            __epi_1xi64 v_mult = __builtin_epi_vmul_1xi64(v_result_1, v_omega_og, max_gvl);
            __epi_1xi64 v_aux = __builtin_epi_vmul_1xi64(v_q, v_coef_mod, max_gvl);
            v_result_1 = __builtin_epi_vsub_1xi64(v_mult, v_aux, max_gvl);
            mask = __builtin_epi_vmsleu_1xi64(v_coef_mod, v_result_1, max_gvl);
            v_result_1 = __builtin_epi_vsub_1xi64_mask(v_result_1, v_result_1, v_coef_mod, mask, max_gvl);
            __epi_1xi64 _v_mult = __builtin_epi_vmul_1xi64(_v_result_1, v_omega_og, max_gvl);
            __epi_1xi64 _v_aux = __builtin_epi_vmul_1xi64(_v_q, v_coef_mod, max_gvl);
            _v_result_1 = __builtin_epi_vsub_1xi64(_v_mult, _v_aux, max_gvl);
            _mask = __builtin_epi_vmsleu_1xi64(v_coef_mod, _v_result_1, max_gvl);
            _v_result_1 = __builtin_epi_vsub_1xi64_mask(_v_result_1, _v_result_1, v_coef_mod, _mask, max_gvl);

            __builtin_epi_vstore_unsigned_1xi64(&(output_stage_array[i >> 1]), v_result_0, max_gvl);
            __builtin_epi_vstore_unsigned_1xi64(&(output_stage_array[(i >> 1) + stopvalue_sft1]), v_result_1, max_gvl);
            __builtin_epi_vstore_unsigned_1xi64(&(output_stage_array[ii >> 1]), _v_result_0, max_gvl);
            __builtin_epi_vstore_unsigned_1xi64(&(output_stage_array[(ii >> 1) + stopvalue_sft1]), _v_result_1, max_gvl);
        }

        element_or_aux ^= 1;
    }

    // for half of the m values,  the results will end up in the auxiliary array rather than the element array
    if (element_or_aux)
    {
        for (uint32_t i = 0; i < p; i += max_gvl)
        {
            __epi_1xi64 v_auxarr = __builtin_epi_vload_unsigned_1xi64(&(output_stage_array[i]), max_gvl);
            __builtin_epi_vstore_unsigned_1xi64(&element[i], v_auxarr, max_gvl);
        }
    }
}