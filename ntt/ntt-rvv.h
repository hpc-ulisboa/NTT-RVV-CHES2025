#ifndef NTT_RVV
#define NTT_RVV

#include <stdint.h>

// Just so a static analyzer does not report a ton of errors
#ifdef __INTELLISENSE__

#define __epi_e64 1
#define __epi_m1 1

typedef int __epi_1xi64;
typedef int __epi_1xi1;

inline int __builtin_epi_vsetvl(unsigned, int, int) { return 256; }
inline __epi_1xi64 __builtin_epi_vmv_v_x_1xi64(long, int) { return 0; }
inline __epi_1xi64 __builtin_epi_vload_unsigned_1xi64(const uint64_t *, int) { return 0; }
inline __epi_1xi64 __builtin_epi_vload_strided_unsigned_1xi64(const uint64_t *, int, int) { return 0; }
inline void __builtin_epi_vstore_unsigned_1xi64(uint64_t *, __epi_1xi64, int) {}
inline void __builtin_epi_vstore_strided_unsigned_1xi64(uint64_t *, __epi_1xi64, int, int) {}
inline void __builtin_epi_vstore_unsigned_1xi64_mask(const uint64_t *, __epi_1xi64, __epi_1xi1, int);

inline __epi_1xi64 __builtin_epi_vadd_1xi64(__epi_1xi64, __epi_1xi64, int) { return 0; }
inline __epi_1xi64 __builtin_epi_vsub_1xi64(__epi_1xi64, __epi_1xi64, int) { return 0; }
inline __epi_1xi64 __builtin_epi_vmul_1xi64(__epi_1xi64, __epi_1xi64, int) { return 0; }
inline __epi_1xi64 __builtin_epi_vmulh_1xi64(__epi_1xi64, __epi_1xi64, int) { return 0; }
inline __epi_1xi64 __builtin_epi_vand_1xi64(__epi_1xi64, __epi_1xi64, int) { return 0; }
inline __epi_1xi64 __builtin_epi_vor_1xi64(__epi_1xi64, __epi_1xi64, int) { return 0; }
inline __epi_1xi64 __builtin_epi_vsll_1xi64(__epi_1xi64, __epi_1xi64, int) { return 0; }
inline __epi_1xi64 __builtin_epi_vsrl_1xi64(__epi_1xi64, __epi_1xi64, int) { return 0; }
inline __epi_1xi64 __builtin_epi_vrgather_1xi64(__epi_1xi64, __epi_1xi64, int) { return 0; }

inline __epi_1xi1 __builtin_epi_vmsleu_1xi64(__epi_1xi64, __epi_1xi64, int) { return 0; }
inline __epi_1xi1 __builtin_epi_vload_1xi1(const uint64_t *) { return 0; }

inline __epi_1xi64 __builtin_epi_vsub_1xi64_mask(__epi_1xi64, __epi_1xi64, __epi_1xi64, __epi_1xi1, int) { return 0; }
inline __epi_1xi64 __builtin_epi_vadd_1xi64_mask(__epi_1xi64, __epi_1xi64, __epi_1xi64, __epi_1xi1, int) { return 0; }

#endif

void ntt_korn_lambiote_vector(
    const uint32_t p, uint32_t t, uint32_t logt, const uint64_t modulus, uint64_t *element,
    const uint64_t *rootOfUnityTable, const uint64_t *preconRootOfUnityTable
);

void intt_pease_vector_mulh(
    const uint32_t p, const uint64_t modulus, uint64_t *element, const uint64_t *rootOfUnityInverseTable,
    const uint64_t *preconRootOfUnityInverseTable, const uint64_t cycloOrderInv, const uint64_t preconCycloOrderInv
);

void free_ntts_mem();

#endif