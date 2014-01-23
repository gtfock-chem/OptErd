#ifndef __ERD_PROFILE_H__
#define __ERD_PROFILE_H__

#include <stdint.h>


#define MAXTHREADS     240
#define NUM_TICKS      21


typedef enum
{
    // general case
    erd__prepare_ctr_ticks              = 0,
    erd__set_ij_kl_pairs_ticks          = 1,
    // VRR
    erd__rys_roots_weights_ticks        = 2,
    erd__2d_coefficients_ticks          = 3,
    erd__2d_pq_integrals_ticks          = 4,
    erd__int2d_to_e0f0_ticks            = 5,
    erd__e0f0_pcgto_block_ticks         = 6,
    
    erd__xyz_to_ry_abcd_ticks           = 7,
    erd__hrr_matrix_ticks               = 8,
    erd__hrr_transform_ticks            = 9,
    erd__move_ry_ticks                  = 10,
    erd__spherical_transform_ticks      = 11,
    erd__csgto_ticks                    = 12,

    // 1111 case
    erd__prepare_ctr_ticks_1111         = 13,
    erd__set_ij_kl_pairs_ticks_1111     = 14,
    erd__ssss_pcgto_block_ticks         = 15,
    erd__sssp_pcgto_block_ticks         = 16,
    erd__sspp_pcgto_block_ticks         = 17,
    erd__sppp_pcgto_block_ticks         = 18,
    erd__pppp_pcgto_block_ticks         = 19,
    erd__1111_csgto_ticks               = 20
} ErdTicks_t;


extern __declspec(align(256)) uint64_t erd_ticks[MAXTHREADS][NUM_TICKS + 8];


void erd_reset_profile (void);

void erd_print_profile (int mode);


#endif /* __ERD_PROFILE_H__ */
