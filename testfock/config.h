#ifndef __CONFIG_H__
#define __CONFIG_H__


#if defined (__MIC__)
#define SIMDW      64
#elif defined (__AVX__)
#define SIMDW      32
#elif defined (__SSE__)
#define SIMDW      16
#else
#define SIMDW      64
#endif

#define MIN(a, b)    ((a) < (b) ? (a) : (b))
#define MAX(a, b)    ((a) > (b) ? (a) : (b))

#define _DEBUG_LEVEL_    2   //  0 to 10,
                             //  0 is no debug print info at all,
                             //  10 is full info


#define PFOCK_MALLOC(size,alignment)    _mm_malloc(size, alignment)
#define PFOCK_FREE(addr)    _mm_free(addr)

#if ( _DEBUG_LEVEL_ == -1 )
#define PFOCK_PRINTF( level, fmt, args... )        {}
#else
#define PFOCK_PRINTF( level, fmt, args... )                                          \
        do                                                                           \
        {                                                                            \
            if ( (unsigned)(level) <= _DEBUG_LEVEL_ )                                \
            {                                                                        \
                sprintf( pfock->str_buf, "%s() line %d: ", __FUNCTION__, __LINE__);  \
                sprintf( pfock->str_buf + strlen(pfock->str_buf), fmt, ##args );     \
                fprintf( stdout, "%s", pfock->str_buf );                             \
                fflush( stdout );                                                    \
            }                                                                        \
        } while ( 0 )
#endif


#if ( _DEBUG_LEVEL_ > 1 )
#define PFOCK_INFO( fmt, args... )                                             \
        do                                                                     \
        {                                                                      \
            sprintf( pfock->str_buf, "    PFock: ");                         \
            sprintf( pfock->str_buf + strlen("    PFock: "), fmt, ##args );  \
            fprintf( stdout, "%s", pfock->str_buf );                           \
            fflush( stdout );                                                  \
        } while ( 0 )
#else
#define PFOCK_INFO( fmt, args... )        {}
#endif


#if defined(_NO_COMM_)
#define MY_GET(ga, lo, hi, ptr, ld)
#define MY_ACC(ga, lo, hi, ptr, ld, one)
#elif defined(_NO_BLOCK_)
#define MY_GET(ga, lo, hi, ptr, ld)
#define MY_ACC(ga, lo, hi, ptr, ld, one)
#else
#define MY_GET(ga, lo, hi, ptr, ld)        NGA_Get(ga, lo, hi, ptr, ld)
#define MY_ACC(ga, lo, hi, ptr, ld, one)   NGA_Acc(ga, lo, hi, ptr, ld, one)
#endif


#endif /* __CONFIG_H__ */
