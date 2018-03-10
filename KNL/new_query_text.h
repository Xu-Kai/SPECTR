#ifndef _NEW_QUERY_TEXT_H_
#define _NEW_QUERY_TEXT_H_
#include "define.h"
//#include <immintrin.h>
#include <emmintrin.h>
#include <tmmintrin.h>
#include <immintrin.h>


#include "new_MurmurHash3.h"
#include <stdint.h>
//#define __ONMIC__ __attribute((target(mic)))
__ONMIC__ bool query_text(char* kmer, const uint8_t* bit_vector, const uint32_t* hash_seed, uint32_t kmer_length, uint32_t num_hash_func_real, uint32_t bit_vector_width);
__ONMIC__ void regions_query_text(char* sequence, const uint8_t* bit_vector, const uint32_t *hash_seed, uint32_t kmer_length, uint32_t num_hash_func, uint32_t bit_vector_width, uint32_t read_length, uint32_t *regions_text);
__ONMIC__ void path_query_text(char *kmer, const uint8_t *bit_vector, const uint32_t *hash_seed, uint32_t kmer_length, uint32_t num_hash_func, uint32_t bit_vector_width, uint32_t tail_flag, char trim_base, bool *pass_result);
__ONMIC__ void pos_path_query_text(char *kmer, const uint8_t *bit_vector, const uint32_t *hash_seed, uint32_t kmer_length, uint32_t num_hash_func, uint32_t bit_vector_width, uint32_t pos, bool *pass_result);
__ONMIC__ bool one_path_query_text(char *kmer, const uint8_t *bit_vector, const uint32_t *hash_seed, uint32_t kmer_length, uint32_t num_hash_func, uint32_t bit_vector_width);
#endif
