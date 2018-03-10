#ifndef _NEW_MURMURHASH3_H
#define _NEW_MURMURHASH3_H
#include <stdint.h>

//#include <immintrin.h>
#include <emmintrin.h>
#include <tmmintrin.h>
#include <immintrin.h>


//#define __ONMIC__ __attribute__((target(mic)))
#define __ONMIC__ 
__ONMIC__ void MurmurHash3_x86_32(const void* key, __m512i vkey_index, __mmask16 vcontinue, 
		uint32_t len, __m512i vhash_seed, __m512i& vhash_value, uint32_t increment);

__ONMIC__ void MurmurHash3_x86_32(const void* key, __m512i vkey_index, __mmask16 vcontinue, 
		uint32_t len, __m512i vhash_seed, __m512i& vhash_value, __m512i vrtail, __mmask16 vtail_mask);
__ONMIC__ void MurmurHash3_x86_32(const void* key, __m512i vkey_index, __mmask16 vcontinue, 
		uint32_t len, __m512i vhash_seed, __m512i& vhash_value, __m512i vrtail, __mmask16 vtail_mask, uint32_t pos);

__ONMIC__ void one_MurmurHash3_x86_32(const void* key, __m512i vkey_index, __mmask16 vcontinue, 
		uint32_t len, __m512i vhash_seed, __m512i& vhash_value);
#endif
