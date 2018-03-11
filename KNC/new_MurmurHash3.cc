#include "new_MurmurHash3.h"
#include <stdio.h>
__ONMIC__ void MurmurHash3_x86_32(const void* key, __m512i vkey_index, __mmask16 vcontinue, 
		uint32_t len, __m512i vhash_seed, __m512i& vhash_value, uint32_t increment){
	__m512i c1 = _mm512_set1_epi32(0xcc9e2d51);
	__m512i c2 = _mm512_set1_epi32(0x1b873593);
	__m512i four = _mm512_set1_epi32(increment);
	__m512i h1 = vhash_seed;
	int nblocks = len >> 2;
	//__m512i vnblocks = _mm512_set1_epi32(nblocks * 4);
	//vkey_index = _mm512_add_epi32(vkey_index, vnblocks);
	/*
	uint32_t *ttt = (uint32_t*)_mm_malloc(16 * 4, 64);
    printf("nblocks %d\n", nblocks);
    */
	for(int i = -nblocks; i; i ++){
		__m512i k1 = _mm512_mask_i32gather_epi32(four, vcontinue, vkey_index, key, 4);
	    /*	
		_mm512_store_epi32(ttt, k1);
		for(uint32_t ii = 0; ii < 16; ii ++){
			printf("%x ", ttt[ii]);	
		}
		printf(" i \n");
	    */	
		k1 = _mm512_mask_mullo_epi32(k1, vcontinue, k1, c1);
		__m512i tmp = _mm512_slli_epi32(k1, 15);
		k1 = _mm512_srli_epi32(k1, 32 - 15);
		k1 = _mm512_or_epi32(tmp, k1);
		k1 = _mm512_mask_mullo_epi32(k1, vcontinue, k1, c2);
		h1 = _mm512_xor_epi32(h1, k1);
		tmp = _mm512_slli_epi32(h1, 13);
		h1 = _mm512_srli_epi32(h1, 32 - 13);
		h1 = _mm512_or_epi32(h1, tmp);
		tmp = _mm512_set1_epi32(5);
		h1 = _mm512_mask_mullo_epi32(h1, vcontinue, h1, tmp);
		tmp = _mm512_set1_epi32(0xe6546b64);
		h1 = _mm512_add_epi32(h1, tmp);
		vkey_index = _mm512_add_epi32(vkey_index, four);	
	
	}
	
	__m512i k1 = _mm512_set1_epi32(0);
	//key_index = _mm512_slli_epi32(key_index, 2);
	__m512i vtail;
	//vkey_index = _mm512_add_epi32(vkey_index, vnblocks);
	/*
	switch(len & 3){
		__m512i vlast = _mm512_mask_i32gather_epi32(four, vcontinue, vkey_index, key, 4);
		case 3:
		vtail = _mm512_and_epi32(vlast, vtail_mask);
		vtail = _mm512_slli_epi32(vtail, 16);
		k1 = _mm512_xor_epi32(k1, vtail);
		vlast = _mm512_srli_epi32(vlast, 8);
		case 2:
		vtail = _mm512_and_epi32(vlast, vtail_mask);
		vtail = _mm512_slli_epi32(vtail, 8);
		k1 = _mm512_xor_epi32(k1, vtail);
		vlast = _mm512_srli_epi32(vlast, 8);
		case 1:
		vtail = _mm512_and_epi32(vlast, vtail_mask);
		k1 = _mm512_xor_epi32(k1, vtail);

		k1 = _mm512_mask_mullo_epi32(k1, vcontinue, k1, c1);
		__m512i tmp = _mm512_slli_epi32(k1, 15);
		k1 = _mm512_srli_epi32(k1, 32 - 15);
		k1 = _mm512_or_epi32(tmp, k1);	
		k1 = _mm512_mask_mullo_epi32(k1, vcontinue, k1, c2);
		h1 = _mm512_xor_epi32(h1, k1);
	}
	*/
		/*	
		_mm512_store_epi32(ttt, vkey_index);
		for(uint32_t i = 0; i < 16; i ++){
			printf("%d ", ttt[i]);	
		}
		printf("***\n");
		uint32_t *kkk = (uint32_t*)key;
		printf("%x \n", kkk[224]);
		*/
	if(len & 3){
		k1 = _mm512_mask_i32gather_epi32(four, vcontinue, vkey_index, key, 4);
		__m512i vtail_mask = _mm512_set1_epi32((1 << ((len & 3) << 3)) - 1);
		k1 = _mm512_and_epi32(vtail_mask, k1);

		k1 = _mm512_mask_mullo_epi32(k1, vcontinue, k1, c1);
		__m512i tmp = _mm512_slli_epi32(k1, 15);
		k1 = _mm512_srli_epi32(k1, 32 - 15);
		k1 = _mm512_or_epi32(tmp, k1);	
		k1 = _mm512_mask_mullo_epi32(k1, vcontinue, k1, c2);
		h1 = _mm512_xor_epi32(h1, k1);
		
	}


	k1 = _mm512_set1_epi32(len);
	h1 = _mm512_xor_epi32(h1, k1);
	c1 = _mm512_set1_epi32(0x85ebca6b);
	c2 = _mm512_set1_epi32(0xc2b2ae35);
	//    h ^= h >> 16;
	//    h *= 0x85ebca6b;
	//    h ^= h >> 13;
	//    h *= 0xc2b2ae35;
	//    h ^= h >> 16;
	k1 = _mm512_srli_epi32(h1, 16);
	h1 = _mm512_xor_epi32(h1, k1);
	h1 = _mm512_mask_mullo_epi32(h1, vcontinue, h1, c1);
	k1 =  _mm512_srli_epi32(h1, 13);
	h1 = _mm512_xor_epi32(h1, k1);
	h1 = _mm512_mask_mullo_epi32(h1, vcontinue, h1, c2);
	k1 = _mm512_srli_epi32(h1, 16);
	h1 = _mm512_xor_epi32(h1, k1);


	//    *(uint32_t*)out = h1;
	vhash_value = h1;


}
__ONMIC__ void one_MurmurHash3_x86_32(const void* key, __m512i vkey_index, __mmask16 vcontinue, 
		uint32_t len, __m512i vhash_seed, __m512i& vhash_value){
	__m512i c1 = _mm512_set1_epi32(0xcc9e2d51);
	__m512i c2 = _mm512_set1_epi32(0x1b873593);
	__m512i one = _mm512_set1_epi32(1);
	__m512i h1 = vhash_seed;
	int nblocks = len >> 2;
	for(int i = -nblocks; i; i ++){
		__m512i k1 = _mm512_mask_i32gather_epi32(one, vcontinue, vkey_index, key, 4);
	    /* 			
        int *tmp1 = (int*)_mm_malloc(16 * sizeof(int), 64);
        _mm512_store_epi32(tmp1, k1);
        for(int i = 0; i < 16; i ++)
            printf("hh %x ", tmp1[i]);
        printf("\n");
        */ 
		k1 = _mm512_mask_mullo_epi32(k1, vcontinue, k1, c1);
		__m512i tmp = _mm512_slli_epi32(k1, 15);
		k1 = _mm512_srli_epi32(k1, 32 - 15);
		k1 = _mm512_or_epi32(tmp, k1);
		k1 = _mm512_mask_mullo_epi32(k1, vcontinue, k1, c2);
		h1 = _mm512_xor_epi32(h1, k1);
		tmp = _mm512_slli_epi32(h1, 13);
		h1 = _mm512_srli_epi32(h1, 32 - 13);
		h1 = _mm512_or_epi32(h1, tmp);
		tmp = _mm512_set1_epi32(5);
		h1 = _mm512_mask_mullo_epi32(h1, vcontinue, h1, tmp);
		tmp = _mm512_set1_epi32(0xe6546b64);
		h1 = _mm512_add_epi32(h1, tmp);
		vkey_index = _mm512_add_epi32(vkey_index, one);	
	
	}
	
	__m512i k1 = _mm512_set1_epi32(0);
	if(len & 3){
		k1 = _mm512_mask_i32gather_epi32(one, vcontinue, vkey_index, key, 4);


		__m512i vtail_mask = _mm512_set1_epi32((1 << ((len & 3) << 3)) - 1);
		k1 = _mm512_and_epi32(vtail_mask, k1);

		k1 = _mm512_mask_mullo_epi32(k1, vcontinue, k1, c1);
		__m512i tmp = _mm512_slli_epi32(k1, 15);
		k1 = _mm512_srli_epi32(k1, 32 - 15);
		k1 = _mm512_or_epi32(tmp, k1);	
		k1 = _mm512_mask_mullo_epi32(k1, vcontinue, k1, c2);
		h1 = _mm512_xor_epi32(h1, k1);
		
	}


	k1 = _mm512_set1_epi32(len);
	h1 = _mm512_xor_epi32(h1, k1);
	c1 = _mm512_set1_epi32(0x85ebca6b);
	c2 = _mm512_set1_epi32(0xc2b2ae35);
	//    h ^= h >> 16;
	//    h *= 0x85ebca6b;
	//    h ^= h >> 13;
	//    h *= 0xc2b2ae35;
	//    h ^= h >> 16;
	k1 = _mm512_srli_epi32(h1, 16);
	h1 = _mm512_xor_epi32(h1, k1);
	h1 = _mm512_mask_mullo_epi32(h1, vcontinue, h1, c1);
	k1 =  _mm512_srli_epi32(h1, 13);
	h1 = _mm512_xor_epi32(h1, k1);
	h1 = _mm512_mask_mullo_epi32(h1, vcontinue, h1, c2);
	k1 = _mm512_srli_epi32(h1, 16);
	h1 = _mm512_xor_epi32(h1, k1);


	//    *(uint32_t*)out = h1;
	vhash_value = h1;


}

__ONMIC__ void MurmurHash3_x86_32(const void* key, __m512i vkey_index, __mmask16 vcontinue, 
		uint32_t len, __m512i vhash_seed, __m512i& vhash_value, __m512i vrtail, __mmask16 vtail_mask){
	__m512i c1 = _mm512_set1_epi32(0xcc9e2d51);
	__m512i c2 = _mm512_set1_epi32(0x1b873593);
	__m512i one = _mm512_set1_epi32(1);
	__m512i h1 = vhash_seed;
	int nblocks = len >> 2;
	//__m512i vnblocks = _mm512_set1_epi32(nblocks * 4);
	//vkey_index = _mm512_add_epi32(vkey_index, vnblocks);
	//uint32_t *ttt = (uint32_t*)_mm_malloc(16 * 4, 64);
	for(int i = -nblocks; i; i ++){
		__m512i k1 = _mm512_mask_i32gather_epi32(one, vcontinue, vkey_index, key, 4);
		/*
		_mm512_store_epi32(ttt, k1);
		for(uint32_t j = 0; j < 16; j ++){
			printf("%x ", ttt[j]);	
		}
		printf("\nkkkkkk***\n");
		*/
		if(i == -nblocks){
			k1 = _mm512_mask_or_epi32(k1, vtail_mask, k1, vrtail);
		}
		
		
		k1 = _mm512_mask_mullo_epi32(k1, vcontinue, k1, c1);
		__m512i tmp = _mm512_slli_epi32(k1, 15);
		k1 = _mm512_srli_epi32(k1, 32 - 15);
		k1 = _mm512_or_epi32(tmp, k1);
		k1 = _mm512_mask_mullo_epi32(k1, vcontinue, k1, c2);
		h1 = _mm512_xor_epi32(h1, k1);
		tmp = _mm512_slli_epi32(h1, 13);
		h1 = _mm512_srli_epi32(h1, 32 - 13);
		h1 = _mm512_or_epi32(h1, tmp);
		tmp = _mm512_set1_epi32(5);
		h1 = _mm512_mask_mullo_epi32(h1, vcontinue, h1, tmp);
		tmp = _mm512_set1_epi32(0xe6546b64);
		h1 = _mm512_add_epi32(h1, tmp);
		vkey_index = _mm512_add_epi32(vkey_index, one);	
	
	}
	
	__m512i k1 = _mm512_set1_epi32(0);
	//key_index = _mm512_slli_epi32(key_index, 2);
	__m512i vtail;
	//vkey_index = _mm512_add_epi32(vkey_index, vnblocks);
	/*
	switch(len & 3){
		__m512i vlast = _mm512_mask_i32gather_epi32(four, vcontinue, vkey_index, key, 4);
		case 3:
		vtail = _mm512_and_epi32(vlast, vtail_mask);
		vtail = _mm512_slli_epi32(vtail, 16);
		k1 = _mm512_xor_epi32(k1, vtail);
		vlast = _mm512_srli_epi32(vlast, 8);
		case 2:
		vtail = _mm512_and_epi32(vlast, vtail_mask);
		vtail = _mm512_slli_epi32(vtail, 8);
		k1 = _mm512_xor_epi32(k1, vtail);
		vlast = _mm512_srli_epi32(vlast, 8);
		case 1:
		vtail = _mm512_and_epi32(vlast, vtail_mask);
		k1 = _mm512_xor_epi32(k1, vtail);

		k1 = _mm512_mask_mullo_epi32(k1, vcontinue, k1, c1);
		__m512i tmp = _mm512_slli_epi32(k1, 15);
		k1 = _mm512_srli_epi32(k1, 32 - 15);
		k1 = _mm512_or_epi32(tmp, k1);	
		k1 = _mm512_mask_mullo_epi32(k1, vcontinue, k1, c2);
		h1 = _mm512_xor_epi32(h1, k1);
	}
	*/
		/*	
		_mm512_store_epi32(ttt, vkey_index);
		for(uint32_t i = 0; i < 16; i ++){
			printf("%d ", ttt[i]);	
		}
		printf("***\n");
		uint32_t *kkk = (uint32_t*)key;
		printf("%x \n", kkk[224]);
		*/
	vtail_mask = _mm512_knot(vtail_mask);
	if(len & 3){
		k1 = _mm512_mask_i32gather_epi32(one, vcontinue, vkey_index, key, 4);

		k1 = _mm512_mask_or_epi32(k1, vtail_mask, k1, vrtail);

		__m512i vtail_mask = _mm512_set1_epi32((1 << ((len & 3) << 3)) - 1);
		k1 = _mm512_and_epi32(vtail_mask, k1);

		k1 = _mm512_mask_mullo_epi32(k1, vcontinue, k1, c1);
		__m512i tmp = _mm512_slli_epi32(k1, 15);
		k1 = _mm512_srli_epi32(k1, 32 - 15);
		k1 = _mm512_or_epi32(tmp, k1);	
		k1 = _mm512_mask_mullo_epi32(k1, vcontinue, k1, c2);
		h1 = _mm512_xor_epi32(h1, k1);
		
	}


	k1 = _mm512_set1_epi32(len);
	h1 = _mm512_xor_epi32(h1, k1);
	c1 = _mm512_set1_epi32(0x85ebca6b);
	c2 = _mm512_set1_epi32(0xc2b2ae35);
	//    h ^= h >> 16;
	//    h *= 0x85ebca6b;
	//    h ^= h >> 13;
	//    h *= 0xc2b2ae35;
	//    h ^= h >> 16;
	k1 = _mm512_srli_epi32(h1, 16);
	h1 = _mm512_xor_epi32(h1, k1);
	h1 = _mm512_mask_mullo_epi32(h1, vcontinue, h1, c1);
	k1 =  _mm512_srli_epi32(h1, 13);
	h1 = _mm512_xor_epi32(h1, k1);
	h1 = _mm512_mask_mullo_epi32(h1, vcontinue, h1, c2);
	k1 = _mm512_srli_epi32(h1, 16);
	h1 = _mm512_xor_epi32(h1, k1);


	//    *(uint32_t*)out = h1;
	vhash_value = h1;


}

__ONMIC__ void MurmurHash3_x86_32(const void* key, __m512i vkey_index, __mmask16 vcontinue, 
		uint32_t len, __m512i vhash_seed, __m512i& vhash_value, __m512i vrtail, __mmask16 vtail_mask, uint32_t pos){
	__m512i c1 = _mm512_set1_epi32(0xcc9e2d51);
	__m512i c2 = _mm512_set1_epi32(0x1b873593);
	__m512i one = _mm512_set1_epi32(1);
	__m512i h1 = vhash_seed;
	int nblocks = len >> 2;
    __mmask16 vntail_mask = _mm512_knot(vtail_mask);
	//__m512i vnblocks = _mm512_set1_epi32(nblocks * 4);
	//vkey_index = _mm512_add_epi32(vkey_index, vnblocks);
	//uint32_t *ttt = (uint32_t*)_mm_malloc(16 * 4, 64);
	for(int i = 0; i < nblocks; i ++){
		__m512i k1 = _mm512_mask_i32gather_epi32(one, vcontinue, vkey_index, key, 4);
	    /*   	
		_mm512_store_epi32(ttt, k1);
		for(uint32_t j = 0; j < 16; j ++){
			printf("%x ", ttt[j]);	
		}
		printf("\n");
        */
        
	
        /*
        int *tmp1 = (int*)_mm_malloc(16 * sizeof(int), 64);
        _mm512_store_epi32(tmp1, k1);
        printf("\n");
        printf("hh %x ", tmp1[0]);
        printf("\n");
        */
		if(pos  / 4  == i){
			k1 = _mm512_mask_or_epi32(k1, vntail_mask, k1, vrtail);
		}
		
		if((len - 1 - pos) / 4  == i){
			k1 = _mm512_mask_or_epi32(k1, vtail_mask, k1, vrtail);
		}

        /* 
	    int *tmp1 = (int*)_mm_malloc(16 * sizeof(int), 64);
        _mm512_store_epi32(tmp1, k1);
        for(int j = 0; j < 16; j ++)
        printf("%x ", tmp1[j]);
        printf("\n");
       */ 
	
		k1 = _mm512_mask_mullo_epi32(k1, vcontinue, k1, c1);
		__m512i tmp = _mm512_slli_epi32(k1, 15);
		k1 = _mm512_srli_epi32(k1, 32 - 15);
		k1 = _mm512_or_epi32(tmp, k1);
		k1 = _mm512_mask_mullo_epi32(k1, vcontinue, k1, c2);
		h1 = _mm512_xor_epi32(h1, k1);
		tmp = _mm512_slli_epi32(h1, 13);
		h1 = _mm512_srli_epi32(h1, 32 - 13);
		h1 = _mm512_or_epi32(h1, tmp);
		tmp = _mm512_set1_epi32(5);
		h1 = _mm512_mask_mullo_epi32(h1, vcontinue, h1, tmp);
		tmp = _mm512_set1_epi32(0xe6546b64);
		h1 = _mm512_add_epi32(h1, tmp);
		vkey_index = _mm512_add_epi32(vkey_index, one);	
	
	}
	
	__m512i k1 = _mm512_set1_epi32(0);
	//key_index = _mm512_slli_epi32(key_index, 2);
	__m512i vtail;
	//vkey_index = _mm512_add_epi32(vkey_index, vnblocks);
	/*
	switch(len & 3){
		__m512i vlast = _mm512_mask_i32gather_epi32(four, vcontinue, vkey_index, key, 4);
		case 3:
		vtail = _mm512_and_epi32(vlast, vtail_mask);
		vtail = _mm512_slli_epi32(vtail, 16);
		k1 = _mm512_xor_epi32(k1, vtail);
		vlast = _mm512_srli_epi32(vlast, 8);
		case 2:
		vtail = _mm512_and_epi32(vlast, vtail_mask);
		vtail = _mm512_slli_epi32(vtail, 8);
		k1 = _mm512_xor_epi32(k1, vtail);
		vlast = _mm512_srli_epi32(vlast, 8);
		case 1:
		vtail = _mm512_and_epi32(vlast, vtail_mask);
		k1 = _mm512_xor_epi32(k1, vtail);

		k1 = _mm512_mask_mullo_epi32(k1, vcontinue, k1, c1);
		__m512i tmp = _mm512_slli_epi32(k1, 15);
		k1 = _mm512_srli_epi32(k1, 32 - 15);
		k1 = _mm512_or_epi32(tmp, k1);	
		k1 = _mm512_mask_mullo_epi32(k1, vcontinue, k1, c2);
		h1 = _mm512_xor_epi32(h1, k1);
	}
	*/
		/*	
		_mm512_store_epi32(ttt, vkey_index);
		for(uint32_t i = 0; i < 16; i ++){
			printf("%d ", ttt[i]);	
		}
		printf("***\n");
		uint32_t *kkk = (uint32_t*)key;
		printf("%x \n", kkk[224]);
		*/
	if(len & 3){
		k1 = _mm512_mask_i32gather_epi32(one, vcontinue, vkey_index, key, 4);
        /* 
	    int *tmp1 = (int*)_mm_malloc(16 * sizeof(int), 64);
        _mm512_store_epi32(tmp1, k1);
        for(int i = 0; i < 16; i ++)
        printf("%x ", tmp1[i]);
        printf("\n");
        */

       /* _mm512_store_epi32(tmp1, vrtail);
        for(int i = 0; i < 16; i ++)
        printf("%x ", tmp1[i]);
        printf("\n");
    */
		if(pos >= nblocks * 4){
			k1 = _mm512_mask_or_epi32(k1, vntail_mask, k1, vrtail);
		}
		
		if((len - 1 - pos) >= nblocks * 4){
			k1 = _mm512_mask_or_epi32(k1, vtail_mask, k1, vrtail);
		}
		vtail = _mm512_set1_epi32((1 << ((len & 3) << 3)) - 1);
		k1 = _mm512_and_epi32(vtail, k1);
    
        /* 
        _mm512_store_epi32(tmp1, k1);
        for(int i = 0; i < 16; i ++)
        printf("%x ", tmp1[i]);
        printf("\n");
        */
    
		k1 = _mm512_mask_mullo_epi32(k1, vcontinue, k1, c1);
		__m512i tmp = _mm512_slli_epi32(k1, 15);
		k1 = _mm512_srli_epi32(k1, 32 - 15);
		k1 = _mm512_or_epi32(tmp, k1);	
		k1 = _mm512_mask_mullo_epi32(k1, vcontinue, k1, c2);
		h1 = _mm512_xor_epi32(h1, k1);
		
	}


	k1 = _mm512_set1_epi32(len);
	h1 = _mm512_xor_epi32(h1, k1);
	c1 = _mm512_set1_epi32(0x85ebca6b);
	c2 = _mm512_set1_epi32(0xc2b2ae35);
	//    h ^= h >> 16;
	//    h *= 0x85ebca6b;
	//    h ^= h >> 13;
	//    h *= 0xc2b2ae35;
	//    h ^= h >> 16;
	k1 = _mm512_srli_epi32(h1, 16);
	h1 = _mm512_xor_epi32(h1, k1);
	h1 = _mm512_mask_mullo_epi32(h1, vcontinue, h1, c1);
	k1 =  _mm512_srli_epi32(h1, 13);
	h1 = _mm512_xor_epi32(h1, k1);
	h1 = _mm512_mask_mullo_epi32(h1, vcontinue, h1, c2);
	k1 = _mm512_srli_epi32(h1, 16);
	h1 = _mm512_xor_epi32(h1, k1);


	//    *(uint32_t*)out = h1;
	vhash_value = h1;


}

