#include "new_query_text.h"
#include "MurmurHash3.h"
//#include <immintrin.h>
#define READ_MAX_LEN 256
#define KMER_MAX_LEN 128
//----------------------------------------------------------------------
// query_text
//----------------------------------------------------------------------
__ONMIC__ bool query_text(char* kmer, const uint8_t* bit_vector,
		const uint32_t* hash_seed, uint32_t kmer_length,
		uint32_t num_hash_func_real, uint32_t bit_vector_width) {
	//--------------------------------------------------
	// reverse complement or not
	//--------------------------------------------------

	//char* kmer_internal = (char*) malloc(kmer_length + 2);
	// lexicographical comparison
	//  1: kmer is faster than its reverse complement in the lexicographical order
	//  0: same
	// -1: reverse complement of kmer is faster than kmer in the lexicographical order
	char kmer_internal[KMER_MAX_LEN] ;

	int comparison_result = 0;
	for (unsigned int it = 0, it_rc = (kmer_length - 1); it < kmer_length;
			it++, it_rc--) {
		if (kmer[it] == 'A') {
			if (kmer[it_rc] == 'T') {
			} else {
				comparison_result = 1;
				break;
			}
		} else if (kmer[it] == 'C') {
			if (kmer[it_rc] == 'G') {
			} else if (kmer[it_rc] == 'T') {
				comparison_result = -1;
				break;
			} else {
				comparison_result = 1;
				break;
			}
		} else if (kmer[it] == 'G') {
			if (kmer[it_rc] == 'A') {
				comparison_result = 1;
				break;
			} else if (kmer[it_rc] == 'C') {
			} else {
				comparison_result = -1;
				break;
			}
		} else if (kmer[it] == 'T') {
			if (kmer[it_rc] == 'A') {
			} else {
				comparison_result = -1;
				break;
			}
		} else {
			//std::cout << std::endl << "ERROR: Illegal character " << kmer[it] << " (query_text)" << std::endl << std::endl;
			printf("ERROR: 1 Illegal charactor %c query_text\n\n", kmer[it]);
		}
	}
	//printf("com %d\n", comparison_result);
	//printf("result %d\n", comparison_result);
	// kmer is faster
	if (comparison_result >= 0) {
		//	      kmer_internal = kmer;
		memcpy(kmer_internal, kmer, kmer_length);
		//printf("hhhhh1\n");
	}
	// reverse complement is faster
	else {
		// reverse complement
		//printf("hhhhh2\n");
		for (unsigned int it = 0, rc_index = kmer_length - 1; it < kmer_length;
				it++, rc_index--) {
			switch (kmer[it]) {
				case 'A':
					kmer_internal[rc_index] = 'T';
					break;
				case 'C':
					kmer_internal[rc_index] = 'G';
					break;
				case 'G':
					kmer_internal[rc_index] = 'C';
					break;
				case 'T':
					kmer_internal[rc_index] = 'A';
					break;
				default:
					//   std::cout << std::endl << "ERROR: Illegal character " << kmer[it] << " (query_text)" << std::endl << std::endl;
					printf("ERROR: 2 Illegal charactor %c query_text\n\n", kmer[it]);
					break;
			}
		}
	}

	//	free(kmer_internal);

	//	   uint64_t original_index1;
	//	   uint64_t original_index2;
	//
	//	   uint64_t hash [2];
	//	 	   char* kmer_internal = kmer;
	uint32_t original_index1;
	//	uint32_t original_index2;
	//	uint32_t original_index3;
	//	uint32_t original_index4;
	uint32_t hash[4];
	//
	unsigned short int bit_index1;
	//	unsigned short int bit_index2;
	//	unsigned short int bit_index3;
	//	unsigned short int bit_index4;
	
	//   for(int i = 0; i < kmer_length; i ++){
	//   printf("%c", kmer_internal[i]);
	//   }
	//   printf("\n");
	   
	//uint32_t *cs = (uint32_t*)kmer_internal;
	//for(int i = 0; i < 8; i ++){
	//printf("%x ", cs[i]);
	//}
	//printf("hhhhhh\n");
	uint32_t or_index; 
	for (uint32_t it_hash_func = 0; it_hash_func < num_hash_func_real; it_hash_func++) {

		MurmurHash3_x86_32(kmer_internal, kmer_length, hash_seed[it_hash_func],
				&or_index);

		//printf("%x\n", or_index);
		original_index1 = or_index % bit_vector_width;
		bit_index1 = original_index1 % BITS_PER_CHAR;

		//printf("%x %x %x\n", original_index1, bit_vector[original_index1 / BITS_PER_CHAR], (bit_vector[original_index1 / BITS_PER_CHAR] & BIT_MASK[bit_index1]));
		//printf("original_index1 %u\n", original_index1);
		if ((bit_vector[original_index1 / BITS_PER_CHAR] & BIT_MASK[bit_index1])
				!= BIT_MASK[bit_index1]) {
			//_mm_free(kmer_internal);
			//printf("******\n");
			return false;
		}
	}

	//printf("******\n");
	//_mm_free(kmer_internal);
	//kmer[0] = 1;
	return true;

}
__ONMIC__ void regions_query_text(char *sequence, const uint8_t *bit_vector, const uint32_t *hash_seed, uint32_t kmer_length, uint32_t num_hash_func, uint32_t bit_vector_width, uint32_t read_length, uint32_t *regions_text){

	//#ifdef __MIC__
	uint32_t *regions_reverse_text = regions_text + READ_MAX_LEN;
	__m512i vseq, vtext, vseq_reverse, vtext_reverse;
	vtext                    = _mm512_set1_epi32(0);
	__m512i vseq_ind         = _mm512_set_epi32(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
	__m512i vseq_reverse_ind = _mm512_set1_epi32(read_length - 1);
	vseq_reverse_ind         = _mm512_sub_epi32(vseq_reverse_ind, vseq_ind);
	__m512i vincrement       = _mm512_set1_epi32(16);
	__m512i vreverse_char    = _mm512_set_epi32(0, 0, 0, 0, 0, 0, 0, 0, 0x43, 0, 0, 0x41, 0x47, 0, 0x54, 0);
	uint32_t iteration       = read_length / 16;
	__m512i vone              = _mm512_set1_epi32(1);
	__m512i vfour = _mm512_set1_epi32(4);
	__m512i vthree = _mm512_set1_epi32(3);
	__m512i chbit = _mm512_set1_epi32(0xff);
	for(uint32_t i = 0; i < iteration; i ++){
		for(uint32_t j = 0; j < 4; j ++){
			//vseq             = _mm512_i32extgather_epi32(vseq_ind, sequence, _MM_UPCONV_EPI32_UINT8, 1, 1);
			/*******************************************/
			__m512i vseq_ind32 = _mm512_srli_epi32(vseq_ind, 2);
			vseq             = _mm512_i32gather_epi32(vseq_ind32, sequence, 4);
			__m512i vseq_ind8 = _mm512_and_epi32(vseq_ind, vthree);
			vseq_ind8 = _mm512_slli_epi32(vseq_ind8, 3);
			vseq 		= _mm512_srlv_epi32(vseq, vseq_ind8);
			vseq		= _mm512_and_epi32(vseq, chbit);
			/******************************************/

			//vseq_reverse     = _mm512_i32extgather_epi32(vseq_reverse_ind, sequence, _MM_UPCONV_EPI32_UINT8, 1, 1);
			/********************************************************/
			vseq_ind32 = _mm512_srli_epi32(vseq_reverse_ind, 2);
			vseq_reverse             = _mm512_i32gather_epi32(vseq_ind32, sequence, 4);
			vseq_ind8 = _mm512_and_epi32(vseq_reverse_ind, vthree);
			vseq_ind8 = _mm512_slli_epi32(vseq_ind8, 3);
			vseq_reverse 	= _mm512_srlv_epi32(vseq_reverse, vseq_ind8);
			vseq_reverse	= _mm512_and_epi32(vseq_reverse, chbit);
			/********************************************************/

			vseq_reverse     = _mm512_permutevar_epi32(vseq_reverse, vreverse_char);
			vseq             = _mm512_slli_epi32(vseq, 24);
			vtext            = _mm512_srli_epi32(vtext, 8);
			vtext            = _mm512_or_epi32(vseq, vtext);
			vseq_reverse     = _mm512_slli_epi32(vseq_reverse, 24);
			vtext_reverse    = _mm512_srli_epi32(vtext_reverse, 8);
			vtext_reverse    = _mm512_or_epi32(vseq_reverse, vtext_reverse);
			vseq_ind         = _mm512_add_epi32(vseq_ind, vone);
			vseq_reverse_ind = _mm512_sub_epi32(vseq_reverse_ind, vone);
		}
		vseq_ind = _mm512_sub_epi32(vseq_ind, vfour);
		vseq_reverse_ind = _mm512_add_epi32(vseq_reverse_ind, vfour);
		_mm512_store_epi32(regions_text + (i << 4), vtext);
		_mm512_store_epi32(regions_reverse_text + (i << 4), vtext_reverse); 
		vseq_ind = _mm512_add_epi32(vseq_ind, vincrement);
		vseq_reverse_ind = _mm512_sub_epi32(vseq_reverse_ind, vincrement);
	}


	__mmask16 vmask = _mm512_int2mask((1 << (read_length & 15)) - 1);
	__m512i zero = _mm512_set1_epi32(0);
	for(uint32_t i = 0; i < 4; i ++){
		//		vseq             = _mm512_mask_i32extgather_epi32(zero, vmask, vseq_ind, sequence, _MM_UPCONV_EPI32_UINT8, 1, 0);
		//		vseq_reverse     = _mm512_mask_i32extgather_epi32(zero, vmask, vseq_reverse_ind, sequence, _MM_UPCONV_EPI32_UINT8, 1, 0);
		/*******************************************/
		__m512i vseq_ind32 = _mm512_srli_epi32(vseq_ind, 2);
		vseq             = _mm512_mask_i32gather_epi32(zero, vmask, vseq_ind32, sequence, 4);
		__m512i vseq_ind8 = _mm512_and_epi32(vseq_ind, vthree);
		vseq_ind8 = _mm512_slli_epi32(vseq_ind8, 3);
		vseq 		= _mm512_srlv_epi32(vseq, vseq_ind8);
		vseq		= _mm512_and_epi32(vseq, chbit);
		/******************************************/

		/********************************************************/
		vseq_ind32 = _mm512_srli_epi32(vseq_reverse_ind, 2);
		vseq_reverse  = _mm512_mask_i32gather_epi32(zero, vmask, vseq_ind32, sequence, 4);
		vseq_ind8 = _mm512_and_epi32(vseq_reverse_ind, vthree);
		vseq_ind8 = _mm512_slli_epi32(vseq_ind8, 3);
		vseq_reverse 	= _mm512_srlv_epi32(vseq_reverse, vseq_ind8);
		vseq_reverse	= _mm512_and_epi32(vseq_reverse, chbit);
		/********************************************************/


		vseq_reverse     = _mm512_permutevar_epi32(vseq_reverse, vreverse_char);
		vseq             = _mm512_slli_epi32(vseq, 24);
		vtext            = _mm512_srli_epi32(vtext, 8);
		vtext            = _mm512_or_epi32(vseq, vtext);
		vseq_reverse     = _mm512_slli_epi32(vseq_reverse, 24);
		vtext_reverse    = _mm512_srli_epi32(vtext_reverse, 8);
		vtext_reverse    = _mm512_or_epi32(vseq_reverse, vtext_reverse);
		vseq_ind         = _mm512_add_epi32(vseq_ind, vone);
		vseq_reverse_ind = _mm512_sub_epi32(vseq_reverse_ind, vone);
		vmask            = _mm512_cmp_epi32_mask(vseq_reverse_ind, zero, _MM_CMPINT_GE);
	}
	_mm512_store_epi32(regions_text + (iteration << 4), vtext);
	_mm512_store_epi32(regions_reverse_text + (iteration << 4), vtext_reverse); 

	//	for(int i = 0; i < read_length * 4){
	//		printf("%c ", regions_text[i]);
	//	}
	//	printf("*************************");
	//
	//	for(int i = 0; i < read_length * 4){
	//		printf("%c ", regions_reverse_text[i]);
	//	}
	//	printf("*************************");




	uint32_t num_kmers        = read_length - kmer_length + 1;
	vseq_ind                  = _mm512_set_epi32(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
	vseq_reverse_ind          = _mm512_set1_epi32(READ_MAX_LEN + num_kmers - 1);
	vseq_reverse_ind          = _mm512_sub_epi32(vseq_reverse_ind, vseq_ind);
	uint32_t trim_kmer_length = kmer_length & 3 ? (kmer_length >> 2) + 1 : (kmer_length >> 2);
	uint32_t *regions_ind     = regions_text + 2 * READ_MAX_LEN;
	__m512i vtext_mask        = _mm512_set1_epi32(0xff);
	uint32_t trim_cmp_length  = (kmer_length + 3)  >> 2;
	trim_cmp_length           = (trim_kmer_length + 1) >> 1;
	__m512i vnum_kmers        = _mm512_set1_epi32(num_kmers - 1);

	for(uint32_t i = 0; i < ((num_kmers + 15) >> 4); i ++){
		//__m512i vseq_kmer         = _mm512_i32gather_epi32(regions_text, );
		//__m512i vseq_reverse_kmer = _mm512_i32gather_epi32(regions_text);
		vmask = _mm512_int2mask(0xffff);
		__mmask16 vseq_mask;
		__mmask16 vseq_reverse_mask, vlt_mask, vgt_mask;
		vlt_mask                  = _mm512_int2mask(0);
		vgt_mask                  = _mm512_int2mask(0);
		__m512i vseq_ind1         = vseq_ind;
		__m512i vseq_reverse_ind1 = vseq_reverse_ind;
		__mmask16 vcontinue       = _mm512_cmp_epi32_mask(vseq_ind, vnum_kmers, _MM_CMPINT_LE);
		vmask                     = _mm512_kand(vmask, vcontinue);

		uint32_t k = 0;
		for(uint32_t m = 0; m < trim_cmp_length; m ++){
			vtext = _mm512_mask_i32gather_epi32(vtext, vmask, vseq_ind1, regions_text, 4);	
			vtext_reverse = _mm512_mask_i32gather_epi32(vtext_reverse, vmask, vseq_reverse_ind1, regions_text, 4);	
			for(uint32_t j = 0; j < 4; j ++){

				vseq              = _mm512_and_epi32(vtext, vtext_mask);
				vseq_reverse      = _mm512_and_epi32(vtext_reverse, vtext_mask);
				vseq_mask         = _mm512_mask_cmp_epi32_mask(vmask, vseq, vseq_reverse, _MM_CMPINT_LT);
				vseq_reverse_mask = _mm512_mask_cmp_epi32_mask(vmask, vseq, vseq_reverse, _MM_CMPINT_GT);
				vlt_mask          = _mm512_kxor(vlt_mask, vseq_mask);
				vgt_mask          = _mm512_kxor(vgt_mask, vseq_reverse_mask);
				vmask             = _mm512_kor(vgt_mask, vlt_mask);
				vmask             = _mm512_kandn(vmask, vcontinue);
				k                 = _mm512_mask2int(vmask);
				if(k == 0)
					break;
				vtext = _mm512_srli_epi32(vtext, 8);
				vtext_reverse = _mm512_srli_epi32(vtext_reverse, 8);
			}

			if(k == 0)
				break;
			vseq_ind1 = _mm512_add_epi32(vseq_ind1, vfour);
			vseq_reverse_ind1 = _mm512_add_epi32(vseq_reverse_ind1, vfour);
		}
		vlt_mask = _mm512_knot(vgt_mask);
		_mm512_mask_store_epi32(regions_ind + (i << 4), vlt_mask, vseq_ind);
		_mm512_mask_store_epi32(regions_ind + (i << 4), vgt_mask, vseq_reverse_ind);
		vseq_ind = _mm512_add_epi32(vseq_ind, vincrement);
		vseq_reverse_ind = _mm512_sub_epi32(vseq_reverse_ind, vincrement);
	}



	__m512i vhash_seed_index = _mm512_set1_epi32(0); 
	__mmask16 vhash_repalce_mask = _mm512_int2mask(0xff);
	__m512i vhash_seed = _mm512_set1_epi32(hash_seed[0]);
	__m512i vhash_num = _mm512_set1_epi32(num_hash_func - 1);
	vnum_kmers = _mm512_set1_epi32(num_kmers);
	__m512i vhash_value;
	__m512i vBIT = _mm512_set_epi32(0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02,
			0x01, 0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01);
	__m512i vbit_vector_width = _mm512_set1_epi32(bit_vector_width);
	__m512i voriginal_index;
	__m512i vbit_vector_index;
	__m512i vbit;
	__m512i vrem_mask = _mm512_set1_epi32(7);
	__m512i vbit_mask;
	__mmask16 vcontinue, vmask_hash_region, vsolid_kmer_mask, vhash_seed_reset_mask, 
		  vhash_seed_increment_mask, vregions_index_increment_mask;
	vcontinue = _mm512_int2mask(0xffff);
	__m512i vregions_index = _mm512_set_epi32(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);


	__m512i vkmer_index = _mm512_load_epi32(regions_text + READ_MAX_LEN * 2);//_mm512_set_epi32(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
	regions_ind[1] = 0;
	vcontinue = _mm512_cmp_epi32_mask(vregions_index, vnum_kmers, _MM_CMPINT_LT);
	/***************************************************/
	//__m512i three = _mm512_set1_epi32(3);
	//__m512i chbit = _mm512_set1_epi32(0xff);
	//__m512i vone = _mm512_set1_epi32(1);
	__m512i vrem32 = _mm512_set1_epi32(0x1f);
	/**************************************************/


	//vcontinue = _mm512_int2mask(0x7f);
	//__m512i mmm = _mm512_set1_epi32(64);
	//vregions_index = _mm512_add_epi32(vregions_index, mmm);
	uint32_t hh = 0;
	while(1){

		//hash value
		//hash function
		MurmurHash3_x86_32(regions_text, vkmer_index, vcontinue, kmer_length, vhash_seed, vhash_value, 4);  



		//voriginal_index = _mm512_rem_epi32(vhash_value, vbit_vector_width);
		//vbit_vector_index = _mm512_srli_epi32(voriginal_index, 3); 
		////vbit_vector_index = _mm512_set1_epi32(1);
		////		vbit = _mm512_mask_i32extgather_epi32(vbit, vcontinue, vbit_vector_index, bit_vector, _MM_UPCONV_EPI32_UINT8, 1, 0);
		////		/***************************************************************/
		//__m512i vseq_ind32 = _mm512_srli_epi32(vbit_vector_index, 2);
		//vbit             = _mm512_mask_i32gather_epi32(zero, vmask, vseq_ind32, bit_vector, 4);
		//__m512i vseq_ind8 = _mm512_and_epi32(vbit_vector_index, three);
		//vseq_ind8 = _mm512_slli_epi32(vseq_ind8, 3);
		//vbit 		= _mm512_srlv_epi32(vbit, vseq_ind8);
		//vbit		= _mm512_and_epi32(vbit, chbit);
		///****************************************************************/

		//vbit_mask = _mm512_permutevar_epi32(voriginal_index, vBIT);
		//vbit = _mm512_and_epi32(vbit, vbit_mask);
		//vmask = _mm512_cmp_epi32_mask(vbit, vbit_mask, _MM_CMPINT_EQ);

/*************************************************************************/
		voriginal_index = _mm512_rem_epu32(vhash_value, vbit_vector_width);
		vbit_vector_index = _mm512_srli_epi32(voriginal_index, 5);
		vbit             = _mm512_mask_i32gather_epi32(vbit, vcontinue, vbit_vector_index, bit_vector, 4);
		vbit_mask = _mm512_and_epi32(voriginal_index, vrem32);
		vbit_mask =  _mm512_sllv_epi32(vone, vbit_mask);
		vbit     = _mm512_and_epi32(vbit, vbit_mask);
		vmask = _mm512_cmp_epi32_mask(vbit, vbit_mask, _MM_CMPINT_EQ);
/*************************************************************************/



		vmask_hash_region = _mm512_cmp_epi32_mask(vhash_seed_index, vhash_num, _MM_CMPINT_EQ);
		vsolid_kmer_mask = _mm512_kand(vmask, vmask_hash_region);
		vsolid_kmer_mask = _mm512_kand(vcontinue, vsolid_kmer_mask);

		_mm512_mask_i32scatter_epi32(regions_ind, vsolid_kmer_mask, vregions_index, vone, 4);
		//
		vhash_seed_reset_mask = _mm512_knot(vmask);
		vhash_seed_reset_mask =	_mm512_kor(vhash_seed_reset_mask, vsolid_kmer_mask);
		vhash_seed_reset_mask = _mm512_kand(vhash_seed_reset_mask, vcontinue);

		vhash_seed_index = _mm512_mask_mov_epi32(vhash_seed_index, vhash_seed_reset_mask, zero);
		vhash_seed_increment_mask = _mm512_kandn(vhash_seed_reset_mask, vcontinue);

		vhash_seed_index = _mm512_mask_add_epi32(vhash_seed_index, vhash_seed_increment_mask, vhash_seed_index, vone);
		vhash_seed = _mm512_mask_i32gather_epi32(vhash_seed, vcontinue, vhash_seed_index, hash_seed, 4);
		//break;	
		//vregions_index_increment_mask = _mm512_kor(vhash_seed_reset_mask, vsolid_kmer_mask);


		vregions_index = _mm512_mask_add_epi32(vregions_index, vhash_seed_reset_mask, vregions_index, vincrement);
		vcontinue = _mm512_cmp_epi32_mask(vregions_index, vnum_kmers, _MM_CMPINT_LT);

		vhash_seed_reset_mask = _mm512_kand(vhash_seed_reset_mask, vcontinue);
		vkmer_index = _mm512_mask_i32gather_epi32(vkmer_index, vhash_seed_reset_mask, vregions_index, regions_ind, 4);
		uint32_t k = _mm512_mask2int(vcontinue);
		if(k == 0)
			break;

	}
}

void pos_path_query_text(char *kmer_buffer, const uint8_t *bit_vector, const uint32_t *hash_seed, uint32_t kmer_length, uint32_t num_hash_func, uint32_t bit_vector_width, uint32_t pos, bool *pass_result){

	//#ifdef __MIC__
	uint32_t kmer_int_size = KMER_MAX_LEN / 4;	

	char save_char = kmer_buffer[kmer_length - 1 - pos];
	switch(save_char){
		case 'A':
			kmer_buffer[pos] = 'T';		
			break;
		case 'C':
			kmer_buffer[pos] = 'G';		
			break;
		case 'G':
			kmer_buffer[pos] = 'C';		
			break;
		case 'T' :
			kmer_buffer[pos] = 'A';		
	}

	int comparison_result = 1;
	int compare_ind = 0;
	for (unsigned int it = 0, it_rc = (kmer_length - 1); it < kmer_length;
			it++, it_rc--) {
		compare_ind = it;
		if (kmer_buffer[it] == 'A') {
			if (kmer_buffer[it_rc] == 'T') {
			} else {
				comparison_result = 1;
				break;
			}
		} else if (kmer_buffer[it] == 'C') {
			if (kmer_buffer[it_rc] == 'G') {
			} else if (kmer_buffer[it_rc] == 'T') {
				comparison_result = -1;
				break;
			} else {
				comparison_result = 1;
				break;
			}
		} else if (kmer_buffer[it] == 'G') {
			if (kmer_buffer[it_rc] == 'A') {
				comparison_result = 1;
				break;
			} else if (kmer_buffer[it_rc] == 'C') {
			} else {
				comparison_result = -1;
				break;
			}
		} else if (kmer_buffer[it] == 'T') {
			if (kmer_buffer[it_rc] == 'A') {
			} else {
				comparison_result = -1;
				break;
			}
		} else {
			//std::cout << std::endl << "ERROR: Illegal character " << kmer[it] << " (query_text)" << std::endl << std::endl;
			printf("ERROR: 1 Illegal charactor %c query_text\n\n", kmer_buffer[it]);
		}
	}
	kmer_buffer[pos] = 0;     

	uint32_t c_ind = kmer_length - 1 - pos;
	uint32_t kk = c_ind < pos ? c_ind : pos;

	//	__m512i vreverse_char   = _mm512_set_epi32(0, 0, 0, 0, 0, 0, 0, 0, 0x43, 0, 0, 0x41, 0x47, 0, 0x54, 0);

	//	__m512i vkmer_length = _mm512_set1_epi32(kmer_length);
	//	__m512i vkmer_index = _mm512_set_epi32(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
	//	__m512i vkmer_reverse_index = _mm512_set1_epi32(KMER_MAX_LEN + kmer_length - 1);
	//	vkmer_reverse_index = _mm512_sub_epi32(vkmer_reverse_index, vkmer_index);
	//	__mmask16 vcontinue = _mm512_cmp_epi32_mask(vkmer_index, vkmer_length, _MM_CMPINT_LT);
	//	__m512i vkmer, vkmer_reverse;
	//	__m512i vincrement = _mm512_set1_epi32(16);
	//    if(compare_ind >= kk || comparison_result == -1){

	//	    for(uint32_t i = 0; i < ((kmer_length + 15) >> 4); i ++){
	//		    vkmer = _mm512_mask_i32extgather_epi32(vkmer, vcontinue, vkmer_index, kmer_buffer, _MM_UPCONV_EPI32_UINT8, 1, 0);
	//		    vkmer_reverse = _mm512_permutevar_epi32(vkmer, vreverse_char);
	//		    _mm512_mask_i32extscatter_epi32(kmer_buffer, vcontinue, vkmer_reverse_index, vkmer_reverse, _MM_DOWNCONV_EPI32_UINT8, 1, 0);
	//		    vkmer_index = _mm512_add_epi32(vkmer_index, vincrement);
	//		    vkmer_reverse_index = _mm512_sub_epi32(vkmer_reverse_index, vincrement);
	//		    vcontinue = _mm512_cmp_epi32_mask(vkmer_index, vkmer_length, _MM_CMPINT_LT);
	//	    }
	//	    }

	/********************************************************************************/
	__m512i vreverse_char   = _mm512_set_epi32(0, 0, 0, 0, 0, 0, 0, 0, 0x43, 0, 0, 0x41, 0x47, 0, 0x54, 0);
	__m512i vkmer_length = _mm512_set1_epi32(kmer_length);
	__m512i vkmer_index8 = _mm512_set_epi32(0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60);
	__m512i vkmer_index = _mm512_set_epi32(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);


	__m512i vkmer_reverse_index = _mm512_set1_epi32(KMER_MAX_LEN + kmer_length - 1);

	vkmer_reverse_index = _mm512_sub_epi32(vkmer_reverse_index, vkmer_index8);
	__mmask16 vcontinue = _mm512_cmp_epi32_mask(vkmer_index8, vkmer_length, _MM_CMPINT_LT);
	__m512i vkmer, vkmer_reverse;
	__m512i vincrement = _mm512_set1_epi32(16);
//	for(uint32_t i = 0; i < kmer_length; i ++){
//		printf("%c", kmer_buffer[i]);
//	}
//	printf("\n");


	if(compare_ind >= kk || comparison_result == -1){

		__m512i vone = _mm512_set1_epi32(1);
		__m512i vbase_addr = _mm512_set1_epi32(KMER_MAX_LEN>>2);
		__m512i vsixty = _mm512_set1_epi32(60);
		__m512i vsixtyfour = _mm512_set1_epi32(64);
		vkmer_reverse_index = _mm512_srli_epi32(vkmer_reverse_index, 2);

		for(int m = 0; m < kmer_length; m += 64){
			vkmer = _mm512_mask_i32gather_epi32(vkmer, vcontinue, vkmer_index, kmer_buffer, 4);

			vkmer_reverse = _mm512_set1_epi32(0);

			for(int j = 0; j < 4; j ++){
				__m512i vkmer_re8 = _mm512_mask_permutevar_epi32(vkmer_reverse, vcontinue, vkmer, vreverse_char); 

//				int *tmp1 = (int*)_mm_malloc(16 * sizeof(int), 64);
//				_mm512_store_epi32(tmp1, vkmer_reverse_index);
//				for(int i = 0; i < 16; i ++)
//					printf("%d ", tmp1[i]);
//				printf("\n");
//
//				_mm512_store_epi32(tmp1, vkmer);
//				for(int i = 0; i < 16; i ++)
//					printf("%c ", tmp1[i]);
//				printf("\n");
//
//				_mm512_store_epi32(tmp1, vkmer_re8);
//				for(int i = 0; i < 16; i ++)
//					printf("%c ", tmp1[i]);
//				printf("\n");
//
				vkmer_reverse = _mm512_mask_slli_epi32(vkmer_reverse, vcontinue, vkmer_reverse, 8); 
				vkmer_reverse = _mm512_mask_or_epi32(vkmer_reverse, vcontinue, vkmer_reverse, vkmer_re8);
//
//				_mm512_store_epi32(tmp1, vkmer_reverse);
//				for(int i = 0; i < 16; i ++)
//					printf("%x ", tmp1[i]);
//				printf("\n");
//				//vkmer_index8 = _mm512_add_epi32(vkmer_index8, vone);
				//vcontinue = _mm512_cmp_epi32_mask(vkmer_index8, vkmer_length, _MM_CMPINT_LT);

				vkmer = _mm512_mask_srli_epi32(vkmer, vcontinue, vkmer, 8); 
			}
			if((kmer_length & 3) != 0){

				__m512i vpermt = _mm512_set_epi32(16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1);
				//__m512i vpermt = _mm512_sub_epi32(vsi, vone);

				vkmer = _mm512_permutevar_epi32(vpermt, vkmer_reverse); 
//				char *tmp1 = (char*)_mm_malloc(16 * sizeof(int), 64);
//				printf("**********************\n");
//
//				_mm512_store_epi32(tmp1, vkmer_reverse);
//				for(int i = 0; i < 64; i ++)
//					printf("%x ", tmp1[i]);
//				printf("\n");

//				_mm512_store_epi32(tmp1, vkmer);
//				for(int i = 0; i < 64; i ++)
//					printf("%x ", tmp1[i]);
//				printf("\n");
//				//__m512i vshift = _mm512_set1_epi32((kmer_length & 3) << 3);  

				vkmer = _mm512_slli_epi32(vkmer, (kmer_length & 3) << 3);
//				_mm512_store_epi32(tmp1, vkmer);
//				for(int i = 0; i < 64; i ++)
//					printf("%x ", tmp1[i]);
//				printf("\n");


				vkmer_reverse = _mm512_srli_epi32(vkmer_reverse, (4 - (kmer_length & 3)) << 3); 
//				_mm512_store_epi32(tmp1, vkmer_reverse);
//				for(int i = 0; i < 64; i ++)
//					printf("%x ", tmp1[i]);
//				printf("\n");

				vkmer_reverse = _mm512_or_epi32(vkmer_reverse, vkmer);

//				_mm512_store_epi32(tmp1, vkmer_reverse);
//				for(int i = 0; i < 64; i ++)
//					printf("%x ", tmp1[i]);
//				printf("**********************\n");

			}

			vcontinue = _mm512_cmp_epi32_mask(vkmer_reverse_index, vbase_addr, _MM_CMPINT_GE);
			_mm512_mask_i32scatter_epi32(kmer_buffer, vcontinue, vkmer_reverse_index, vkmer_reverse, 4);
			vkmer_reverse_index = _mm512_sub_epi32(vkmer_reverse_index, vincrement);

			vkmer_index = _mm512_add_epi32(vkmer_index, vincrement);
			vkmer_index8 = _mm512_add_epi32(vkmer_index8, vsixty);
			vcontinue = _mm512_cmp_epi32_mask(vkmer_index8, vkmer_length, _MM_CMPINT_LT);
		}		

	}

	//	    /********************************************************************************/

//	for(uint32_t i = 0; i < kmer_length; i ++){
//		printf("%x ", kmer_buffer[KMER_MAX_LEN + i]);
//	}
//	printf("\n");
//
//	for(uint32_t i = 0; i < kmer_length; i ++){
//		printf("%c ", kmer_buffer[KMER_MAX_LEN + i]);
//	}
//	printf("\n");



	uint32_t kmer_pos[4];
	kmer_pos[0] = kmer_pos[1] = kmer_pos[2] = kmer_pos[3] = 1;
	char cmp_base;
	cmp_base = kmer_buffer[c_ind];
	uint32_t rest_kmer = 0;
	if(compare_ind > kk){ 
		switch(cmp_base){
			case 'A':
				kmer_pos[0] = 0;
				kmer_pos[1] = 0;
				kmer_pos[2] = 0;
				rest_kmer   = 3;
				break;
			case 'C':
				kmer_pos[0] = 0;
				kmer_pos[1] = 0;
				rest_kmer   = 2;
				kmer_pos[3] = kmer_int_size;
				break;
			case 'G':
				kmer_pos[0] = 0;
				rest_kmer   = 1;
				kmer_pos[2] = kmer_int_size;
				kmer_pos[3] = kmer_int_size;
				break;
			case 'T' :
				rest_kmer   = 0;
				kmer_pos[1] = kmer_int_size;
				kmer_pos[2] = kmer_int_size;
				kmer_pos[3] = kmer_int_size;
		}
		if(comparison_result == 1){
			kmer_pos[rest_kmer] = 0;
		}else{
			kmer_pos[rest_kmer] = kmer_int_size;
		}
	}else {
		if(comparison_result == 1){
			kmer_pos[0] = kmer_pos[1] = kmer_pos[2] = kmer_pos[3] = 0;
		}else{
			kmer_pos[0] = kmer_pos[1] = kmer_pos[2] = kmer_pos[3] = kmer_int_size;
		}
	}



	//__m512i four = _mm512_set1_epi32(4);
	__m512i vhash_seed_index = _mm512_set_epi32(3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0);
	__m512i vhash_seed;	
	__m512i vtail = _mm512_set_epi32('T', 'T', 'T', 'T', 'G', 'G', 'G', 'G', 'C', 'C', 'C', 'C', 'A', 'A', 'A', 'A');	
	uint32_t increment = 4;
	uint32_t result_mask = 0xf;
	uint32_t num_trim_base = 4;
	vkmer_index = _mm512_set_epi32(kmer_pos[3], kmer_pos[3], kmer_pos[3], kmer_pos[3], kmer_pos[2], kmer_pos[2], kmer_pos[2], kmer_pos[2], 
			kmer_pos[1], kmer_pos[1], kmer_pos[1], kmer_pos[1], kmer_pos[0], kmer_pos[0], kmer_pos[0], kmer_pos[0]);

	vincrement = _mm512_set1_epi32(4);
	//printf("%d %d %d %d \n", kmer_pos[0], kmer_pos[1], kmer_pos[2], kmer_pos[3]);


	__m512i vtext, vtext_reverse;
	__m512i vhash_seed_num = _mm512_set1_epi32(num_hash_func - 1);
	vhash_seed = _mm512_set1_epi32(kmer_int_size);
	__mmask16 vtail_mask = _mm512_cmp_epi32_mask(vhash_seed, vkmer_index, _MM_CMPINT_EQ);

	vtail = _mm512_mask_permutevar_epi32(vtail, vtail_mask, vtail, vreverse_char);
	__m512i vpos = _mm512_set1_epi32(pos);
	__m512i vindex = _mm512_set1_epi32(kmer_length - 1);
	vpos = _mm512_mask_sub_epi32(vpos, vtail_mask, vindex, vpos);  
	vindex = _mm512_set1_epi32(3);
	vpos = _mm512_and_epi32(vpos, vindex); 
	vpos = _mm512_slli_epi32(vpos,  3);
	vtail = _mm512_sllv_epi32(vtail, vpos);

	vcontinue  = _mm512_cmp_epi32_mask(vhash_seed_index, vhash_seed_num, _MM_CMPINT_LE);
	vhash_seed = _mm512_mask_i32gather_epi32(vhash_seed, vcontinue, vhash_seed_index, hash_seed, 4);
	__m512i vhash_value;
	__m512i voriginal_index;
	__m512i vBIT = _mm512_set_epi32(0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01, 0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01);
      	__m512i vbit_vector_width = _mm512_set1_epi32(bit_vector_width);
	__m512i vbit_vector_index;
	__m512i vbit;
	//__m512i vrem_mask = _mm512_set1_epi32(7);
	__m512i vbit_mask;

	uint32_t result = num_hash_func >= increment ? result_mask : (1 << num_hash_func) - 1;
	uint32_t hash_remainder = num_hash_func;
	pass_result[0] = 1;
	pass_result[1] = 1;
	pass_result[2] = 1;
	pass_result[3] = 1;
	//__m512i four = _mm512_set1_epi32(4);
	__mmask16 vmask;
	/***************************************************/
	__m512i three = _mm512_set1_epi32(3);
	__m512i chbit = _mm512_set1_epi32(0xff);
	__m512i vone = _mm512_set1_epi32(1);
	__m512i vrem32 = _mm512_set1_epi32(0x1f);
	/**************************************************/
	/*
	   int *tmp1 = (int*)_mm_malloc(16 * sizeof(int), 64);
	   _mm512_store_epi32(tmp1, vkmer_index);
	   for(int i = 0; i < 16; i ++)
	   printf("%d ", tmp1[i]);
	   printf("\n");
	   */
	/*
	   uint32_t *tt = (uint32_t*)_mm_malloc(16 * sizeof(uint32_t), 64);
	   printf("tail mask %x \n", vtail_mask);
	//vtail = _mm512_slli_epi32(vtail, 24);
	_mm512_store_epi32(tt, vkmer_index);
	for(uint32_t i = 0; i < 16; i ++){
	printf("%x ", tt[i]);
	}
	printf("\n");
	_mm512_store_epi32(tt, vtail);
	for(uint32_t i = 0; i < 16; i ++){
	printf("%x ", tt[i]);
	}
	printf("\n");
	for(uint32_t i = 0; i < kmer_length; i ++){
	printf("%c", kmer_buffer[36 + i]);
	}
	printf("\n");
	uint32_t *mm = (uint32_t*)kmer_buffer;
	for(uint32_t i = 0; i < 72/4; i ++){
	printf("%x  ", mm[i]);
	}	
	printf("hhhhhh\n");
	*/

	while(1){
		//compute hash value	

		MurmurHash3_x86_32(kmer_buffer, vkmer_index, vcontinue, kmer_length, vhash_seed, vhash_value, vtail, vtail_mask, pos);  


//		uint32_t *tmp1 = (uint32_t*)_mm_malloc(16 * sizeof(int), 64);
//		uint32_t *tmp2 = (uint32_t*)_mm_malloc(16 * sizeof(int), 64);

//		//printf("%x \n", k);
		//vhash_value = _mm512_set1_epi32(0xffffffff);

		voriginal_index = _mm512_rem_epu32(vhash_value, vbit_vector_width);

//		uint32_t *tmp1 = (uint32_t*)_mm_malloc(16 * sizeof(int), 64);
//		_mm512_store_epi32(tmp1, vbit_vector_width);
//		for(int i = 0; i < 16; i ++)
//			printf("%x ", tmp1[i]);
//		printf("\n");
//
//
//		_mm512_store_epi32(tmp1, vhash_value);
//		for(int i = 0; i < 16; i ++)
//			printf("%x ", tmp1[i]);
//		printf("\n");
//
//		printf("original index \n");
//		_mm512_store_epi32(tmp2, voriginal_index);
//		for(int i = 0; i < 16; i ++)
//			printf("%x ", tmp2[i]);
//		printf("\n");

		//vbit_vector_index = _mm512_srli_epi32(voriginal_index, 3); 
		/****************************************************************************/
		//		vbit = _mm512_mask_i32extgather_epi32(vbit, vcontinue, vbit_vector_index, bit_vector, _MM_UPCONV_EPI32_UINT8, 1, 0);
		//
		
		vbit_vector_index = _mm512_srli_epi32(voriginal_index, 5);
		vbit             = _mm512_mask_i32gather_epi32(vbit, vcontinue, vbit_vector_index, bit_vector, 4);

//		printf("mask %x \n", vcontinue);
//
//		_mm512_store_epi32(tmp1, vbit);
//		for(int i = 0; i < 16; i ++)
//			printf("%x ", tmp1[i]);
//		printf("\n");
//
//		_mm512_store_epi32(tmp1, vbit_vector_index);
//		for(int i = 0; i < 16; i ++){
//			printf("%x ", tmp1[i]);
//			if(tmp1[i] * 4 >  tmp2[i]/8){
//				printf("wrong !!! %d %u %u\n", i, tmp1[i]*4, tmp2[i]);
//			}
//		}
//		printf("\n");

//		printf("%x %x %x %x\n", bit_vector[tmp1[0] * 4], bit_vector[tmp1[0] * 4 + 1], 
//			bit_vector[tmp1[0] * 4 + 2 ], bit_vector[tmp1[0] * 4 + 3]);
//

//		_mm512_store_epi32(tmp1, vbit);
//		for(int i = 0; i < 16; i ++)
//			printf("%x ", tmp1[i]);
//		printf("\n");



		/*****************************************************************/	
//		vbit_mask = _mm512_permutevar_epi32(voriginal_index, vBIT);
//		vbit = _mm512_and_epi32(vbit, vbit_mask);
//		printf("%x %x %x %x\n", bit_vector[tmp1[0] * 4], bit_vector[tmp1[0] * 4 + 1], 
//			bit_vector[tmp1[0] * 4 + 2 ], bit_vector[tmp1[0] * 4 + 3]);


		vbit_mask = _mm512_and_epi32(voriginal_index, vrem32);
		vbit_mask =  _mm512_sllv_epi32(vone, vbit_mask);
		vbit     = _mm512_and_epi32(vbit, vbit_mask);
		vmask = _mm512_cmp_epi32_mask(vbit, vbit_mask, _MM_CMPINT_EQ);

//		_mm512_store_epi32(tmp1, vbit);
//		for(int i = 0; i < 16; i ++)
//			printf("%x ", tmp1[i]);
//		printf("\n");



		/*****************************************************************/	
//		vbit_mask = _mm512_permutevar_epi32(voriginal_index, vBIT);
//		vbit = _mm512_and_epi32(vbit, vbit_mask);
//		vmask = _mm512_cmp_epi32_mask(vbit, vbit_mask, _MM_CMPINT_EQ);
		uint32_t k = _mm512_mask2int(vmask);
		uint32_t candidate_solid_kmer_mask = 0;



		for(uint32_t i = 0; i < num_trim_base; i ++){
			if(((k & result) == result) && pass_result[i]){
				pass_result[i] = true;
				candidate_solid_kmer_mask += (result_mask << (increment * i));
			}else{
				pass_result[i] = false;
			}
			k >>= increment;			
		}	
		hash_remainder -= increment;
		result = hash_remainder >= increment ? result_mask : (1 << hash_remainder) - 1;
		vhash_seed_index = _mm512_add_epi32(vhash_seed_index, vincrement);	
		vcontinue = _mm512_cmp_epi32_mask(vhash_seed_index, vhash_seed_num, _MM_CMPINT_LE);
		__mmask16 vcandidate_solid_kmer_mask = _mm512_int2mask(candidate_solid_kmer_mask);
		vcontinue = _mm512_kand(vcontinue, vcandidate_solid_kmer_mask);
		k = _mm512_mask2int(vcontinue);
		if(k == 0)
			break;
		vhash_seed = _mm512_mask_i32gather_epi32(vhash_seed, vcontinue, vhash_seed_index, hash_seed, 4);

	}

	//_mm512_store_epi32(tt, vhash_value);
	/*	
		for(uint32_t i = 0; i < 4; i ++){
		printf("%x ", pass_result[i]);
		}
		printf("\n");
		*/

	//break;
	//_mm_free(tt);

	//#endif

}
__ONMIC__ void path_query_text(char *kmer_buffer, const uint8_t *bit_vector, const uint32_t *hash_seed, uint32_t kmer_length, uint32_t num_hash_func, uint32_t bit_vector_width, uint32_t tail_flag, char trim_base, bool *pass_result){

	//#ifdef __MIC__
	uint32_t kmer_int_size = KMER_MAX_LEN / 4;	
	char save_char;
	if(tail_flag == 0){
		save_char = kmer_buffer[0];
		kmer_buffer[0] = 0;
	}else if(tail_flag == 1){
		save_char = kmer_buffer[kmer_length - 1];
		kmer_buffer[kmer_length - 1] = 0;
	}

	/*	
		for(uint32_t i = 0; i < kmer_length - 1; i ++){
		printf("%c", kmer_buffer[i]);
		}
		printf("\n");
		*/	

	//	__m512i vreverse_char    = _mm512_set_epi32(0, 0, 0, 0, 0, 0, 0, 0, 0x43, 0, 0, 0x41, 0x47, 0, 0x54, 0);
	//	__m512i vkmer_length = _mm512_set1_epi32(kmer_length);
	//	__m512i vkmer_index = _mm512_set_epi32(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
	//	__m512i vkmer_reverse_index = _mm512_set1_epi32(KMER_MAX_LEN + kmer_length - 1);
	//	vkmer_reverse_index = _mm512_sub_epi32(vkmer_reverse_index, vkmer_index);
	//	__mmask16 vcontinue = _mm512_cmp_epi32_mask(vkmer_index, vkmer_length, _MM_CMPINT_LT);
	//	__m512i vkmer, vkmer_reverse;
	//	__m512i vincrement = _mm512_set1_epi32(16);
	//	for(uint32_t i = 0; i < ((kmer_length + 15) >> 4); i ++){
	//		vkmer = _mm512_mask_i32extgather_epi32(vkmer, vcontinue, vkmer_index, kmer_buffer, _MM_UPCONV_EPI32_UINT8, 1, 0);
	//		vkmer_reverse = _mm512_permutevar_epi32(vkmer, vreverse_char);
	//		_mm512_mask_i32extscatter_epi32(kmer_buffer, vcontinue, vkmer_reverse_index, vkmer_reverse, _MM_DOWNCONV_EPI32_UINT8, 1, 0);
	//		vkmer_index = _mm512_add_epi32(vkmer_index, vincrement);
	//		vkmer_reverse_index = _mm512_sub_epi32(vkmer_reverse_index, vincrement);
	//		vcontinue = _mm512_cmp_epi32_mask(vkmer_index, vkmer_length, _MM_CMPINT_LT);
	//	}
	/********************************************************************************/
	//__m512i vreverse_char    = _mm512_set_epi32(0, 0, 0, 0, 0, 0, 0, 0, 0x43, 0, 0, 0x41, 0x47, 0, 0x54, 0);
	//__m512i vkmer_length = _mm512_set1_epi32(kmer_length);
	//__m512i vkmer_index8 = _mm512_set_epi32(0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60);
	//__m512i vkmer_index = _mm512_set_epi32(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);


	//__m512i vkmer_reverse_index = _mm512_set1_epi32(KMER_MAX_LEN + kmer_length - 1);

	//vkmer_reverse_index = _mm512_sub_epi32(vkmer_reverse_index, vkmer_index8);
	//__mmask16 vcontinue = _mm512_cmp_epi32_mask(vkmer_index8, vkmer_length, _MM_CMPINT_LT);
	//__m512i vkmer, vkmer_reverse;
	//__m512i vincrement = _mm512_set1_epi32(16);
	//__m512i vone = _mm512_set1_epi32(1);
	//__m512i vzero = _mm512_set1_epi32(0);
	//__m512i vsixty = _mm512_set1_epi32(60);
	//__m512i vsixtyfour = _mm512_set1_epi32(64);
	//vkmer_reverse_index = _mm512_srli_epi32(vkmer_reverse_index, 2);
	//for(int i = 0; i < kmer_length; i += 64){
	//	vkmer = _mm512_mask_i32gather_epi32(vkmer, vcontinue, vkmer_index, kmer_buffer, 4);
	//	for(int j = 0; j < 4; j ++){
	//		__m512i vkmer_re8 = _mm512_mask_permutevar_epi32(vkmer_reverse, vcontinue, vkmer, vreverse_char); 
	//		vkmer_reverse = _mm512_mask_slli_epi32(vkmer_reverse, vcontinue, vkmer_reverse, 8); 
	//		vkmer_reverse = _mm512_mask_or_epi32(vkmer_reverse, vcontinue, vkmer_reverse, vkmer_re8);
	//		vkmer_index8 = _mm512_add_epi32(vkmer_index8, vone);
	//		vcontinue = _mm512_cmp_epi32_mask(vkmer_index8, vkmer_length, _MM_CMPINT_LT);
	//	}

	//	vcontinue = _mm512_cmp_epi32_mask(vkmer_reverse_index, vzero, _MM_CMPINT_GE);
	//	_mm512_mask_i32scatter_epi32(kmer_buffer, vcontinue, vkmer_reverse_index, vkmer_reverse, 4);
	//	vkmer_reverse_index = _mm512_sub_epi32(vkmer_reverse_index, vincrement);

	//	vkmer_index = _mm512_add_epi32(vkmer_index, vincrement);
	//	vkmer_index8 = _mm512_add_epi32(vkmer_index8, vsixty);
	//	vcontinue = _mm512_cmp_epi32_mask(vkmer_index8, vkmer_length, _MM_CMPINT_LT);
	//}		

	/*******************************************************************************/
	__m512i vreverse_char   = _mm512_set_epi32(0, 0, 0, 0, 0, 0, 0, 0, 0x43, 0, 0, 0x41, 0x47, 0, 0x54, 0);
	__m512i vkmer_length = _mm512_set1_epi32(kmer_length);
	__m512i vkmer_index8 = _mm512_set_epi32(0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60);
	__m512i vkmer_index = _mm512_set_epi32(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);


	__m512i vkmer_reverse_index = _mm512_set1_epi32(KMER_MAX_LEN + kmer_length - 1);

	vkmer_reverse_index = _mm512_sub_epi32(vkmer_reverse_index, vkmer_index8);
	__mmask16 vcontinue = _mm512_cmp_epi32_mask(vkmer_index8, vkmer_length, _MM_CMPINT_LT);
	__m512i vkmer, vkmer_reverse;
	__m512i vincrement = _mm512_set1_epi32(16);



	__m512i vone = _mm512_set1_epi32(1);
	__m512i vbase_addr = _mm512_set1_epi32(KMER_MAX_LEN>>2);
	__m512i vsixty = _mm512_set1_epi32(60);
	__m512i vsixtyfour = _mm512_set1_epi32(64);
	vkmer_reverse_index = _mm512_srli_epi32(vkmer_reverse_index, 2);

	for(int m = 0; m < kmer_length; m += 64){
		vkmer = _mm512_mask_i32gather_epi32(vkmer, vcontinue, vkmer_index, kmer_buffer, 4);
		vkmer_reverse = _mm512_set1_epi32(0);
		for(int j = 0; j < 4; j ++){
			__m512i vkmer_re8 = _mm512_mask_permutevar_epi32(vkmer_reverse, vcontinue, vkmer, vreverse_char); 
			vkmer_reverse = _mm512_mask_slli_epi32(vkmer_reverse, vcontinue, vkmer_reverse, 8); 
			vkmer_reverse = _mm512_mask_or_epi32(vkmer_reverse, vcontinue, vkmer_reverse, vkmer_re8);
			vkmer = _mm512_mask_srli_epi32(vkmer, vcontinue, vkmer, 8); 
		}
		if((kmer_length & 3) != 0){
			__m512i vpermt = _mm512_set_epi32(16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1);
			vkmer = _mm512_permutevar_epi32(vpermt, vkmer_reverse); 
			vkmer = _mm512_slli_epi32(vkmer, (kmer_length & 3) << 3);
			vkmer_reverse = _mm512_srli_epi32(vkmer_reverse, (4 - (kmer_length & 3)) << 3); 
			vkmer_reverse = _mm512_or_epi32(vkmer_reverse, vkmer);
		}

		vcontinue = _mm512_cmp_epi32_mask(vkmer_reverse_index, vbase_addr, _MM_CMPINT_GE);
		_mm512_mask_i32scatter_epi32(kmer_buffer, vcontinue, vkmer_reverse_index, vkmer_reverse, 4);
		vkmer_reverse_index = _mm512_sub_epi32(vkmer_reverse_index, vincrement);

		vkmer_index = _mm512_add_epi32(vkmer_index, vincrement);
		vkmer_index8 = _mm512_add_epi32(vkmer_index8, vsixty);
		vcontinue = _mm512_cmp_epi32_mask(vkmer_index8, vkmer_length, _MM_CMPINT_LT);
	}		


    /********************************************************************************/


	/*	
		for(uint32_t i = 1; i < kmer_length; i ++){
		printf("%c", kmer_buffer[KMER_MAX_LEN + i]);
		}
		printf("\n");
		*/
	uint32_t kmer_pos[4];
	kmer_pos[0] = kmer_pos[1] = kmer_pos[2] = kmer_pos[3] = 1;
	char cmp_base;
	if(tail_flag == 0){
		cmp_base = kmer_buffer[kmer_length - 1];
	}else if(tail_flag == 1){
		cmp_base = kmer_buffer[0];
	}
	uint32_t rest_kmer = 0;

	switch(cmp_base){
		case 'A':
			kmer_pos[0] = 0;
			kmer_pos[1] = 0;
			kmer_pos[2] = 0;
			rest_kmer   = 3;
			break;
		case 'C':
			kmer_pos[0] = 0;
			kmer_pos[1] = 0;
			rest_kmer   = 2;
			kmer_pos[3] = kmer_int_size;
			break;
		case 'G':
			kmer_pos[0] = 0;
			rest_kmer   = 1;
			kmer_pos[2] = kmer_int_size;
			kmer_pos[3] = kmer_int_size;
			break;
		case 'T' :
			rest_kmer   = 0;
			kmer_pos[1] = kmer_int_size;
			kmer_pos[2] = kmer_int_size;
			kmer_pos[3] = kmer_int_size;
	}

	uint32_t trim_flag = 4;
	switch(trim_base){
		case 'A':
			trim_flag = 0;
			break;
		case 'C':
			trim_flag = 1;
			break;
		case 'G':
			trim_flag = 2;
			break;
		case 'T' :
			trim_flag = 3;
			break;
		default:
			trim_flag = 4;
	}


	if(trim_flag != rest_kmer){
		uint32_t ii = 1, jj = 1;
		if(tail_flag == 0){
			ii = 1;
			jj = 1;
		}
		for(uint32_t i = 0; i < kmer_length / 2; i ++){
			if(kmer_buffer[ii + i] < kmer_buffer[KMER_MAX_LEN + jj + i]){ 
				kmer_pos[rest_kmer] = 0;
				break;
			}else if(kmer_buffer[ii + i] > kmer_buffer[KMER_MAX_LEN + jj + i]){ 
				kmer_pos[rest_kmer] = kmer_int_size;
				break;
			}
		}
		/*
		   vkmer_index = _mm512_set_epi32(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);
		   __m512i one = _mm512_set1_epi32(1);
		   vkmer_index = _mm512_add_epi32(vkmer_index, one);
		   vkmer_reverse_index = _mm512_set1_epi32(KMER_MAX_LEN);
		   vkmer_reverse_index = _mm512_add_epi32(vkmer_reverse_index, vkmer_index);
		   vcontinue = _mm512_cmp_epi32_mask(vkmer_index, vkmer_length, _MM_CMPINT_LT);
		   kmer_pos[rest_kmer] = 0;
		   for(uint32_t i = 0; i < (kmer_length + 15) >> 4; i ++){	
		   vkmer = _mm512_mask_i32extgather_epi32(vkmer, vcontinue, vkmer_index, kmer_buffer, _MM_UPCONV_EPI32_UINT8, 1, 0);
		   vkmer_reverse = _mm512_mask_i32extgather_epi32(vkmer_reverse, vcontinue, vkmer_reverse_index, kmer_buffer, _MM_UPCONV_EPI32_UINT8, 1, 0);
		   __mmask16 vlt_mask = _mm512_mask_cmp_epi32_mask(vcontinue, vkmer, vkmer_reverse, _MM_CMPINT_LT);	
		   __mmask16 vgt_mask = _mm512_mask_cmp_epi32_mask(vcontinue, vkmer, vkmer_reverse, _MM_CMPINT_GT);	
		   if(vlt_mask > vgt_mask){
		   kmer_pos[rest_kmer] = 0;
		   break;
		   }if(vlt_mask < vgt_mask){
		   kmer_pos[rest_kmer] = 8;
		   break;
		   }
		   }
		   */
	}

	__m512i vhash_seed_index = _mm512_set_epi32(3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0);
	__m512i vhash_seed;	
	__m512i vtail = _mm512_set_epi32('T', 'T', 'T', 'T', 'G', 'G', 'G', 'G', 'C', 'C', 'C', 'C', 'A', 'A', 'A', 'A');	
	uint32_t increment = 4;
	uint32_t result_mask = 0xf;
	uint32_t num_trim_base = 4;
	if(trim_flag != 4){
		char trim_char[4];
		trim_char[0] = 'A';
		trim_char[1] = 'C';
		trim_char[2] = 'G';
		trim_char[3] = 'T';
		for(uint32_t i = trim_flag; i < 3; i ++){
			kmer_pos[i] = kmer_pos[i + 1];
			trim_char[i] = trim_char[i + 1];
		} 
		vkmer_index = _mm512_set_epi32(num_hash_func, kmer_pos[2], kmer_pos[2], kmer_pos[2], kmer_pos[2], kmer_pos[2], kmer_pos[1], kmer_pos[1], 
				kmer_pos[1], kmer_pos[1], kmer_pos[1], kmer_pos[0], kmer_pos[0], kmer_pos[0], kmer_pos[0], kmer_pos[0]);
		vhash_seed_index = _mm512_set_epi32(num_hash_func, 4, 3, 2, 1, 0, 4, 3, 2 ,1, 0, 4, 3, 2, 1, 0);
		vtail = _mm512_set_epi32(0, trim_char[2], trim_char[2], trim_char[2], trim_char[2], trim_char[2], trim_char[1], trim_char[1], trim_char[1], trim_char[1], trim_char[1], 
				trim_char[0], trim_char[0], trim_char[0], trim_char[0], trim_char[0]);	
		vincrement = _mm512_set1_epi32(5);
		increment = 5;
		result_mask = 0x1f;
		num_trim_base = 3;
	}else{
		vkmer_index = _mm512_set_epi32(kmer_pos[3], kmer_pos[3], kmer_pos[3], kmer_pos[3], kmer_pos[2], kmer_pos[2], kmer_pos[2], kmer_pos[2], 
				kmer_pos[1], kmer_pos[1], kmer_pos[1], kmer_pos[1], kmer_pos[0], kmer_pos[0], kmer_pos[0], kmer_pos[0]);
		vincrement = _mm512_set1_epi32(4);
	}
	//printf("%d %d %d %d \n", kmer_pos[0], kmer_pos[1], kmer_pos[2], kmer_pos[3]);


	__m512i vtext, vtext_reverse;
	__m512i vhash_seed_num = _mm512_set1_epi32(num_hash_func - 1);
	vhash_seed = _mm512_set1_epi32(kmer_int_size);
	__mmask16 vtail_mask = _mm512_cmp_epi32_mask(vhash_seed, vkmer_index, _MM_CMPINT_EQ);
	vtail = _mm512_mask_permutevar_epi32(vtail, vtail_mask, vtail, vreverse_char);
	if(tail_flag == 1){
		vtail_mask = _mm512_knot(vtail_mask);
	}
	vtail      = _mm512_mask_slli_epi32(vtail, vtail_mask, vtail, ((kmer_length & 3) - 1) << 3);
	vtail_mask = _mm512_knot(vtail_mask);
	vcontinue  = _mm512_cmp_epi32_mask(vhash_seed_index, vhash_seed_num, _MM_CMPINT_LE);
	vhash_seed = _mm512_mask_i32gather_epi32(vhash_seed, vcontinue, vhash_seed_index, hash_seed, 4);
	__m512i vhash_value;
	__m512i voriginal_index;
	__m512i vBIT = _mm512_set_epi32(0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01, 0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01);
	__m512i vbit_vector_width = _mm512_set1_epi32(bit_vector_width);
	__m512i vbit_vector_index;
	__m512i vbit;
	//__m512i vrem_mask = _mm512_set1_epi32(7);
	__m512i vbit_mask;

	uint32_t result = num_hash_func >= increment ? result_mask : (1 << num_hash_func) - 1;
	uint32_t hash_remainder = num_hash_func;
	pass_result[0] = 1;
	pass_result[1] = 1;
	pass_result[2] = 1;
	pass_result[3] = 1;
	//__m512i four = _mm512_set1_epi32(4);
	__mmask16 vmask;
	/***************************************************/
	__m512i three = _mm512_set1_epi32(3);
	__m512i chbit = _mm512_set1_epi32(0xff);
	//__m512i vone = _mm512_set1_epi32(1);
	__m512i vrem32 = _mm512_set1_epi32(0x1f);
	
	/**************************************************/

	/*
	   uint32_t *tt = (uint32_t*)_mm_malloc(16 * sizeof(uint32_t), 64);
	   printf("tail mask %x \n", vtail_mask);
	//vtail = _mm512_slli_epi32(vtail, 24);
	_mm512_store_epi32(tt, vkmer_index);
	for(uint32_t i = 0; i < 16; i ++){
	printf("%x ", tt[i]);
	}
	printf("\n");
	_mm512_store_epi32(tt, vtail);
	for(uint32_t i = 0; i < 16; i ++){
	printf("%x ", tt[i]);
	}
	printf("\n");
	for(uint32_t i = 0; i < kmer_length; i ++){
	printf("%c", kmer_buffer[36 + i]);
	}
	printf("\n");
	uint32_t *mm = (uint32_t*)kmer_buffer;
	for(uint32_t i = 0; i < 72/4; i ++){
	printf("%x  ", mm[i]);
	}	
	printf("hhhhhh\n");
	*/
	while(1){
		//compute hash value	

		MurmurHash3_x86_32(kmer_buffer, vkmer_index, vcontinue, kmer_length, vhash_seed, vhash_value, vtail, vtail_mask);  
		//voriginal_index = _mm512_rem_epi32(vhash_value, vbit_vector_width);
		//vbit_vector_index = _mm512_srli_epi32(voriginal_index, 3); 
		////		vbit = _mm512_mask_i32extgather_epi32(vbit, vcontinue, vbit_vector_index, bit_vector, _MM_UPCONV_EPI32_UINT8, 1, 0);
		///****************************************************************************/
		////		vbit = _mm512_mask_i32extgather_epi32(vbit, vcontinue, vbit_vector_index, bit_vector, _MM_UPCONV_EPI32_UINT8, 1, 0);
		//__m512i vseq_ind32 = _mm512_srli_epi32(vbit_vector_index, 2);
		//vbit             = _mm512_mask_i32gather_epi32(vbit, vmask, vseq_ind32, bit_vector, 4);
		//__m512i vseq_ind8 = _mm512_and_epi32(vbit_vector_index, three);
		//vseq_ind8 = _mm512_slli_epi32(vseq_ind8, 3);
		//vbit 		= _mm512_srlv_epi32(vbit, vseq_ind8);
		//vbit		= _mm512_and_epi32(vbit, chbit);
		///*****************************************************************/	


		//vbit_mask = _mm512_permutevar_epi32(voriginal_index, vBIT);
		//vbit = _mm512_and_epi32(vbit, vbit_mask);
		//vmask = _mm512_cmp_epi32_mask(vbit, vbit_mask, _MM_CMPINT_EQ);
/*************************************************************************/
		voriginal_index = _mm512_rem_epu32(vhash_value, vbit_vector_width);
		vbit_vector_index = _mm512_srli_epi32(voriginal_index, 5);
		vbit             = _mm512_mask_i32gather_epi32(vbit, vcontinue, vbit_vector_index, bit_vector, 4);
		vbit_mask = _mm512_and_epi32(voriginal_index, vrem32);
		vbit_mask =  _mm512_sllv_epi32(vone, vbit_mask);
		vbit     = _mm512_and_epi32(vbit, vbit_mask);
		vmask = _mm512_cmp_epi32_mask(vbit, vbit_mask, _MM_CMPINT_EQ);
/*************************************************************************/

		uint32_t k = _mm512_mask2int(vmask);
		uint32_t candidate_solid_kmer_mask = 0;
		for(uint32_t i = 0; i < num_trim_base; i ++){
			if((k & result) == result && (pass_result[i] == 1)){
				pass_result[i] = 1;
				candidate_solid_kmer_mask += (result_mask << (increment * i));
			}else{
				pass_result[i] = 0;
			}
			k >>= increment;			
		}	
		hash_remainder -= increment;
		result = hash_remainder >= increment ? result_mask : (1 << hash_remainder) - 1;
		vhash_seed_index = _mm512_add_epi32(vhash_seed_index, vincrement);	
		vcontinue = _mm512_cmp_epi32_mask(vhash_seed_index, vhash_seed_num, _MM_CMPINT_LE);
		__mmask16 vcandidate_solid_kmer_mask = _mm512_int2mask(candidate_solid_kmer_mask);
		vcontinue = _mm512_kand(vcontinue, vcandidate_solid_kmer_mask);
		k = _mm512_mask2int(vcontinue);
		if(k == 0)
			break;
		vhash_seed = _mm512_mask_i32gather_epi32(vhash_seed, vcontinue, vhash_seed_index, hash_seed, 4);

	}

	//_mm512_store_epi32(tt, vhash_value);
	/*	
		for(uint32_t i = 0; i < 4; i ++){
		printf("%x ", pass_result[i]);
		}
		printf("\n");
		*/

	//break;
	//_mm_free(tt);
	if(tail_flag == 0){
		kmer_buffer[0] = save_char;
	}else if(tail_flag == 1){
		kmer_buffer[kmer_length - 1] = save_char;
	}
	//#endif

}

__ONMIC__ bool one_path_query_text(char *kmer_buffer, const uint8_t *bit_vector, const uint32_t *hash_seed, uint32_t kmer_length, uint32_t num_hash_func, uint32_t bit_vector_width){
	//#ifdef __MIC__	
	int comparison_result = 0;
	uint32_t kmer_int_size = KMER_MAX_LEN / 4;	
	uint32_t kmer_pos = 0;
	for (unsigned int it = 0, it_rc = (kmer_length - 1); it < kmer_length;
			it++, it_rc--) {
		if (kmer_buffer[it] == 'A') {
			if (kmer_buffer[it_rc] == 'T') {
			} else {
				comparison_result = 1;
				break;
			}
		} else if (kmer_buffer[it] == 'C') {
			if (kmer_buffer[it_rc] == 'G') {
			} else if (kmer_buffer[it_rc] == 'T') {
				comparison_result = -1;
				break;
			} else {
				comparison_result = 1;
				break;
			}
		} else if (kmer_buffer[it] == 'G') {
			if (kmer_buffer[it_rc] == 'A') {
				comparison_result = 1;
				break;
			} else if (kmer_buffer[it_rc] == 'C') {
			} else {
				comparison_result = -1;
				break;
			}
		} else if (kmer_buffer[it] == 'T') {
			if (kmer_buffer[it_rc] == 'A') {
			} else {
				comparison_result = -1;
				break;
			}
		} else {
			//std::cout << std::endl << "ERROR: Illegal character " << kmer[it] << " (query_text)" << std::endl << std::endl;
			printf("ERROR: 1 Illegal charactor %c query_text\n\n", kmer_buffer[it]);
		}
	}
	//printf("compare result %d \n", comparison_result);

	//	__m512i vreverse_char   = _mm512_set_epi32(0, 0, 0, 0, 0, 0, 0, 0, 0x43, 0, 0, 0x41, 0x47, 0, 0x54, 0);
	//
	//	__m512i vkmer_length = _mm512_set1_epi32(kmer_length);
	//	__m512i vkmer_index = _mm512_set_epi32(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
	//	__m512i vkmer_reverse_index = _mm512_set1_epi32(KMER_MAX_LEN + kmer_length - 1);
	//	vkmer_reverse_index = _mm512_sub_epi32(vkmer_reverse_index, vkmer_index);
	//	__mmask16 vcontinue = _mm512_cmp_epi32_mask(vkmer_index, vkmer_length, _MM_CMPINT_LT);
	//	__m512i vkmer, vkmer_reverse;
	//	__m512i vincrement = _mm512_set1_epi32(16);
	//    if(comparison_result == -1){
	//        kmer_pos = kmer_int_size;
	//	    for(uint32_t i = 0; i < ((kmer_length + 15) >> 4); i ++){
	//		    vkmer = _mm512_mask_i32extgather_epi32(vkmer, vcontinue, vkmer_index, kmer_buffer, _MM_UPCONV_EPI32_UINT8, 1, 0);
	//		    vkmer_reverse = _mm512_permutevar_epi32(vkmer, vreverse_char);
	//		    _mm512_mask_i32extscatter_epi32(kmer_buffer, vcontinue, vkmer_reverse_index, vkmer_reverse, _MM_DOWNCONV_EPI32_UINT8, 1, 0);
	//		    vkmer_index = _mm512_add_epi32(vkmer_index, vincrement);
	//		    vkmer_reverse_index = _mm512_sub_epi32(vkmer_reverse_index, vincrement);
	//		    vcontinue = _mm512_cmp_epi32_mask(vkmer_index, vkmer_length, _MM_CMPINT_LT);
	//	    }
	//    }
	//
	//
	/********************************************************************************/
	//__m512i vreverse_char   = _mm512_set_epi32(0, 0, 0, 0, 0, 0, 0, 0, 0x43, 0, 0, 0x41, 0x47, 0, 0x54, 0);
	//__m512i vkmer_length = _mm512_set1_epi32(kmer_length);
	//__m512i vkmer_index8 = _mm512_set_epi32(0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60);
	//__m512i vkmer_index = _mm512_set_epi32(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);


	//__m512i vkmer_reverse_index = _mm512_set1_epi32(KMER_MAX_LEN + kmer_length - 1);

	//vkmer_reverse_index = _mm512_sub_epi32(vkmer_reverse_index, vkmer_index8);
	//__mmask16 vcontinue = _mm512_cmp_epi32_mask(vkmer_index8, vkmer_length, _MM_CMPINT_LT);
	//__m512i vkmer, vkmer_reverse;
	//__m512i vincrement = _mm512_set1_epi32(16);
	//if(comparison_result == -1){
	//	__m512i vone = _mm512_set1_epi32(1);
	//	__m512i vzero = _mm512_set1_epi32(0);
	//	__m512i vsixty = _mm512_set1_epi32(60);
	//	__m512i vsixtyfour = _mm512_set1_epi32(64);
	//	vkmer_reverse_index = _mm512_srli_epi32(vkmer_reverse_index, 2);
	//	for(int i = 0; i < kmer_length; i += 64){
	//		vkmer = _mm512_mask_i32gather_epi32(vkmer, vcontinue, vkmer_index, kmer_buffer, 4);
	//		for(int j = 0; j < 4; j ++){
	//			__m512i vkmer_re8 = _mm512_mask_permutevar_epi32(vkmer_reverse, vcontinue, vkmer, vreverse_char); 
	//			vkmer_reverse = _mm512_mask_slli_epi32(vkmer_reverse, vcontinue, vkmer_reverse, 8); 
	//			vkmer_reverse = _mm512_mask_or_epi32(vkmer_reverse, vcontinue, vkmer_reverse, vkmer_re8);
	//			vkmer_index8 = _mm512_add_epi32(vkmer_index8, vone);
	//			vcontinue = _mm512_cmp_epi32_mask(vkmer_index8, vkmer_length, _MM_CMPINT_LT);
	//		}

	//		vcontinue = _mm512_cmp_epi32_mask(vkmer_reverse_index, vzero, _MM_CMPINT_GE);
	//		_mm512_mask_i32scatter_epi32(kmer_buffer, vcontinue, vkmer_reverse_index, vkmer_reverse, 4);
	//		vkmer_reverse_index = _mm512_sub_epi32(vkmer_reverse_index, vincrement);

	//		vkmer_index = _mm512_add_epi32(vkmer_index, vincrement);
	//		vkmer_index8 = _mm512_add_epi32(vkmer_index8, vsixty);
	//		vcontinue = _mm512_cmp_epi32_mask(vkmer_index8, vkmer_length, _MM_CMPINT_LT);
	//	}		

	//}
	__m512i vreverse_char   = _mm512_set_epi32(0, 0, 0, 0, 0, 0, 0, 0, 0x43, 0, 0, 0x41, 0x47, 0, 0x54, 0);
	__m512i vkmer_length = _mm512_set1_epi32(kmer_length);
	__m512i vkmer_index8 = _mm512_set_epi32(0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60);
	__m512i vkmer_index = _mm512_set_epi32(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);


	__m512i vkmer_reverse_index = _mm512_set1_epi32(KMER_MAX_LEN + kmer_length - 1);

	vkmer_reverse_index = _mm512_sub_epi32(vkmer_reverse_index, vkmer_index8);
	__mmask16 vcontinue = _mm512_cmp_epi32_mask(vkmer_index8, vkmer_length, _MM_CMPINT_LT);
	__m512i vkmer, vkmer_reverse;
	__m512i vincrement = _mm512_set1_epi32(16);

	if(comparison_result == -1){

	        kmer_pos = kmer_int_size;
		__m512i vone = _mm512_set1_epi32(1);
		__m512i vbase_addr = _mm512_set1_epi32(KMER_MAX_LEN>>2);
		__m512i vsixty = _mm512_set1_epi32(60);
		__m512i vsixtyfour = _mm512_set1_epi32(64);
		vkmer_reverse_index = _mm512_srli_epi32(vkmer_reverse_index, 2);

		for(int m = 0; m < kmer_length; m += 64){
			vkmer = _mm512_mask_i32gather_epi32(vkmer, vcontinue, vkmer_index, kmer_buffer, 4);

			vkmer_reverse = _mm512_set1_epi32(0);

			for(int j = 0; j < 4; j ++){
				__m512i vkmer_re8 = _mm512_mask_permutevar_epi32(vkmer_reverse, vcontinue, vkmer, vreverse_char); 
				vkmer_reverse = _mm512_mask_slli_epi32(vkmer_reverse, vcontinue, vkmer_reverse, 8); 
				vkmer_reverse = _mm512_mask_or_epi32(vkmer_reverse, vcontinue, vkmer_reverse, vkmer_re8);
				vkmer = _mm512_mask_srli_epi32(vkmer, vcontinue, vkmer, 8); 
			}
			if((kmer_length & 3) != 0){

				__m512i vpermt = _mm512_set_epi32(16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1);

				vkmer = _mm512_permutevar_epi32(vpermt, vkmer_reverse); 
				vkmer = _mm512_slli_epi32(vkmer, (kmer_length & 3) << 3);
				vkmer_reverse = _mm512_srli_epi32(vkmer_reverse, (4 - (kmer_length & 3)) << 3); 
				vkmer_reverse = _mm512_or_epi32(vkmer_reverse, vkmer);

			}

			vcontinue = _mm512_cmp_epi32_mask(vkmer_reverse_index, vbase_addr, _MM_CMPINT_GE);
			_mm512_mask_i32scatter_epi32(kmer_buffer, vcontinue, vkmer_reverse_index, vkmer_reverse, 4);
			vkmer_reverse_index = _mm512_sub_epi32(vkmer_reverse_index, vincrement);

			vkmer_index = _mm512_add_epi32(vkmer_index, vincrement);
			vkmer_index8 = _mm512_add_epi32(vkmer_index8, vsixty);
			vcontinue = _mm512_cmp_epi32_mask(vkmer_index8, vkmer_length, _MM_CMPINT_LT);
		}		

	}

	//	    /********************************************************************************/



	/************************************************************************************************/

	 
	   //uint32_t kk = 0;
	   //if(comparison_result < 0)
	   //kk = KMER_MAX_LEN;
	   //for(int i = 0; i < 31; i ++)
	   //printf("%c", kmer_buffer[kk + i]); 
	   //printf("\n");
	    
	__m512i vhash_seed_index = _mm512_set_epi32(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1,0); 
	__m512i vhash_seed;// = _mm512_set1_epi32(hash_seed[0]);
	__m512i vhash_num = _mm512_set1_epi32(num_hash_func - 1);
	__m512i vhash_value;
	__m512i vBIT = _mm512_set_epi32(0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02,
			0x01, 0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01);
	__m512i vbit_vector_width = _mm512_set1_epi32(bit_vector_width);
	__m512i voriginal_index;
	__m512i vbit_vector_index;
	__m512i vbit;
	//__m512i vrem_mask = _mm512_set1_epi32(7);
	__m512i vbit_mask;
	__mmask16 vmask; //vmask_hash_region, vsolid_kmer_mask, vhash_seed_reset_mask, 
	//vhash_seed_increment_mask, vregions_index_increment_mask;
	//vcontinue = _mm512_int2mask(0xffff);
	vcontinue = _mm512_cmp_epi32_mask(vhash_seed_index, vhash_num, _MM_CMPINT_LE);		
	vkmer_index = _mm512_set1_epi32(kmer_pos);
	uint32_t result = 0;
	vhash_seed = _mm512_mask_i32gather_epi32(vhash_seed, vcontinue, vhash_seed_index, hash_seed, 4);
	uint32_t k, result_continue;
	/***************************************************/
	__m512i three = _mm512_set1_epi32(3);
	__m512i chbit = _mm512_set1_epi32(0xff);
	__m512i vone = _mm512_set1_epi32(1);
	__m512i vrem32 = _mm512_set1_epi32(0x1f);
	/**************************************************/


	while(1){
		one_MurmurHash3_x86_32(kmer_buffer, vkmer_index, vcontinue, kmer_length, vhash_seed, vhash_value);  

		//uint32_t *tmp1 = (uint32_t*)_mm_malloc(16 * sizeof(int), 64);
		//_mm512_store_epi32(tmp1, vhash_value);
		//for(int i = 0; i < 16; i ++)
		//	printf("%x ", tmp1[i]);
		//printf("\n");

		/*    
		      int *tmp1 = (int*)_mm_malloc(16 * sizeof(int), 64);
		      _mm512_store_epi32(tmp1, vhash_value);
		      for(int i = 0; i < 16; i ++)
		      printf("%x ", tmp1[i]);
		      printf("\n");
		      */

		/*	
			_mm512_store_epi32(ttt, vhash_value);
			for(uint32_t i = 0; i < 16; i ++){
			printf("%x ", ttt[i]);
			}
			printf("\n %x \n", vcontinue);

			hh += 1;
			if(hh == 6)
			break;

			_mm512_store_epi32(ttt, vhash_value);
			for(uint32_t i = 0; i < 16; i ++){
			printf("%x ", ttt[i]);
			}
			printf("\n");
			*/
		//voriginal_index = _mm512_rem_epi32(vhash_value, vbit_vector_width);
		//vbit_vector_index = _mm512_srli_epi32(voriginal_index, 3); 
		////		vbit = _mm512_mask_i32extgather_epi32(vbit, vcontinue, vbit_vector_index, bit_vector, _MM_UPCONV_EPI32_UINT8, 1, 0);
		///****************************************************************************/
		////		vbit = _mm512_mask_i32extgather_epi32(vbit, vcontinue, vbit_vector_index, bit_vector, _MM_UPCONV_EPI32_UINT8, 1, 0);
		//__m512i vseq_ind32 = _mm512_srli_epi32(vbit_vector_index, 2);
		//vbit             = _mm512_mask_i32gather_epi32(vbit, vmask, vseq_ind32, bit_vector, 4);
		//__m512i vseq_ind8 = _mm512_and_epi32(vbit_vector_index, three);
		//vseq_ind8 = _mm512_slli_epi32(vseq_ind8, 3);
		//vbit 		= _mm512_srlv_epi32(vbit, vseq_ind8);
		//vbit		= _mm512_and_epi32(vbit, chbit);
		///*****************************************************************/	

		//vbit_mask = _mm512_permutevar_epi32(voriginal_index, vBIT);
		//vbit = _mm512_and_epi32(vbit, vbit_mask);
		//vmask = _mm512_cmp_epi32_mask(vbit, vbit_mask, _MM_CMPINT_EQ);
		
/*************************************************************************/
		voriginal_index = _mm512_rem_epu32(vhash_value, vbit_vector_width);
		vbit_vector_index = _mm512_srli_epi32(voriginal_index, 5);
		vbit             = _mm512_mask_i32gather_epi32(vbit, vcontinue, vbit_vector_index, bit_vector, 4);
		vbit_mask = _mm512_and_epi32(voriginal_index, vrem32);
		vbit_mask =  _mm512_sllv_epi32(vone, vbit_mask);
		vbit     = _mm512_and_epi32(vbit, vbit_mask);
		vmask = _mm512_cmp_epi32_mask(vbit, vbit_mask, _MM_CMPINT_EQ);
/*************************************************************************/


		vmask = _mm512_kand(vcontinue, vmask);		
		
		k = _mm512_mask2int(vmask);
		result_continue = _mm512_mask2int(vcontinue);
		vhash_seed_index = _mm512_add_epi32(vhash_seed_index, vincrement);
		vcontinue = _mm512_cmp_epi32_mask(vhash_seed_index, vhash_num, _MM_CMPINT_LE);
		if((k != result_continue) || (vcontinue == 0))
			break;		
		vhash_seed = _mm512_mask_i32gather_epi32(vhash_seed, vcontinue, vhash_seed_index, hash_seed, 4);
	}	
	if(k == result_continue)
		return true;	
	//#endif
	return false;
}

