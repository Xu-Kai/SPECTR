/*
 * correct_errors.h
 *
 *  Created on: 2015年12月14日
 *      Author: xk
 */

#ifndef MIC_CORRECT_ERRORS_H_
#define MIC_CORRECT_ERRORS_H_
#include "define.h"
#include "parse_args.h"
#include <stdint.h>

 void correct_errors_on_mic(init_args arg,C_arg& c_inst_arg,uint8_t* bit_vector);
__ONMIC__ void mic_correct_errors(char* read_content, char* read_quality, char* sequence_modification,
        uint32_t* hash_seed, uint8_t* bit_vector, init_args& arg,
        uint32_t& num_corrected_errors_local,uint32_t& trim_5_end, uint32_t& trim_3_end,
        uint32_t max_trimmed_bases,uint32_t min_check_length, uint32_t read_length, int ssss, func_mem& funcm);
	


#endif /* CORRECT_ERRORS_H_ */
