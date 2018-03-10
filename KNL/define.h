/*
 * define.h
 *
 *  Created on: 2015年12月8日
 *      Author: xk
 */

#ifndef DEFINE_H_
#define DEFINE_H_

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<time.h>
#include <stdint.h>
#include <mpi.h>
// openmp
#include <omp.h>
#include <unistd.h>

#include<sys/stat.h>
#include <sys/types.h>

#include "kmc_file.h"

// definitions
#define QS_HISTOGRAM_MAX          126
#define BITS_PER_CHAR             8
#define NUM_NEOCLEOTIDE           4
#define A                         0
#define C                         1
#define G                         2
#define T                         3
#define DEFAULT_SEED              0
#define DEFAULT_FPR               0.001
#define MAX_EXTENSION             5
#define MIN_MAX_MEM               4
#define MAX_LOW_QS_BASES          3
#define MIN_KMER_LENGTH           10
#define MAX_KMER_LENGTH           255 // max value of KMC
#define MAX_KMER_THRESHOLD        65535
#define MIN_MAX_MEM               4
#define MAX_MAX_MEM               1024
#define PHRED33                   33
#define PHRED64                   64
#define MAX_PHRED                 41
#define QS_CUTOFF_RATIO           0.05
#define QS_EXTREMELY_LOW_RATIO    0.01
#define MMAP_FILE_SIZE            209715200 // 200 * 1024 * 1024 = 200 MB
#define KMC_DEFAULT_MIN_COUNT     2
#define KMER_HISTOGRAM_SIZE       256 // KMC default
#define KMC_DIR                   "kmc/bin"
#define KMC_BINARY                "kmc"
#define PIGZ_DIR                  "pigz/pigz-2.3.3"
#define PIGZ_BINARY               "pigz"
#define KMER_BLOCK_SIZE           100000
#define OPENMP_CHUNK_SIZE         10
#define READ_BLOCK_SIZE_RATIO     1000
#define BLOOM_RCV_BUFFER_SIZE     104857600 // 100 * 1024 * 1024 = 100 MB
#define CHECK_RANGE_RATIO         0.07
#define MAX_N_RATIO               0.1
#define MIN_BASES_AFTER_TRIMMING  30
#define MAX_TRIMMING_RATE         0.6
#define SUBST_CHAR                'A'
#define MAX_ERROR_RATE            0.25
#define BIT_VEC_INC_RATIO         1.2
#define MIN_SOLID_LENGTH          2
#define MIN_NON_SOLID_LENGTH      2
#define FP_SUSPECT_LENGTH         1
#define SOLID_REGION_ADJUST_RANGE 4
#define MAX_CANDIDATE_PATHS       74
#define NUM_ALLOWABLE_FAILS       2
#define INIT_MIN_QS               1000000
#define MAX_MODIFICATION          4
#define MIN_QS_DIFF               10


//#define __ONMIC__ __attribute__((target(mic)))
#define __ONMIC__ 

#define bloom_type uint32_t
//#pragma offload_attribute (push,target(mic))
 __ONMIC__ static const char NEOCLEOTIDE[NUM_NEOCLEOTIDE] = {'A', 'C', 'G', 'T'};

 __ONMIC__ static const unsigned char BIT_MASK[BITS_PER_CHAR] = {
                                                      0x01, //00000001
                                                      0x02, //00000010
                                                      0x04, //00000100
                                                      0x08, //00001000
                                                      0x10, //00010000
                                                      0x20, //00100000
                                                      0x40, //01000000
                                                      0x80  //10000000
                                                     };



// __ONMIC__ typedef unsigned int bloom_type;

__ONMIC__ typedef struct temporary_memory{
	char* kmer;
	uint32_t kmer_size;
	uint32_t kmer_ind;
	uint32_t  kmer_length;
	char* sequence;
	uint32_t sequence_size;
	uint32_t sequence_ind;
	uint32_t candidate_path_size;
	uint32_t candidate_path_ind;
	uint32_t* kmer_index;
	uint32_t *ubuffer;
	char *cbuffer;
	uint32_t ubuffer_offset;
	uint32_t cbuffer_offset;
	uint32_t *ukmer_buffer;
}func_mem;

 __ONMIC__ typedef struct init_argsss{
	uint32_t  bit_vector_width ;
	uint32_t  bit_vector_width_byte ;
	uint32_t  num_hash_func;
	uint32_t  kmer_length ;
	uint32_t  kmer_occurrence_threshold ;
	uint32_t  random_seed ;
	uint32_t  num_unique_solid_kmers;
	uint32_t read_len;
	uint32_t quality_score_offset;
	uint32_t quality_score_cutoff;
	uint32_t extremely_low_quality_score;
	uint32_t max_extension;
	uint32_t read_num;
    uint32_t* hash_seed;
    uint8_t* bit_vector;
	bool notrim;
} init_args;
 __ONMIC__ typedef struct computer_result{
    uint64_t num_corrected_errors;
    uint64_t num_corrected_reads;
    uint64_t num_trimmed_bases; 
} stats_info;


//#pragma offload_attribute (pop)
typedef struct c_const_file{
    char* read_list;
    char **read_file_name;
    uint32_t read_list_size;
    uint32_t read_file_num;
    uint32_t *read_group;
}C_file_name;
typedef struct Count_solid_kmers{
    uint64_t num_unique_solid_kmers;
    uint64_t *num_occurrences_histogram;
    uint64_t *reduced_num_occurrences_histogram;
    uint64_t kmer_occurrence_threshold;
    uint64_t valey_point;
    uint32_t kmer_length;
    int rank_node;
    int size_node;
    bool get_valey_point;
    CKMCFile kmer_database;
    std::string rank_node_text;
    std::string kmc_prefix;
}C_count_solid_kmers2;
#endif /* DEFINE_H_ */
