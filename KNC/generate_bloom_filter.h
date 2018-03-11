/*
 * generate_bloom_filter.h
 *
 *  Created on: 2016年2月27日
 *      Author: xk
 */

#ifndef GENERATE_BLOOM_FILTER_H_
#define GENERATE_BLOOM_FILTER_H_
#include "define.h"
#include "parse_args.h"
#include "time.h"
#include "kmc_file.h"
#include "new_query_text.h"
#include "MurmurHash3.h"
//----------------------------------------------------------------------
// C_correct_errors
//----------------------------------------------------------------------
class C_generate_bloom {
public:
   // variables
   bloom_type num_unique_solid_kmers;
   bloom_type random_seed;

   std::size_t bit_vector_width;
   std::size_t bit_vector_width_byte;
   std::size_t extremely_low_quality_score;
   std::size_t kmer_occurrence_threshold;
   std::size_t max_extension;
   std::size_t num_reads;
   std::size_t num_trimmed_bases;
   std::size_t quality_score_cutoff;
   std::size_t quality_score_offset;


   unsigned int kmer_length;

   int rank_node;
   int rank_smp;
   int size_node;
   int group_rank;

   unsigned short int num_hash_func;
   unsigned short int num_hash_func_real;


   // constructors
   C_generate_bloom() :
                       num_unique_solid_kmers(0),
                       random_seed(0),
                       bit_vector_width(0),
                       bit_vector_width_byte(0),
                       extremely_low_quality_score(0),
                       kmer_occurrence_threshold(0),
                       max_extension(0),
                       num_reads(0),
                       num_trimmed_bases(0),
                       quality_score_cutoff(0),
                       quality_score_offset(0),
                       kmer_length(0),
                       rank_node(-1),
                       rank_smp(-1),
                       size_node(-1),
                       num_hash_func(0),
                       num_hash_func_real(0),
                       num_corrected_errors(0),
                       num_corrected_reads(0),
                       read_block_size(0)
                      {};


   void determine_bloom_filter_parameters(const C_arg& c_inst_args, C_time& c_inst_time);
   void program_kmers(init_args& arg, const C_arg& c_inst_args, C_time& c_inst_time, MPI_Comm& group_comm);
   void read_bloom_filter_parameters(const C_arg& c_inst_args, C_time& c_inst_time);
   void generate_hash_seed(const bloom_type random_seed, unsigned int** hash_seed);
private:

   std::size_t num_corrected_errors;
   std::size_t num_corrected_reads;
   std::size_t read_block_size;

   std::string rank_node_text;
   std::string kmc_prefix;


};


#endif /* GENERATE_BLOOM_FILTER_H_ */
