#ifndef _COUNT_SOLID_KMERS2_H
#define _COUNT_SOLID_KMERS2_H



#include "parse_args.h"
#include "check_inputs.h"
#include "kmc_file.h"
#include "define.h"

   void init_count_solid_kmers(C_count_solid_kmers2& c_inst_count_solid_kmers2);
   void count_kmers2(const C_arg& c_inst_args, C_time& c_inst_time, C_count_solid_kmers2& c_inst_count_solid_kmers2);

   //functions
   void run_kmc(const C_arg& c_inst_args, C_count_solid_kmers2& c_inst_count_solid_kmers2);
   void determine_kmer_occurrence_threshold(const C_arg& c_inst_args, C_count_solid_kmers2& c_inst_count_solid_kmers2);
   void generate_kmer_occurrence_histogram(const C_arg& c_inst_args, C_count_solid_kmers2& c_inst_count_solid_kmers2);
   void count_unique_solid_kmers(const C_arg& c_inst_args, C_count_solid_kmers2& c_inst_count_solid_kmers2);
   void write_kmer_histogram(const C_arg& c_inst_args, C_count_solid_kmers2& c_inst_count_solid_kmers2);



#endif
