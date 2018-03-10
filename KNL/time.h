#ifndef _TIME_LOCAL_H
#define _TIME_LOCAL_H

#include "define.h"
#include<string>

//----------------------------------------------------------------------
// C_time
//----------------------------------------------------------------------
class C_time {
public:
   // variables
   // parse arguments
   std::string start_parse_args;
   std::string end_parse_args;

   // check read files
   std::string start_check_read_file;
   std::string end_check_read_file;

   // count k-mers
   std::string start_count_kmers;
   std::string end_count_kmers;

   // program k-mers into the bloom filter
   std::string start_program_kmers_into_bloom_filter;
   std::string end_program_kmers_into_bloom_filter;

   // correct errors in reads
   std::string start_correct_errors_in_reads;
   std::string end_correct_errors_in_reads;

   // merge output files
   std::string start_merge_output_files;
   std::string end_merge_output_files;

   // constructors
   C_time() {};
};



#endif
