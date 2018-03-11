#ifndef _PARSE_READ_H
#define _PARSE_READ_H



#include "parse_args.h"
#include "time.h"



//----------------------------------------------------------------------
// C_check_read
//----------------------------------------------------------------------
class C_check_read {
public:
   // variables
   std::size_t extremely_low_quality_score;
   std::size_t init_num_elements;
   std::size_t interval;
   std::size_t longest_line_length;
   std::size_t min_num_processed_reads;
   std::size_t num_reads;
   std::size_t quality_score_cutoff;
   std::size_t quality_score_offset;

   long long int read_file_size_byte;
   long long int read_file_size_byte1;
   long long int read_file_size_byte2;

   int max_quality_score;
   int min_quality_score;
   int read_file_unit_size_byte;
   int read_file_unit_size_byte1;
   int read_file_unit_size_byte2;
   int size_node;

   std::vector<std::size_t> num_reads_vector;
   std::vector<std::size_t> num_reads_vector1;
   std::vector<std::size_t> num_reads_vector2;
   std::vector<std::size_t> starting_point_vector;
   std::vector<std::size_t> starting_point_vector1;
   std::vector<std::size_t> starting_point_vector2;
   uint64_t*  qs_histo_vec;

   // constructor
   C_check_read() :
                   extremely_low_quality_score(0),
                   init_num_elements(0),
                   interval(0),
                   longest_line_length(0),
                   min_num_processed_reads(0),
                   num_reads(0),
                   quality_score_cutoff(0),
                   quality_score_offset(0),
                   read_file_size_byte(0),
                   read_file_size_byte1(0),
                   read_file_size_byte2(0),
                   max_quality_score(-1000),
                   min_quality_score(1000),
                   read_file_unit_size_byte(0),
                   read_file_unit_size_byte1(0),
                   read_file_unit_size_byte2(0),
                   size_node(-1)
                  {};

   // functions
   void check_read_file(const C_arg& c_inst_args, C_time& c_inst_time, int group_rank, MPI_Comm &group_comm, int node_num);

private:
   // variables

   // functions
   void check_read_file_fastq_paired(const C_arg& c_inst_args, int group_rank, MPI_Comm& group_comm, int node_num);
   void check_read_file_fastq_single(const C_arg& c_inst_args);
   void construct_quality_score_histogram_single_unzipped(const C_arg& c_inst_args);
   void construct_quality_score_histogram_single_gzipped(const C_arg& c_inst_args);
   void construct_quality_score_histogram_paired_unzipped(const C_arg& c_inst_args, int group_rank, MPI_Comm& group_comm, int node_num);
   void construct_quality_score_histogram_paired_gzipped(const C_arg& c_inst_args);
};



#endif
