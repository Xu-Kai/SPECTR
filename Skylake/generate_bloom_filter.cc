/*
 * generate_bloom_filter.cc
 *
 *  Created on: 2016年2月27日
 *      Author: xk
 */

#include "generate_bloom_filter.h"
#include "define.h"
//----------------------------------------------------------------------
// determine_bloom_filter_parameters
//----------------------------------------------------------------------
void C_generate_bloom::determine_bloom_filter_parameters(const C_arg& c_inst_args, C_time& c_inst_time) {
   if (rank_node == 0) {
      std::cout << "Determining Bloom filter parameters" << std::endl;
   }

   // find optimal parameters of for this Bloom filter
   // initialize variables
   double min_bit_vector_width_element = std::numeric_limits<double>::infinity();
   double min_num_hash_func = 1.0;
   double current_bit_vector_width_element = 0.0;

   // find the number of hash functions for minimum bit-vector size
   for (double j = 1.0; j < 1000.0; ++j) {
      if ((current_bit_vector_width_element = ((- j * num_unique_solid_kmers) / std::log(1.0 - std::pow(c_inst_args.target_false_positive_prob, 1.0 / j)))) < min_bit_vector_width_element) {
         min_bit_vector_width_element = current_bit_vector_width_element;
         min_num_hash_func = j;
      }
   }

   // find the tradeoff point between the bit vector size and the number of hash functions
   for (double j = 1.0; j < min_num_hash_func; ++j) {
      if ((current_bit_vector_width_element = ((- j * num_unique_solid_kmers) / std::log(1.0 - std::pow(c_inst_args.target_false_positive_prob, 1.0 / j)))) < (min_bit_vector_width_element * BIT_VEC_INC_RATIO)) {
         min_bit_vector_width_element = current_bit_vector_width_element;
         min_num_hash_func = j;
         break;
      }
   }

   // determine the size of the bloom filter
   // the 128-bit murmur hash function can generate two 64-bit outputs at the same time
   num_hash_func = static_cast<unsigned int>(min_num_hash_func);

   // make num_hash_func an even number
   if ((num_hash_func % 2) == 1) {
      num_hash_func++;
   }
   num_hash_func_real = num_hash_func / 2;

    
   bit_vector_width      = static_cast<bloom_type>(min_bit_vector_width_element);
   bit_vector_width *= 2;
   bit_vector_width      += (((bit_vector_width % BITS_PER_CHAR) != 0) ? (BITS_PER_CHAR - (bit_vector_width % BITS_PER_CHAR)) : 0);
   bit_vector_width_byte = bit_vector_width / BITS_PER_CHAR;

   bloom_type bit_vector_size_mb(bit_vector_width_byte / (1024 * 1024));
   if (bit_vector_size_mb < 1) {
      bit_vector_size_mb = 1;
   }

   if (rank_node == 0) {
      std::cout << "     Bloom filter" << std::endl;
      std::cout << "     Number of keys           : " << num_unique_solid_kmers << std::endl;
      std::cout << "     Bit-vector size          : " << bit_vector_size_mb << " MB" << std::endl;
      std::cout << "     Number of hash functions : " << num_hash_func << std::endl;
      std::cout << "     Determining Bloom filter parameters: done" << std::endl;
      std::cout << std::endl;
   }
}

//----------------------------------------------------------------------
// generate_hash_seed
//----------------------------------------------------------------------
void C_generate_bloom::generate_hash_seed(const bloom_type random_seed, unsigned int** hash_seed) {
   //hash_seed.resize(num_hash_func_real * 2);
    *hash_seed = (unsigned int*)malloc(num_hash_func_real * 2 * sizeof(unsigned int));
   srand(static_cast<unsigned int>(random_seed));

   for (unsigned short int it_seed = 0; it_seed < num_hash_func_real * 2; it_seed++) {
      (*hash_seed)[it_seed] = static_cast<unsigned int>(rand());
   }
}
//----------------------------------------------------------------------
// program_kmers
//----------------------------------------------------------------------
void C_generate_bloom::program_kmers(init_args& arg, const C_arg& c_inst_args, C_time& c_inst_time, MPI_Comm& group_comm) {
   time_t rawtime;
   time(&rawtime);
   c_inst_time.start_program_kmers_into_bloom_filter = asctime(localtime(&rawtime));

   if (rank_node == 0) {
      std::cout << "Programming k-mers to the Bloom filter" << std::endl;
   }

   CKMCFile kmer_database;

   //--------------------------------------------------
   // open the database
   //--------------------------------------------------
   // set rank_node_text
   std::stringstream rank_node_stream;
   rank_node_stream << std::setw(5) << std::setfill('0') << rank_node;
   rank_node_text = rank_node_stream.str();

   kmc_prefix = c_inst_args.kmc_prefix + "." + rank_node_text;

   // open the database
   if (!kmer_database.OpenForListing(kmc_prefix)) {
      std::cout << std::endl << "ERROR: Cannot open " << kmc_prefix << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 200);
   }

   //--------------------------------------------------
   // generate a temporary bit vector
   //--------------------------------------------------
   // generate unique hash seeds
   //std::vector<unsigned int> hash_seed;
   //printf("start hash seed \n");
   uint32_t* hash_seed;
   generate_hash_seed(c_inst_args.random_seed, &hash_seed);
    //for(int i = 0; i < 6; i ++){
      //  printf("%d\n", hash_seed[i]);
    //}
   // generate a bloom filter
   //std::vector<unsigned char> bit_vector(bit_vector_width_byte, 0);
    uint8_t* bit_vector = (uint8_t*)malloc(bit_vector_width_byte * sizeof(uint8_t));
    memset(bit_vector, 0, bit_vector_width_byte * sizeof(uint8_t));

   //--------------------------------------------------
   // iterate the database
   //--------------------------------------------------
   CKmerAPI kmer((uint32)kmer_length);

   float num_occurrences;

   std::vector<CKmerAPI> kmer_vector;

   std::size_t kmer_vector_index(0);

   std::string kmer_string;

   bloom_type original_index1;
   bloom_type original_index2;
   bloom_type hash [2];

//   uint64_t original_index1;
//   uint64_t original_index2;
//   uint64_t hash [2];

   unsigned short int bit_index1;
   unsigned short int bit_index2;

   kmer_vector.resize(KMER_BLOCK_SIZE);
   // iterate k-mers
   while (kmer_database.ReadNextKmer(kmer, num_occurrences)) {
      if (num_occurrences >= kmer_occurrence_threshold) {
         // add the k-mer to the vector
         kmer_vector[kmer_vector_index] = kmer;

         kmer_vector_index++;

         // kmer_vector is full
         if (kmer_vector_index == KMER_BLOCK_SIZE) {
            //--------------------------------------------------
            // correct reads in a block
            //--------------------------------------------------
            #pragma omp parallel shared(bit_vector, hash_seed, kmer_vector) private(kmer_string, original_index1, original_index2, hash, bit_index1, bit_index2)
            {
               // iterate reads
               #pragma omp for schedule(dynamic, OPENMP_CHUNK_SIZE)
               for (std::size_t it_kmer = 0; it_kmer < KMER_BLOCK_SIZE; it_kmer++) {
                  // convert the k-mer in the db to a string
                  (kmer_vector[it_kmer]).to_string(kmer_string);
                  kmer_string = kmer_string.substr(0, kmer_length);

                  // program the k-mer into the bloom filter
                  for (unsigned short int it_hash_func = 0; it_hash_func < num_hash_func_real * 2; it_hash_func++) {
                     MurmurHash3_x86_32(kmer_string.c_str(), kmer_length, hash_seed[it_hash_func], hash);

                     original_index1 = hash[0] % bit_vector_width;
//                     original_index2 = hash[1] % bit_vector_width;
                     bit_index1      = original_index1 % BITS_PER_CHAR;
//                     bit_index2      = original_index2 % BITS_PER_CHAR;

                     #pragma omp atomic
                     bit_vector[original_index1 / BITS_PER_CHAR] |= BIT_MASK[bit_index1];

//                     #pragma omp atomic
//                     bit_vector[original_index2 / BITS_PER_CHAR] |= BIT_MASK[bit_index2];
                  }
               }
            }

            kmer_vector_index = 0;
         }
      }
   }

   // program remaining k-mers
   if (kmer_vector_index > 0) {
      #pragma omp parallel shared(bit_vector, hash_seed, kmer_vector, kmer_vector_index) private(kmer_string, original_index1, original_index2, hash, bit_index1, bit_index2)
      {
         // iterate reads
         #pragma omp for schedule(dynamic, OPENMP_CHUNK_SIZE)
         for (std::size_t it_kmer = 0; it_kmer < kmer_vector_index; it_kmer++) {
            // convert the k-mer in the db to a string
            kmer_vector[it_kmer].to_string(kmer_string);
            kmer_string = kmer_string.substr(0, kmer_length);

            // program the k-mer into the bloom filter
            for (unsigned short int it_hash_func = 0; it_hash_func < num_hash_func_real * 2; it_hash_func++) {
               MurmurHash3_x86_32(kmer_string.c_str(), kmer_length, hash_seed[it_hash_func], hash);

               original_index1 = hash[0] % bit_vector_width;
//               original_index2 = hash[1] % bit_vector_width;
               bit_index1      = original_index1 % BITS_PER_CHAR;
//               bit_index2      = original_index2 % BITS_PER_CHAR;

               #pragma omp atomic
               bit_vector[original_index1 / BITS_PER_CHAR] |= BIT_MASK[bit_index1];

//               #pragma omp atomic
//               bit_vector[original_index2 / BITS_PER_CHAR] |= BIT_MASK[bit_index2];
            }
         }
      }
   }

   kmer_vector.clear();
    

   //--------------------------------------------------
   // do bit-wise or for the bloom filter data
   //--------------------------------------------------
   std::size_t receive_buffer_size_byte;
   std::size_t num_iterations;
   std::size_t buffer_size_residue;

   // determine the size of the receive buffer
   if (bit_vector_width_byte >= BLOOM_RCV_BUFFER_SIZE) {
      receive_buffer_size_byte = BLOOM_RCV_BUFFER_SIZE;
   }
   else {
      receive_buffer_size_byte = bit_vector_width_byte;
   }

   num_iterations      = bit_vector_width_byte / receive_buffer_size_byte;
   buffer_size_residue = bit_vector_width_byte % receive_buffer_size_byte;


   unsigned char* receive_buffer(new unsigned char[static_cast<std::size_t>(receive_buffer_size_byte)]);

   //MPI_Barrier(comm_node);

   // bit-wise or
   for (std::size_t it_iter = 0; it_iter < num_iterations; it_iter++) {
      MPI_Allreduce(&bit_vector[it_iter * receive_buffer_size_byte], receive_buffer, receive_buffer_size_byte, MPI_UNSIGNED_CHAR, MPI_BOR, group_comm);
      memcpy(&bit_vector[it_iter * receive_buffer_size_byte], receive_buffer, receive_buffer_size_byte);
   }

   MPI_Allreduce(&bit_vector[num_iterations * receive_buffer_size_byte], receive_buffer, buffer_size_residue, MPI_UNSIGNED_CHAR, MPI_BOR, group_comm);
   memcpy(&bit_vector[num_iterations * receive_buffer_size_byte], receive_buffer, buffer_size_residue);

//   boost::filesystem::path path_kmc_pre(kmc_prefix + ".kmc_pre");
//   boost::filesystem::path path_kmc_suf(kmc_prefix + ".kmc_suf");
//   boost::filesystem::remove(path_kmc_pre);
//   boost::filesystem::remove(path_kmc_suf);
   std::string path_kmc_pre = kmc_prefix + ".kmc_pre";
   std::string path_kmc_suf = kmc_prefix  + ".kmc_suf";

   remove(path_kmc_pre.c_str());
   remove(path_kmc_suf.c_str());
   //--------------------------------------------------
   // dump the bloom filter data
   //--------------------------------------------------
   if (group_rank == 0) {
      // open the bloom filter dump file
      std::ofstream f_bf_dump_size;
      std::ofstream f_bf_dump_data;
      f_bf_dump_size.open(c_inst_args.bf_size_file_name.c_str());
      f_bf_dump_data.open(c_inst_args.bf_data_file_name.c_str(), std::ios::binary);

      // check the bloom filter dump file
      if (f_bf_dump_size.is_open() == false) {
         std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.bf_size_file_name << std::endl << std::endl;
         MPI_Abort(MPI_COMM_WORLD, 201);
      }

      if (f_bf_dump_data.is_open() == false) {
         std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.bf_data_file_name << std::endl << std::endl;
         MPI_Abort(MPI_COMM_WORLD, 202);
      }

      // write the size file
      // order: <bit vector size in bits> <bit vector size in bytes> <# of hash functions> <k-mer occurrence threshold> <hash function random seed> <# of unique solid k-mers>
      f_bf_dump_size << bit_vector_width << " ";
      f_bf_dump_size << bit_vector_width_byte << " ";
      f_bf_dump_size << num_hash_func << " ";
      f_bf_dump_size << kmer_length << " ";
      f_bf_dump_size << c_inst_args.kmer_occurrence_threshold << " ";
      f_bf_dump_size << c_inst_args.random_seed << " ";
      f_bf_dump_size << num_unique_solid_kmers << std::endl;

      // dump the bit vector
      f_bf_dump_data.write((const char*)&bit_vector[0], bit_vector_width_byte);

      f_bf_dump_data.close();
      f_bf_dump_size.close();
   }
   free(hash_seed);
   free(bit_vector);
   // purge allocated memory
   //bit_vector.clear();

   delete[] receive_buffer;

   time(&rawtime);
   c_inst_time.end_program_kmers_into_bloom_filter = asctime(localtime(&rawtime));

   if (rank_node == 0) {
      std::cout << "     Programming k-mers to the Bloom filter: done" << std::endl;
      std::cout << std::endl;
   }
}


//----------------------------------------------------------------------
// read_bloom_filter_parameters
//----------------------------------------------------------------------
void C_generate_bloom::read_bloom_filter_parameters(const C_arg& c_inst_args, C_time& c_inst_time) {
   // open the bloom filter dump file
   std::ifstream f_bf_dump_size;
   f_bf_dump_size.open(c_inst_args.loaded_bf_size_file_name.c_str());

   // check the bloom filter dump file
   if (f_bf_dump_size.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.loaded_bf_size_file_name << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 203);
   }

   // read the size file
   // order: <bit vector size> <# of hash functions> <k-mer occurrence threshold> <# of unique solid k-mers>
   std::size_t kmer_length_read;

   f_bf_dump_size
                  >> bit_vector_width
                  >> bit_vector_width_byte
                  >> num_hash_func
                  >> kmer_length_read
                  >> kmer_occurrence_threshold
                  >> random_seed
                  >> num_unique_solid_kmers;
    num_hash_func_real = num_hash_func / 2;
   if (f_bf_dump_size.good() == false) {
      std::cout << std::endl << "ERROR: Wrong number of items in " << c_inst_args.loaded_bf_size_file_name << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 204);
   }

   // compare the k-mer length between the loaded data and the command line input
   if (kmer_length_read != kmer_length) {
      std::cout << std::endl << "ERROR: k-mer length in the loaded data: " << kmer_length << " vs in the command line: " << kmer_length << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 205);
   }

   // set rank_node_text
   // this will be needed when errors are corrected
   std::stringstream rank_node_stream;
   rank_node_stream << std::setw(5) << std::setfill('0') << rank_node;
   rank_node_text = rank_node_stream.str();

   f_bf_dump_size.close();
}


