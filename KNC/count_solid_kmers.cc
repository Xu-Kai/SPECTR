#include "count_solid_kmers.h"



//----------------------------------------------------------------------
// count_kmers
//----------------------------------------------------------------------
void C_count_solid_kmers::count_kmers(const C_arg& c_inst_args, C_time& c_inst_time, int group_rank, MPI_Comm& group_comm, int node_num) {
   time_t rawtime;
   time(&rawtime);
   c_inst_time.start_count_kmers = asctime(localtime(&rawtime));

   if (rank_node == 0) {
      std::cout << "Counting the number of k-mers" << std::endl;
   }
    /*
   std::cout << num_unique_solid_kmers << 
      " " << kmer_occurrence_threshold << " " 
     << valey_point << " "
     << kmer_length <<" " 
     << rank_node << " "
    << size_node << " " 
    <<get_valey_point << std::endl; 
*/
   // set rank_node_text
   std::stringstream rank_node_stream;
   rank_node_stream << std::setw(5) << std::setfill('0') << rank_node;
   rank_node_text = rank_node_stream.str();

   kmc_prefix = c_inst_args.kmc_prefix + "." + rank_node_text;
   //std::cout << kmc_prefix << std::endl;
    
   // run kmc
   run_kmc(c_inst_args, group_rank, group_comm, node_num);
    
   // the k-mer occurrence threshold should be determined
   if (c_inst_args.set_kmer_occurrence_threshold == false) {
      // determine the threshold
      // the k-mer occurrence histogram is also generated
      determine_kmer_occurrence_threshold(c_inst_args, group_rank, group_comm, node_num);
   }
   // the k-mer occurrence threshold is pre-determined
   else {
      generate_kmer_occurrence_histogram(c_inst_args, group_rank, group_comm, node_num);
   }
    
   if (group_rank == 0) {
      // count the number of unique solid k-mers
      count_unique_solid_kmers(c_inst_args, group_rank, group_comm, node_num);

      // write the k-mer occurrence histogram
      write_kmer_histogram(c_inst_args, group_rank, group_comm, node_num);

      std::cout << "     Number of unique solid k-mers: " << num_unique_solid_kmers << std::endl;
      std::cout << "     Counting the number of k-mers: done" << std::endl;
      std::cout << std::endl;
   }
   delete[] num_occurrences_histogram;
   delete[] reduced_num_occurrences_histogram;
    
   time(&rawtime);
   c_inst_time.end_count_kmers = asctime(localtime(&rawtime));
}



//----------------------------------------------------------------------
// run_kmc
//----------------------------------------------------------------------
void C_count_solid_kmers::run_kmc(const C_arg& c_inst_args, int group_rank, MPI_Comm& group_comm, int node_num) {
   // command for running kmc
   std::string cmd(c_inst_args.kmc_binary);

   // stream buffer for numbers
   std::stringstream sstream_tmp;

   // add k
   sstream_tmp << c_inst_args.kmer_length;
   cmd += (" -k" + sstream_tmp.str());

   // add max memory usage
   sstream_tmp.str("");
   sstream_tmp << c_inst_args.max_mem;
   cmd += (" -m" + sstream_tmp.str() + " -sm -fq -v");

   size_node = 1;
   // add the number of nodes
   sstream_tmp.str("");
   sstream_tmp << size_node;
   cmd += (" -d" + sstream_tmp.str());

   // add the rank of a current node
   sstream_tmp.str("");
   sstream_tmp << rank_node;
   cmd += (" -a" + sstream_tmp.str());

   // add the number of threads
   // if specified
   if (c_inst_args.smpthread > 0) {
      sstream_tmp.str("");
      sstream_tmp << c_inst_args.smpthread;
      cmd += (" -t" + sstream_tmp.str());
   }
   else {
      std::cout << std::endl << "ERROR: The number of threads is 0" << std::endl << std::endl;
//      MPI_Abort(MPI_COMM_WORLD, 300);
   }

   // the k-mer occurrence threshold is already given
   if (c_inst_args.set_kmer_occurrence_threshold == true) {
      sstream_tmp.str("");
      sstream_tmp << c_inst_args.kmer_occurrence_threshold;

      cmd += (" -ci" + sstream_tmp.str());
   }

   // two input read files
   if (c_inst_args.paired_read == true && node_num == 1) {
      std::ofstream f_input_file_list;
      f_input_file_list.open(c_inst_args.input_list_file.c_str());
      //std::cout << "file input list: " << c_inst_args.input_list_file << std::endl;

      if (f_input_file_list.is_open() == false) {
         std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.input_list_file << std::endl << std::endl;
         MPI_Abort(MPI_COMM_WORLD, 1200);
      }
      else {
         f_input_file_list << c_inst_args.read_file_name1 << std::endl;
         f_input_file_list << c_inst_args.read_file_name2 << std::endl;
      }

      f_input_file_list.close();

      cmd += (" @" + c_inst_args.input_list_file);
   }
   // one input read file
   else if(c_inst_args.paired_read == true){
      if(group_rank == 0)
        cmd += (" " + c_inst_args.read_file_name1);
      else if(group_rank == 1) 
        cmd += (" " + c_inst_args.read_file_name2);
   }else{
        cmd += (" " + c_inst_args.read_file_name);
    }


   cmd += (" " + kmc_prefix);
   cmd += (" " + c_inst_args.prefix + "." + rank_node_text);
   cmd += ("> " + kmc_prefix + ".log 2>&1");


   //std::cout << cmd << std::endl;
   //std::cout << "rank_node_text" << c_inst_args.prefix << std::endl;
   //std::cout << "rank_node_text" << rank_node_text << std::endl;
   std::string kmc_tmp = c_inst_args.prefix + "."+rank_node_text ;
   char buf[360];
   getcwd(buf, sizeof(buf));
   std::string kmc_tmp_absolute(buf);
   kmc_tmp_absolute = kmc_tmp_absolute + "/" +  kmc_tmp;
   //std::cout << kmc_tmp_absolute<<std::endl;
   //std::cout << cmd <<std::endl;
   mkdir(kmc_tmp_absolute.c_str(), 0700);
//   // run kmc
   if (system(cmd.c_str()) != 0) {
      std::cout << std::endl << "ERROR: KMC is abnormally terminated. See " << c_inst_args.prefix << "*.log" << std::endl << std::endl;
//      MPI_Abort(MPI_COMM_WORLD, 301);
   }
//   // remove the kmc temporary directory
   rmdir(kmc_tmp_absolute.c_str());
}



//----------------------------------------------------------------------
// determine_kmer_occurrence_threshold
//----------------------------------------------------------------------
void C_count_solid_kmers::determine_kmer_occurrence_threshold(const C_arg& c_inst_args, int group_rank, MPI_Comm& group_comm, int node_num) {
   // generate the k-mer occurrence histogram
   generate_kmer_occurrence_histogram(c_inst_args, group_rank, group_comm, node_num);

   if (group_rank == 0) {
      // determine the k-mer occurrence threshold
      bloom_type prev_value(0);

      for (std::size_t it_histo = KMC_DEFAULT_MIN_COUNT - 1; it_histo < KMER_HISTOGRAM_SIZE - 1; it_histo++) {
         // find the valey point
         if ((reduced_num_occurrences_histogram[it_histo] > prev_value) && (it_histo > 1) && (get_valey_point == false)) {
            valey_point = it_histo - 1;
            get_valey_point = true;
         }

         prev_value = reduced_num_occurrences_histogram[it_histo];
      }

      if (get_valey_point == true) {
         if (valey_point <= MAX_KMER_THRESHOLD) {
            std::cout << "     k-mer occurrence threshold   : " << valey_point << std::endl;

            kmer_occurrence_threshold = valey_point;
         }
         else {
            std::cout << "WARNING: The automatically determined k-mer threshold " << valey_point << " is larger than the maximum value " << MAX_KMER_THRESHOLD << std::endl;
            std::cout << "         The k-mer threshold is set to " << MAX_KMER_THRESHOLD << std::endl;

            kmer_occurrence_threshold = MAX_KMER_THRESHOLD;
         }
      }
      else {
         std::cout << std::endl << "ERROR: No valey point exists in the histogram and no k-mer occurrence threshold is given. Please, try different k values." << valey_point << std::endl << std::endl;
//         MPI_Abort(MPI_COMM_WORLD, 302);
      }
   }
}



//----------------------------------------------------------------------
// generate_kmer_occurrence_histogram
//----------------------------------------------------------------------
void C_count_solid_kmers::generate_kmer_occurrence_histogram(const C_arg& c_inst_args, int group_rank, MPI_Comm& group_comm, int node_num) {
   // initialize the histgram of k-mer occurrences
   num_occurrences_histogram         = new bloom_type_count[KMER_HISTOGRAM_SIZE];
   reduced_num_occurrences_histogram = new bloom_type_count[KMER_HISTOGRAM_SIZE];

   for (std::size_t it_histo = 0; it_histo < KMER_HISTOGRAM_SIZE; it_histo++) {
      num_occurrences_histogram[it_histo] = 0;
   }

   // open the k-mer database
   if (!kmer_database.OpenForListing(kmc_prefix)) {
      std::cout << std::endl << "ERROR: Cannot open " << kmc_prefix << std::endl << std::endl;
//      MPI_Abort(MPI_COMM_WORLD, 303);
   }

   CKmerAPI kmer((uint32)c_inst_args.kmer_length);

   float num_occurrences;

   // iterate the database
   // index: occurrences
   while (kmer_database.ReadNextKmer(kmer, num_occurrences)) {
      // update the histogram
      if (num_occurrences > (KMER_HISTOGRAM_SIZE - 1)) {
         num_occurrences_histogram[KMER_HISTOGRAM_SIZE - 1]++;
      }
      else {
         num_occurrences_histogram[(std::size_t)num_occurrences]++;
      }
   }

   //MPI_Barrier(comm_node);
//   reduced_num_occurrences_histogram = 0;
   // reduce the k-mer histogram
   MPI_Reduce(num_occurrences_histogram, reduced_num_occurrences_histogram, KMER_HISTOGRAM_SIZE, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, group_comm);
  //for(int i = 0; i < KMER_HISTOGRAM_SIZE; i++){
//	   reduced_num_occurrences_histogram[i] = num_occurrences_histogram[i];
 // }

   // fill num_occurrences_histogram[1]
   // KMC_DEFAULT_MIN_COUNT: 2
   reduced_num_occurrences_histogram[KMC_DEFAULT_MIN_COUNT - 1] = reduced_num_occurrences_histogram[KMC_DEFAULT_MIN_COUNT] + 1000;
}



//----------------------------------------------------------------------
// count_unique_solid_kmers
//----------------------------------------------------------------------
void C_count_solid_kmers::count_unique_solid_kmers(const C_arg& c_inst_args, int group_rank, MPI_Comm& group_comm, int node_num) {
   num_unique_solid_kmers = 0;

   for (std::size_t it_histo = kmer_occurrence_threshold; it_histo < KMER_HISTOGRAM_SIZE; it_histo++) {
      num_unique_solid_kmers += reduced_num_occurrences_histogram[it_histo];
   }
}



//----------------------------------------------------------------------
// write_kmer_histogram
//----------------------------------------------------------------------
void C_count_solid_kmers::write_kmer_histogram(const C_arg& c_inst_args, int group_rank, MPI_Comm& group_comm, int node_num) {
   std::ofstream f_histo;
   f_histo.open(c_inst_args.kmer_histo_file_name.c_str());

   // set the minimum histogram index
   std::size_t min_histogram_index;

   // the k-mer occurrence threshold is not pre-determined
   if (c_inst_args.set_kmer_occurrence_threshold == false) {
      min_histogram_index = KMC_DEFAULT_MIN_COUNT;
   }
   // the k-mer occurrence threshold is pre-determined
   else {
      min_histogram_index = kmer_occurrence_threshold;
   }

   if (f_histo.is_open()) {
      for (std::size_t it_histo = min_histogram_index; it_histo < (KMER_HISTOGRAM_SIZE - 1); it_histo++) {
         f_histo << std::setw(7) << it_histo << ": " << std::setw(10) << reduced_num_occurrences_histogram[it_histo] << std::endl;
      }

      f_histo << "<=" << std::setw(5) << KMER_HISTOGRAM_SIZE - 1 << ": " << std::setw(10) << reduced_num_occurrences_histogram[KMER_HISTOGRAM_SIZE - 1] << std::endl;
   }
   else {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.kmer_histo_file_name << std::endl << std::endl;
//      MPI_Abort(MPI_COMM_WORLD, 304);
   }

   f_histo.close();
}
