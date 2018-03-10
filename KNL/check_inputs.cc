#include "check_inputs.h"







//----------------------------------------------------------------------
// check_read_file
//----------------------------------------------------------------------
void C_check_read::check_read_file(const C_arg& c_inst_args, C_time& c_inst_time, int group_rank, MPI_Comm &group_comm, int node_num) {
   // start measuring run-time
   time_t rawtime;
   time(&rawtime);
   c_inst_time.start_check_read_file = asctime(localtime(&rawtime));

   int my_rank = 0;
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
   if(my_rank == 0)
       std::cout << "Checking input read files" << std::endl;

   // paired
   if (c_inst_args.paired_read == true) {
      check_read_file_fastq_paired(c_inst_args, group_rank, group_comm, node_num);
   }
   // single
   else {
      check_read_file_fastq_single(c_inst_args);
   }

   time(&rawtime);
   c_inst_time.end_check_read_file = asctime(localtime(&rawtime));
}



//----------------------------------------------------------------------
// check_read_file_fastq_single
//----------------------------------------------------------------------
void C_check_read::check_read_file_fastq_single(const C_arg& c_inst_args) {
   //----------------------------------------------------------------------
   // construct a quality score histogram
   //----------------------------------------------------------------------
   //qs_histo_vec.resize(QS_HISTOGRAM_MAX + 1, 0);
    qs_histo_vec = (uint64_t*)malloc((QS_HISTOGRAM_MAX + 1) * sizeof(uint64_t));
    memset(qs_histo_vec, 0, (QS_HISTOGRAM_MAX+1) * sizeof(uint64_t));

   if (c_inst_args.gzipped_input_read) {
      construct_quality_score_histogram_single_gzipped(c_inst_args);
   }
   else {
      construct_quality_score_histogram_single_unzipped(c_inst_args);
   }

   //----------------------------------------------------------------
   // find the max/min quality scores
   //----------------------------------------------------------------
   max_quality_score = -1000;
   min_quality_score =  1000;

   // max
   for (int it_histo = 0; it_histo <= QS_HISTOGRAM_MAX; it_histo++) {
      if (qs_histo_vec[it_histo] > 0) {
         max_quality_score = it_histo;
      }
   }

   // min
   for (int it_histo = QS_HISTOGRAM_MAX; it_histo >= 0; it_histo--) {
      if (qs_histo_vec[it_histo] > 0) {
         min_quality_score = it_histo;
      }
   }

   //----------------------------------------------------------------
   // determine the quality score offset
   //----------------------------------------------------------------
   // http://en.wikipedia.org/wiki/FASTQ_format
   if (min_quality_score < 33) {
      std::cout << std::endl << "ERROR: There is a quality score < 33" << std::endl << std::endl;
//      MPI_Abort(MPI_COMM_WORLD, 100);
   }
   // 33 <= min_quality_score <= 58
   else if (min_quality_score <= 58) {
      if (max_quality_score <= 74) {
         quality_score_offset = PHRED33;
      }
      else {
         std::cout << std::endl << "ERROR: Irregular quality score range " << min_quality_score << "-" << max_quality_score << std::endl << std::endl;
//         MPI_Abort(MPI_COMM_WORLD, 101);
      }
   }
   // 59 <= min_quality_score <= 74
   else if (min_quality_score <= 74) {
      if (max_quality_score <= 74) {
         std::cout << std::endl << "WARNING: Hard to determine the quality score offset (min: " << min_quality_score << ", max: " << max_quality_score << ")" << std::endl;
         std::cout <<              "         Phred+33 will be used" << std::endl << std::endl;
         quality_score_offset = PHRED33;
      }
      // max_quality_score >= 75
      else {
         quality_score_offset = PHRED64;
      }
   }
   // min_quality_score >= 75
   else {
      quality_score_offset = PHRED64;
   }

   //----------------------------------------------------------------------
   // write a histogram file
   //----------------------------------------------------------------------
   std::size_t total_bases(0);

   std::ofstream f_histo;

   f_histo.open(c_inst_args.qs_histo_file_name.c_str());

   if (f_histo.is_open()) {
      for (std::size_t it_histo = quality_score_offset; it_histo <= quality_score_offset + MAX_PHRED; it_histo++) {
         f_histo << std::setw(4) << it_histo - quality_score_offset << ": " << std::setw(14) << qs_histo_vec[it_histo] << std::endl;

         total_bases += qs_histo_vec[it_histo];
      }
   }

   f_histo.close();

   // find the threshold
   std::size_t partial_sum(0);

   bool set_quality_score_cutoff(false);
   bool set_extremely_low_quality_score(false);

   for (std::size_t it_histo = 0; it_histo <= QS_HISTOGRAM_MAX; it_histo++) {
      partial_sum += qs_histo_vec[it_histo];

      if (((1.0 * partial_sum / total_bases) >= QS_CUTOFF_RATIO) && (set_quality_score_cutoff == false)) {
         quality_score_cutoff = it_histo - quality_score_offset;
         set_quality_score_cutoff = true;
      }

      if (((1.0 * partial_sum / total_bases) >= QS_EXTREMELY_LOW_RATIO) && (set_extremely_low_quality_score == false)) {
         extremely_low_quality_score = it_histo - quality_score_offset;
         set_extremely_low_quality_score = true;

         if ((1.0 * partial_sum / total_bases) >= (QS_EXTREMELY_LOW_RATIO + 0.05)) {
            std::cout << std::endl;
            std::cout << "WARNING: Some quality score thresholds are set to a high value" << std::endl;
            std::cout << "         because overall quality scores are too bad." << std::endl;
            std::cout << "         It may cause long runtime and large memory usage." << std::endl << std::endl;
         }
      }
   }
    
   int my_rank = 0;
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
   if(my_rank == 0){
   std::cout << "     Number of reads           : " << num_reads << std::endl;
   std::cout << "     Quality score offset      : " << quality_score_offset << std::endl;
   std::cout << "     Quality score threshold   : " << quality_score_cutoff << std::endl;
   std::cout << "     Low quality score threhold: " << extremely_low_quality_score << std::endl;
   std::cout << "     Checking input read files : done" << std::endl;
   std::cout << std::endl;
   }
}



//----------------------------------------------------------------------
// check_read_file_fastq_paired
//----------------------------------------------------------------------
void C_check_read::check_read_file_fastq_paired(const C_arg& c_inst_args, int group_rank, MPI_Comm& group_comm, int node_num) {
   //----------------------------------------------------------------------
   // construct a quality score histogram
   //----------------------------------------------------------------------
   //qs_histo_vec.resize(QS_HISTOGRAM_MAX + 1, 0);
    qs_histo_vec = (uint64_t*)malloc((QS_HISTOGRAM_MAX + 1) * sizeof(uint64_t));
    memset(qs_histo_vec, 0, (QS_HISTOGRAM_MAX+1) * sizeof(uint64_t));

   if (c_inst_args.gzipped_input_read) {
      construct_quality_score_histogram_paired_gzipped(c_inst_args);
   }
   else {
      construct_quality_score_histogram_paired_unzipped(c_inst_args, group_rank, group_comm, node_num);
   }

   //----------------------------------------------------------------
   // find the max/min quality scores
   //----------------------------------------------------------------
   max_quality_score = -1000;
   min_quality_score =  1000;
if(group_rank == 0){
    std::cout << "determine quality_offset" << std::endl;
   // max
   for (int it_histo = 0; it_histo <= QS_HISTOGRAM_MAX; it_histo++) {
    //std::cout << it_histo <<" " <<qs_histo_vec[it_histo] << std::endl;
      if (qs_histo_vec[it_histo] > 0) {
         max_quality_score = it_histo;
      }
   }

   // min
   for (int it_histo = QS_HISTOGRAM_MAX; it_histo >= 0; it_histo--) {
      if (qs_histo_vec[it_histo] > 0) {
         min_quality_score = it_histo;
      }
   }

   //----------------------------------------------------------------
   // determine the quality score offset
   //----------------------------------------------------------------
   // http://en.wikipedia.org/wiki/FASTQ_format
   if (min_quality_score < 33) {
      std::cout << std::endl << "ERROR: There is a quality score < 33" << std::endl << std::endl;
//      MPI_Abort(MPI_COMM_WORLD, 102);
   }
   // 33 <= min_quality_score <= 58
   
   else if (min_quality_score <= 58) {
      if (max_quality_score <= 78) {
         quality_score_offset = PHRED33;
      }
      else {
         std::cout << std::endl << "ERROR: Irregular quality score range " << min_quality_score << "-" << max_quality_score << std::endl << std::endl;
//         MPI_Abort(MPI_COMM_WORLD, 103);
      }
   }
   
   // 59 <= min_quality_score <= 74
   else if (min_quality_score <= 75) {
      if (max_quality_score <= 75) {
         std::cout << std::endl << "WARNING: Hard to determine the quality score offset (min: " << min_quality_score << ", max: " << max_quality_score << ")" << std::endl;
         std::cout <<              "         Phred+33 will be used" << std::endl << std::endl;
         quality_score_offset = PHRED33;
      }
      // max_quality_score >= 75
      else {
         quality_score_offset = PHRED64;
      }
   }
   // min_quality_score >= 75
   else {
      quality_score_offset = PHRED64;
   }

   //----------------------------------------------------------------------
   // write a histogram file
   //----------------------------------------------------------------------
   std::size_t total_bases(0);

   std::ofstream f_histo;

   f_histo.open(c_inst_args.qs_histo_file_name.c_str());

   if (f_histo.is_open()) {
      for (std::size_t it_histo = quality_score_offset; it_histo <= quality_score_offset + MAX_PHRED; it_histo++) {
         f_histo << std::setw(4) << it_histo - quality_score_offset << ": " << std::setw(14) << qs_histo_vec[it_histo] << std::endl;

         total_bases += qs_histo_vec[it_histo];
      }
   }

   f_histo.close();

   // find the threshold
   std::size_t partial_sum(0);

   bool set_quality_score_cutoff(false);
   bool set_extremely_low_quality_score(false);

   for (std::size_t it_histo = 0; it_histo <= QS_HISTOGRAM_MAX; it_histo++) {
      partial_sum += qs_histo_vec[it_histo];

      if (((1.0 * partial_sum / total_bases) >= QS_CUTOFF_RATIO) && (set_quality_score_cutoff == false)) {
         quality_score_cutoff = it_histo - quality_score_offset;
         set_quality_score_cutoff = true;
      }

      if (((1.0 * partial_sum / total_bases) >= QS_EXTREMELY_LOW_RATIO) && (set_extremely_low_quality_score == false)) {
         extremely_low_quality_score = it_histo - quality_score_offset;
         set_extremely_low_quality_score = true;

         if ((1.0 * partial_sum / total_bases) >= (QS_EXTREMELY_LOW_RATIO + 0.05)) {
            std::cout << std::endl;
            std::cout << "WARNING: Some quality score thresholds are set to a high value" << std::endl;
            std::cout << "         because overall quality scores are too bad." << std::endl;
            std::cout << "         It may cause long runtime and large memory usage." << std::endl << std::endl;
         }
      }
   }

   std::cout << "     Number of pairs           : " << num_reads << std::endl;
   std::cout << "     Quality score offset      : " << quality_score_offset << std::endl;
   std::cout << "     Quality score threshold   : " << quality_score_cutoff << std::endl;
   std::cout << "     Low quality score threhold: " << extremely_low_quality_score << std::endl;
   std::cout << "     Checking input read files : done" << std::endl << std::endl;
   }
    
   free(qs_histo_vec);
}



//----------------------------------------------------------------------
// construct_quality_score_histogram_single_unzipped
//----------------------------------------------------------------------
void C_check_read::construct_quality_score_histogram_single_unzipped(const C_arg& c_inst_args) {

   // open an input read file
//   boost::iostreams::mapped_file_source f_read;
	   // find out the size of read files
//	   f_read.open(c_inst_args.read_file_name.c_str());

//	   if (!f_read.is_open()) {
//	      std::cout << std::endl << "ERROR: Cannot open " << f_read << std::endl << std::endl;
//	      MPI_Abort(MPI_COMM_WORLD, 104);
//	   }

//	   read_file_size_byte      = f_read.size();
//	   read_file_unit_size_byte = f_read.alignment();
		FILE * f_read1;

	  f_read1 = fopen (c_inst_args.read_file_name.c_str(),"rb");
	  if (f_read1==NULL) perror ("Error opening file");
	  else
	  {
	    fseek (f_read1, 0, SEEK_END);   // non-portable
	    read_file_size_byte=ftell (f_read1);
	    read_file_unit_size_byte = 100*1024*1024;
	    fclose (f_read1);
	  }







   std::size_t num_units(ceil(1.0 * read_file_size_byte / read_file_unit_size_byte));

   // check whether the size of the input file is 0
   if (read_file_size_byte == 0) {
      std::cout << std::endl << "ERROR: The size of " << c_inst_args.read_file_name << " is 0" << std::endl << std::endl;
//      MPI_Abort(MPI_COMM_WORLD, 105);
   }

   // check whether the alignment unit offset is larger than MMAP_FILE_SIZE
   if (MMAP_FILE_SIZE <= read_file_unit_size_byte) {
      std::cout << std::endl << "ERROR: The block size for "  << " should be larger than "  << std::endl << std::endl;
//      MPI_Abort(MPI_COMM_WORLD, 106);
   }

//   // check whether the size of the input file < read_file_unit_size_byte X <number of nodes>
//   // if it is true, some cores may not have one alignment unit in the input file
//   if (num_units < (std::size_t)size_node) {
//      std::cout << std::endl << "ERROR: The size of " << c_inst_args.read_file_name << " is too small. Run BLESS again without MPI." << std::endl << std::endl;
////      MPI_Abort(MPI_COMM_WORLD, 107);
//   }


   FILE* pfile = fopen (c_inst_args.read_file_name.c_str(),"rb");
	int buffer_size = 100 * 1024 * 1024; //100M
	char* buffer = (char*) malloc(buffer_size);
	uint32_t read_size = 1;
	uint32_t flag = 0;
	uint32_t iter = 0;
	while (read_size != 0) {
		read_size = fread(buffer, 1, buffer_size, pfile);
		int pointer = 0;
		for (iter = 0; iter < read_size; iter++) {
			if (flag == 0) {
				if (buffer[iter] == '\n') {
					flag = 1;
					pointer = 0;
				}

			} else if (flag == 1) {

				if (buffer[iter] == '\n') {
					pointer = 0;
					flag = 2;
				}
			} else if (flag == 2) {
				if (buffer[iter] == '\n') {
					pointer = 0;
					flag = 3;
				}
			} else if (flag == 3) {
				if (buffer[iter] == '\n') {
					flag = 0;
					pointer = 0;
					num_reads ++;
				}else
					qs_histo_vec[buffer[iter]] ++;
			}
		}
	}

	fclose(pfile);
}



//----------------------------------------------------------------------
// construct_quality_score_histogram_paired_unzipped
//----------------------------------------------------------------------
void C_check_read::construct_quality_score_histogram_paired_unzipped(const C_arg& c_inst_args, int group_rank, MPI_Comm& group_comm, int node_num) {


		  FILE * f_read_1;

		  f_read_1 = fopen (c_inst_args.read_file_name1.c_str(),"rb");
		  if (f_read_1==NULL) perror ("Error opening file");
		  else
		  {
		    fseek (f_read_1, 0, SEEK_END);   // non-portable
		    read_file_size_byte1=ftell (f_read_1);
		    read_file_unit_size_byte1 = 100*1024*1024;
		    fclose (f_read_1);
		  }
			FILE * f_read_2;

		  f_read_2 = fopen (c_inst_args.read_file_name2.c_str(),"rb");
		  if (f_read_2==NULL) perror ("Error opening file");
		  else
		  {
		    fseek (f_read_2, 0, SEEK_END);   // non-portable
		    read_file_size_byte2=ftell (f_read_2);
		    read_file_unit_size_byte2 = 100*1024*1024;
		    fclose (f_read_2);
		  }


		   std::size_t num_units1(ceil(1.0 * read_file_size_byte1 / read_file_unit_size_byte1));
		   std::size_t num_units2(ceil(1.0 * read_file_size_byte2 / read_file_unit_size_byte2));





	   std::size_t num_units(ceil(1.0 * read_file_size_byte / read_file_unit_size_byte));

	   // check whether the size of the input file is 0
	   if (read_file_size_byte1 == 0) {
	      std::cout << std::endl << "ERROR: The size of " << c_inst_args.read_file_name << " is 0" << std::endl << std::endl;
	//      MPI_Abort(MPI_COMM_WORLD, 105);
	   }
	   if (read_file_size_byte2 == 0) {
	      std::cout << std::endl << "ERROR: The size of " << c_inst_args.read_file_name << " is 0" << std::endl << std::endl;
	//      MPI_Abort(MPI_COMM_WORLD, 105);
	   }
	   // check whether the alignment unit offset is larger than MMAP_FILE_SIZE
	   if (MMAP_FILE_SIZE <= read_file_unit_size_byte1) {
	      std::cout << std::endl << "ERROR: The block size for "  << " should be larger than "  << std::endl << std::endl;
	//      MPI_Abort(MPI_COMM_WORLD, 106);
	   }
	   // check whether the alignment unit offset is larger than MMAP_FILE_SIZE
	   if (MMAP_FILE_SIZE <= read_file_unit_size_byte2) {
	      std::cout << std::endl << "ERROR: The block size for "  << " should be larger than "  << std::endl << std::endl;
	//      MPI_Abort(MPI_COMM_WORLD, 106);
	   }
//	   // check whether the size of the input file < read_file_unit_size_byte X <number of nodes>
//	   // if it is true, some cores may not have one alignment unit in the input file
//	   if (num_units < (std::size_t)size_node) {
//	      std::cout << std::endl << "ERROR: The size of " << c_inst_args.read_file_name << " is too small. Run BLESS again without MPI." << std::endl << std::endl;
//	//      MPI_Abort(MPI_COMM_WORLD, 107);
//	   }
	   uint64_t num_reads_1 = (0);
	   uint64_t num_reads_2 = (0);
        if(group_rank == 0 || node_num == 1){
	   FILE* pfile = fopen (c_inst_args.read_file_name1.c_str(),"r");
        char line[512];
        uint64_t line_count = 0;
        while(fgets(line, sizeof(line), pfile)){
            line_count += 1;
            if(line_count % 4 == 0){
                int line_c = 0;
                while(line[line_c] != '\n'){ 
					qs_histo_vec[line[line_c]] ++;
                    line_c ++;
                }
                num_reads_1 += 1;
            }


        } 
		fclose(pfile);
        }

        if(group_rank == 1 || node_num == 1){
	   FILE* pfile = fopen (c_inst_args.read_file_name2.c_str(),"r");
        char line[512];
        uint64_t line_count = 0;
        while(fgets(line, sizeof(line), pfile)){
            line_count += 1;
            if(line_count % 4 == 0){
                int line_c = 0;
                while(line[line_c] != '\n'){ 
					qs_histo_vec[line[line_c]] ++;
                    line_c ++;
                }
                num_reads_2 += 1;
            }
        } 
		fclose(pfile);
        }
        if(node_num > 1)
            MPI_Bcast(&num_reads_2, 1, MPI_LONG_LONG_INT, 1, group_comm);
		// compair the number of reads in two input files
        if(group_rank == 0){
		   if (num_reads_1 == num_reads_2) {
		      num_reads = num_reads_1;
		   }
		   else {
		      std::cout << std::endl << "ERROR: Number of lines in two input files are not same" << std::endl << std::endl;
		    }
        }
    uint64_t *reduced_histo_vec = (uint64_t*)malloc((QS_HISTOGRAM_MAX + 1) * sizeof(uint64_t));
    memset(reduced_histo_vec, 0, (QS_HISTOGRAM_MAX+1) * sizeof(uint64_t));
    MPI_Reduce(qs_histo_vec, reduced_histo_vec, QS_HISTOGRAM_MAX + 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, group_comm);   

    for(int i = 0; i < QS_HISTOGRAM_MAX + 1; i ++)
      qs_histo_vec[i] = reduced_histo_vec[i]; 
    free(reduced_histo_vec);




}



//----------------------------------------------------------------------
// construct_quality_score_histogram_single_gzipped
//----------------------------------------------------------------------
void C_check_read::construct_quality_score_histogram_single_gzipped(const C_arg& c_inst_args) {
}



//----------------------------------------------------------------------
// construct_quality_score_histogram_paired_gzipped
//----------------------------------------------------------------------
void C_check_read::construct_quality_score_histogram_paired_gzipped(const C_arg& c_inst_args) {
}
