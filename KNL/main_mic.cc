//#include<stdio.h>
//#include<stdlib.h>
//#include<stdint.h>
#include "parse_args.h"
#include "define.h"
#include "time.h"
#include "mic_correct_errors.h"
//#include "cpu_correct_errors.h"
#include "check_inputs.h"
#include "count_solid_kmers.h"
//#include "count_solid_kmers2.h"
#include "generate_bloom_filter.h"
#include "data_manager.h"
#include <mpi.h>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <semaphore.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/time.h>

#define GETTIME(x) gettimeofday(&x, NULL)
#define DIFFTIME(x,y) y.tv_sec-x.tv_sec + (y.tv_usec - x.tv_usec)/1000000.0; 
std::string remove_new_line(std::string in_string) {
    in_string.erase(std::remove(in_string.begin(), in_string.end(), '\n'), in_string.end());
    return in_string;

}

void summarize_outputs(const C_arg& c_inst_args, C_time& c_inst_time) {
    std::string last_ending_time;

    //--------------------------------------------------
    // stdout
    //--------------------------------------------------
    /*
    std::cout << "Running Time" << std::endl;
    std::cout << "     Checking reads" << std::endl;
    std::cout << "          Start: " << remove_new_line(c_inst_time.start_parse_args) << std::endl;
    std::cout << "          End  : " << remove_new_line(c_inst_time.end_check_read_file) << std::endl;

    if (c_inst_args.load_bf == false) {
        std::cout << "     Counting the number of unique solid k-mers" << std::endl;
        std::cout << "          Start: " << remove_new_line(c_inst_time.end_check_read_file) << std::endl;
        std::cout << "          End  : " << remove_new_line(c_inst_time.end_count_kmers) << std::endl;

        std::cout << "     Programming k-mers into the Bloom filter" << std::endl;
        std::cout << "          Start: " << remove_new_line(c_inst_time.end_count_kmers) << std::endl;
        std::cout << "          End  : " << remove_new_line(c_inst_time.end_program_kmers_into_bloom_filter) << std::endl;
        last_ending_time = c_inst_time.end_program_kmers_into_bloom_filter;
    }
    else {
        last_ending_time = c_inst_time.end_check_read_file;
    }

    std::cout << "     Correcting errors in reads" << std::endl;
    std::cout << "          Start: " << remove_new_line(last_ending_time) << std::endl;
    std::cout << "          End  : " << remove_new_line(c_inst_time.end_correct_errors_in_reads) << std::endl;

    std::cout << "     Merging output files" << std::endl;
    std::cout << "          Start: " << remove_new_line(c_inst_time.end_correct_errors_in_reads) << std::endl;
    std::cout << "          End  : " << remove_new_line(c_inst_time.end_merge_output_files) << std::endl;

    std::cout << std::endl;
    */
    std::cout << "BLESS is successfully completed" << std::endl << std::endl;
}
void ReadFileName(C_file_name& readname, C_arg& c_inst_args){


    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    uint32_t readlist_size = 0;
    FILE* readlist;
    if(my_rank == 0){
        std::cout << c_inst_args.read_list << std::endl;
        readlist = fopen(c_inst_args.read_list.c_str(), "r");
        if(readlist == NULL){
            std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.read_list << std::endl << std::endl;
            exit(EXIT_SUCCESS);
        }
        fseek(readlist, 0, SEEK_END);
        readname.read_list_size = ftell(readlist);
        fseek(readlist, 0, SEEK_SET);
    }
    MPI_Bcast(&(readname.read_list_size), 1, MPI_INT, 0, MPI_COMM_WORLD);
    readname.read_list = (char*)malloc((readname.read_list_size  + 10)* sizeof(char));
    if(my_rank == 0){
        fread(readname.read_list, sizeof(char), readname.read_list_size, readlist);
        fclose(readlist);     
    }

    MPI_Bcast(readname.read_list, readname.read_list_size + 4, MPI_CHAR, 0, MPI_COMM_WORLD);
    readname.read_file_num = 0;
    for(int i = 0; i < readname.read_list_size; i ++){
        if(readname.read_list[i] == '\n')
            readname.read_file_num += 1;
    } 
    readname.read_file_name = (char**)malloc((readname.read_file_num + 1) * sizeof(char*));
    readname.read_file_name[0] = readname.read_list;
    int read_file_num = 1;
    for(int i = 1; i < readname.read_list_size; i ++){ 
        if(readname.read_list[i] == '\n'){
            readname.read_list[i] = 0;
            readname.read_file_name[read_file_num ++] = readname.read_list + i + 1;
        }
    }
    readname.read_list[readname.read_list_size] = 0;
}

void correct_manager(init_args& arg, C_time& c_inst_time, C_arg& c_inst_args, C_check_read& c_inst_check_reads, C_count_solid_kmers& c_inst_count_solid_kmers, C_generate_bloom& c_inst_generate_bloom){




    if (c_inst_args.load_bf == false) {
        arg.bit_vector_width = c_inst_generate_bloom.bit_vector_width;
        arg.bit_vector_width_byte = c_inst_generate_bloom.bit_vector_width_byte;
        arg.num_hash_func = c_inst_generate_bloom.num_hash_func;
        arg.kmer_length = c_inst_generate_bloom.kmer_length;
        arg.kmer_occurrence_threshold = c_inst_generate_bloom.kmer_occurrence_threshold;
        arg.random_seed = c_inst_generate_bloom.random_seed;
        arg.num_unique_solid_kmers = c_inst_generate_bloom.num_unique_solid_kmers;
        arg.quality_score_offset = c_inst_generate_bloom.quality_score_offset;
        arg.quality_score_cutoff = c_inst_generate_bloom.quality_score_cutoff;
        arg.extremely_low_quality_score = c_inst_generate_bloom.extremely_low_quality_score;
        arg.max_extension = c_inst_generate_bloom.max_extension;
        arg.read_num = c_inst_generate_bloom.num_reads;
        arg.notrim = c_inst_args.notrim;
    }
    // reuse existing bloom filter data
    else {
        // read the bloom filter parameters
        // data will be loaded in load_split_bloom_filter later
        // but this is still needed to load some parameters

        arg.bit_vector_width = c_inst_generate_bloom.bit_vector_width;
        arg.bit_vector_width_byte = c_inst_generate_bloom.bit_vector_width_byte;
        arg.num_hash_func = c_inst_generate_bloom.num_hash_func;
        arg.kmer_length = c_inst_generate_bloom.kmer_length;
        arg.kmer_occurrence_threshold = c_inst_generate_bloom.kmer_occurrence_threshold;
        arg.random_seed = c_inst_generate_bloom.random_seed;
        arg.num_unique_solid_kmers = c_inst_generate_bloom.num_unique_solid_kmers;
        arg.quality_score_offset = c_inst_generate_bloom.quality_score_offset;
        arg.quality_score_cutoff = c_inst_generate_bloom.quality_score_cutoff;
        arg.extremely_low_quality_score = c_inst_generate_bloom.extremely_low_quality_score;
        arg.max_extension = c_inst_generate_bloom.max_extension;
        arg.read_num = c_inst_generate_bloom.num_reads;
        arg.notrim = c_inst_args.notrim;
    }

    arg.bit_vector = (uint8_t*)malloc(arg.bit_vector_width_byte * sizeof(uint8_t));

    std::ifstream f_bf_data;
    f_bf_data.open(c_inst_args.bf_data_file_name.c_str(), std::ios::binary);
    FILE* pfile;
    char file_name[200];
    sprintf(file_name, "%s", c_inst_args.bf_data_file_name.c_str());
    pfile = fopen(file_name, "rb");
    if(pfile == NULL)
    printf("file is null\n");
    //uint32_t file_size = fread(arg.bit_vector, 1, arg.bit_vector_width_byte, pfile);
    f_bf_data.read((char*)(arg.bit_vector), arg.bit_vector_width_byte);
    f_bf_data.close();
    fclose(pfile);
    c_inst_generate_bloom.generate_hash_seed(c_inst_args.random_seed, &(arg.hash_seed));

}
void gen_prefix(C_arg& c_inst_args, C_file_name& readname,  string& init_prefix, int my_rank, int node_num, int process_iter){
        if(c_inst_args.paired_read == true){
            int rind =  my_rank / 2;
            rind *= 2;
            if(node_num > 1){
                rind += process_iter * node_num;
            }else{
                
                rind += process_iter * 2;
            }
            std::string name1="";
            c_inst_args.read_file_name1.clear(); 
            int cind = 0;
            int dotind = 0;
            //printf("%s\n", readname.read_file_name[rind]);
            while(1){
                char ct = readname.read_file_name[rind][cind++];
                if(ct == 0)
                    break;
                //printf("%c\n",ct);
                if(ct == '.')
                    dotind = cind;
                name1 += ct;
                c_inst_args.read_file_name1 += ct;
            }
            if(dotind == 0)
                dotind = name1.length();
            //std::cout << "name1 " << name1 << std::endl;
            //std::cout << name1.length() << std::endl;
            std::size_t found1 = name1.find_last_of('/');
            c_inst_args.prefix.clear();
            c_inst_args.prefix = init_prefix;
            //std::cout << init_prefix << std::endl; 
            //std::cout << dotind << std::endl;
            //std::cout << found1 << std::endl;
            //std::cout << name1.substr(found1, 60) << std::endl;
            if(found1 == name1.length())
                found1 = 0;
            else
                found1 += 1; 
            c_inst_args.prefix = c_inst_args.prefix + "."; 
            for(int i = found1; i < dotind; i ++)
                c_inst_args.prefix += readname.read_file_name[rind][i];
            std::string name2 = "";
            c_inst_args.read_file_name2.clear(); 
            cind = 0;
            dotind = 0; 
            while(1){
                char ct = readname.read_file_name[rind + 1][cind++];
                if(ct == 0)
                    break;
                if(ct == '.')
                    dotind = cind;
                name2 += ct;
                c_inst_args.read_file_name2 += ct;
            }
            if(dotind == 0)
                dotind = name2.length();
            found1 = name2.find_last_of('/');
            if(found1 == name2.length())
                found1 = 0;
            else
                found1 += 1; 
            for(int i = found1; i < dotind - 1; i ++)
                c_inst_args.prefix += readname.read_file_name[rind + 1][i];
            //c_inst_args.prefix = c_inst_args.prefix + "." + name2.substr(found1, dotind); 
            //std::cout << "args prefix：" << c_inst_args.prefix << std::endl;
            
        }else{
            int rind =  my_rank;
            rind += process_iter * node_num;
                
            //std::string name1(readname.read_file_name[rind]);
            std::string name1="";
            c_inst_args.read_file_name.clear(); 
            int cind = 0;
            int dotind = 0;
            while(1){
                char ct = readname.read_file_name[rind][cind++];
                if(ct == 0)
                    break;
                if(ct == '.')
                    dotind = cind;
                name1 += ct;
                c_inst_args.read_file_name += ct;
            }

            if(dotind == 0)
                dotind = name1.length();
            std::size_t found1 = name1.find_last_of("/");
            c_inst_args.prefix.clear();
            c_inst_args.prefix = init_prefix;
            if(found1 == name1.length())
                found1 = 0;
            else
                found1 += 1; 
            c_inst_args.prefix = c_inst_args.prefix + "." + name1.substr(found1, dotind); 
 

        }

//            //std::string name1(readname.read_file_name[rind]);
//            std::string name1="";
//            c_inst_args.read_file_name1.clear(); 
//            int cind = 0;
//            while(1){
//                char ct = readname.read_file_name[rind][cind++];
//                if(ct == 0)
//                    break;
//                name1 += ct;
//                c_inst_args.read_file_name1 += ct;
//            }
//
//            std::size_t found = name1.rfind(".");
//            c_inst_args.prefix.clear();
//            c_inst_args.prefix = init_prefix;
//            if(found != std::string::npos)
//                c_inst_args.prefix = c_inst_args.prefix + "." + name1.substr(0, found); 
//            std::string name2 = "";
//            c_inst_args.read_file_name2.clear(); 
//            cind = 0;
//            while(1){
//                char ct = readname.read_file_name[rind + 1][cind++];
//                if(ct == 0)
//                    break;
//                name2 += ct;
//                c_inst_args.read_file_name2 += ct;
//            }
//            found = name2.rfind(".");
//            if(found != std::string::npos)
//                c_inst_args.prefix = c_inst_args.prefix + "." + name2.substr(0, found); 
//            std::cout << "args prefix：" << c_inst_args.prefix << std::endl;
//            
 
}
int main(int argc, char* argv[]) {

    MPI_Init(&argc, &argv);	
    C_time c_inst_time;
    C_file_name readname;
    C_arg c_inst_args(argc, argv, c_inst_time, readname);



    int my_rank;
    int node_num;

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &node_num);    

    int group_id = my_rank / 2;    
    MPI_Comm group_comm;
    MPI_Comm_split(MPI_COMM_WORLD, group_id, my_rank, &group_comm);
    int group_rank;
    MPI_Comm_rank(group_comm, &group_rank);


    ReadFileName(readname, c_inst_args);
    
    int process_num = (readname.read_file_num + node_num - 1)/ node_num;
    if(node_num == 1)
        process_num = (readname.read_file_num + 1 )/ 2;

    double mpi_total_time = 0;
    std::string init_prefix = c_inst_args.prefix;
    //std::cout << "process_num " << process_num << std::endl;
    double mpi_time = 0;
    for(int i = 0; i < process_num; i ++){

        gen_prefix(c_inst_args, readname, init_prefix, my_rank, node_num, i);
       
        std::cout << c_inst_args.prefix << std::endl;

        //std::string test(readname.read_file_name[0]);
        //std::cout << test << std::endl;

        std::cout << c_inst_args.read_file_name1 << std::endl;
        std::cout << c_inst_args.read_file_name2 << std::endl;

        c_inst_args.prepare_args(readname);

        std::cout << "prepare done"  << std::endl;
        
        C_check_read c_inst_check_reads;
        c_inst_check_reads.check_read_file(c_inst_args, c_inst_time, group_rank, group_comm, node_num);
         
        MPI_Bcast(&(c_inst_check_reads.num_reads), 1, MPI_UNSIGNED_LONG_LONG, 0, group_comm);
        MPI_Bcast(&(c_inst_check_reads.quality_score_offset), 1, MPI_UNSIGNED_LONG_LONG, 0, group_comm);
        MPI_Bcast(&(c_inst_check_reads.quality_score_cutoff), 1, MPI_UNSIGNED_LONG_LONG, 0, group_comm);
        MPI_Bcast(&(c_inst_check_reads.extremely_low_quality_score), 1, MPI_UNSIGNED_LONG_LONG, 0, group_comm);
        C_generate_bloom c_inst_generate_bloom;
        C_count_solid_kmers c_inst_count_solid_kmers;

        c_inst_count_solid_kmers.kmer_length = c_inst_args.kmer_length;
        c_inst_count_solid_kmers.rank_node = my_rank;
        
        if (c_inst_args.load_bf == false) {
             c_inst_count_solid_kmers.kmer_occurrence_threshold =
                   c_inst_args.kmer_occurrence_threshold;
            // count k-mers
            c_inst_count_solid_kmers.count_kmers(c_inst_args, c_inst_time, group_rank, group_comm, node_num);

        }else{
            c_inst_args.bf_prefix = c_inst_args.prefix;
        }

        MPI_Bcast(&(c_inst_count_solid_kmers.num_unique_solid_kmers), 1, MPI_UNSIGNED_LONG_LONG, 0,  group_comm);    
        MPI_Bcast(&(c_inst_count_solid_kmers.kmer_occurrence_threshold), 1, MPI_UNSIGNED_LONG_LONG, 0,  group_comm);    
        
        c_inst_generate_bloom.rank_node                   = my_rank;
        c_inst_generate_bloom.rank_smp                    = 1;
        c_inst_generate_bloom.size_node                   = node_num;
        c_inst_generate_bloom.group_rank                  = group_rank;
        c_inst_generate_bloom.max_extension               = c_inst_args.extend;
        c_inst_generate_bloom.quality_score_offset        = c_inst_check_reads.quality_score_offset;
        c_inst_generate_bloom.quality_score_cutoff        = c_inst_check_reads.quality_score_cutoff;
        c_inst_generate_bloom.extremely_low_quality_score = c_inst_check_reads.extremely_low_quality_score;
        c_inst_generate_bloom.num_reads                   = c_inst_check_reads.num_reads;

        c_inst_generate_bloom.kmer_length = c_inst_args.kmer_length;

        struct timeval s_time, e_time;
        GETTIME(s_time);
        init_args arg;
        // newly generate bloom filter data
        if (c_inst_args.load_bf == false) {
            c_inst_generate_bloom.kmer_occurrence_threshold = c_inst_count_solid_kmers.kmer_occurrence_threshold;
            c_inst_generate_bloom.num_unique_solid_kmers    = (uint32_t)(c_inst_count_solid_kmers.num_unique_solid_kmers);
            // program k-mers into a Bloom filter
            c_inst_generate_bloom.determine_bloom_filter_parameters(c_inst_args, c_inst_time);
            // program solid k-mers into the bloom filter
            c_inst_generate_bloom.program_kmers(arg, c_inst_args, c_inst_time, group_comm);
            //std::cout << "bit_vector_witdth" << std::endl;
            //std::cout << c_inst_generate_bloom.bit_vector_width << std::endl;
        }else{
            c_inst_generate_bloom.read_bloom_filter_parameters(c_inst_args, c_inst_time);
        }

        correct_manager(arg, c_inst_time, c_inst_args, c_inst_check_reads, c_inst_count_solid_kmers, c_inst_generate_bloom);
	//for(int i = 0; i < arg.num_hash_func; i ++){
	//	printf("hash func %d\n", arg.hash_seed[i]);
	//}
	
        
        
        if(group_rank == 0 || node_num == 1){ 
            uint8_t* bit_vector = arg.bit_vector;
            //printf("%c\n", bit_vector[50000]); 
            std::cout << c_inst_args.read_file_name1 << std::endl;    
            //printf("%d %d %d %d %d %d %d %d %d %d %d %d %d\n", arg.bit_vector_width, arg.bit_vector_width_byte, arg.num_hash_func, arg.kmer_length, arg.kmer_occurrence_threshold, arg.random_seed,arg.num_unique_solid_kmers, arg.read_num, arg.read_len, arg.quality_score_offset, arg.quality_score_cutoff, arg.extremely_low_quality_score, arg.max_extension); 
            FILE* pfile = fopen(c_inst_args.read_file_name1.c_str(), "r");
            FILE* f_corrected_read = fopen(c_inst_args.corrected_read_file_name1.c_str(), "w");
            //about file
            struct stat sb;
            int fmapread, fmapwrite;
            fmapread = open(c_inst_args.read_file_name1.c_str(),O_RDONLY);
            fstat(fmapread, &sb);
            close(fmapread);

            manage_data(arg, c_inst_time, c_inst_args, bit_vector, my_rank, sb, pfile, f_corrected_read, mpi_total_time);
            fclose(pfile);
            fclose(f_corrected_read);
        }
        if(group_rank == 1 || node_num == 1){ 
            //printf("%d %d %d %d %d %d %d %d %d %d %d %d %d\n", arg.bit_vector_width, arg.bit_vector_width_byte, arg.num_hash_func, arg.kmer_length, arg.kmer_occurrence_threshold, arg.random_seed,arg.num_unique_solid_kmers, arg.read_num, arg.read_len, arg.quality_score_offset, arg.quality_score_cutoff, arg.extremely_low_quality_score, arg.max_extension); 
            uint8_t* bit_vector = arg.bit_vector;
            //printf("%c\n", bit_vector[50000]); 
            std::cout << c_inst_args.read_file_name2 << std::endl;    
            FILE* pfile = fopen(c_inst_args.read_file_name2.c_str(), "r");
            FILE* f_corrected_read = fopen(c_inst_args.corrected_read_file_name2.c_str(), "w");
            //about file
            struct stat sb;
            int fmapread, fmapwrite;
            fmapread = open(c_inst_args.read_file_name2.c_str(),O_RDONLY);
            fstat(fmapread, &sb);
            close(fmapread);

            manage_data(arg, c_inst_time, c_inst_args, bit_vector, my_rank, sb, pfile, f_corrected_read, mpi_total_time);
            fclose(pfile);
            fclose(f_corrected_read);
        }
        
         
        free(arg.bit_vector);
        free(arg.hash_seed);
        
        MPI_Barrier(MPI_COMM_WORLD); 

        GETTIME(e_time);
        mpi_time += DIFFTIME(s_time, e_time);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    //printf("my rank == %d, mpi total time %lf %lf\n", my_rank, mpi_total_time, mpi_time);

    if(my_rank == 0)
        summarize_outputs(c_inst_args, c_inst_time);
    MPI_Finalize();
    return 0;

}
