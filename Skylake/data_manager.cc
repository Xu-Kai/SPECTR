#include "data_manager.h"
#include "mic_correct_errors.h"
#include "cpu_correct_errors.h"
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <semaphore.h>
#include <unistd.h>
#include <fcntl.h>
#include "new_query_text.h"
#include <pthread.h>
#define READ_MAX_LEN 256
#define BUFFER_NAME(x,y) x##y
#define MIC_COMPUTE_INIT(x) nocopy(BUFFER_NAME(mic_buffer, x):length(mic_buffer_size) alloc_if(0) free_if(0))\
                in(read_offset:length(0) alloc_if(0) free_if(0))\
                in(bit_vector:length(0) alloc_if(0) free_if(0))\
                in(hash_seed:length(0) alloc_if(0) free_if(0))\
                in(hash_array, notrim)\
                in(mic_offload_size, read_offset_size)\
                inout(mic_out_size, statistics)
#define GETTIME(x) gettimeofday(&x, NULL)
#define DIFFTIME(x,y) y.tv_sec-x.tv_sec + (y.tv_usec - x.tv_usec)/1000000.0; 
pthread_mutex_t malloc_lock;
void offset_data(char* buffer_list, uint64_t buffer_size, uint64_t& buffer_list_size){
 
    if(buffer_size < buffer_list_size + 4096){
            uint32_t n_num = 0;
            uint64_t offset_pos[8];
            uint64_t i = buffer_size - 1;
            if(buffer_size > buffer_list_size)
                i = buffer_list_size - 1;
            uint64_t j;
            uint32_t plus_num = 0;
            for(j = 0; j < 2048; j ++){
                if(buffer_list[i - j] == '\n'){
                     offset_pos[n_num] = j;
                     n_num += 1;
                     if(buffer_list[i - j + 1] == '+'){
                        plus_num = n_num;
                     }
                     if(n_num - plus_num == 2 && buffer_list[i - j + 1] == '@'){ 
                        if(n_num >= 5)
                           break;
                     }
                 }   
             }
            buffer_list_size = buffer_size - offset_pos[n_num - 5];

    }
     
            
}
__ONMIC__ void mic_correct(char* read_content , uint32_t data_size, uint32_t* out_size, uint32_t* read_offset, uint32_t read_offset_size, init_args& arg, uint8_t* bit_vector, uint32_t* hash_seed, stats_info& statistics){
    
    uint32_t num_threads = omp_get_num_procs() - 4;

    uint32_t read_offset_index[256];
    uint32_t read_num[256];
  while(2048 * num_threads > data_size){
    num_threads /= 2; 
  }
  if(num_threads == 0)
      num_threads = 1;

#pragma omp parallel for 	
    for(uint32_t i = 0; i < num_threads; i ++){
        read_offset_index[i] = i * (data_size / num_threads);
        if(i != 0){
            uint32_t n_num = 0;
            uint64_t offset_pos[8];
            uint64_t ii = read_offset_index[i];
            uint64_t j;
            uint32_t plus_num = 0;
            for(j = 0; j < 2048; j ++){
                if(read_content[ii - j] == '\n'){
                     offset_pos[n_num] = j;
                     n_num += 1;
                     if(read_content[ii - j + 1] == '+'){
                        plus_num = n_num;
                     }
                     if(n_num - plus_num == 2 && read_content[ii - j + 1] == '@'){ 
                        if(n_num >= 5)
                           break;
                     }
                 }   
             }
            read_offset_index[i] = ii - offset_pos[n_num - 5] + 1;
        }
    }
    

   uint32_t read_offset_thread_size = read_offset_size / num_threads;
#pragma omp parallel for
    for(uint32_t i = 0; i < num_threads; i ++){
        read_num[i] = 0;
        uint32_t offset_begin = read_offset_index[i];
        uint32_t offset_end = data_size;
        uint32_t *read_offset_tmp = read_offset +  i * read_offset_thread_size;;
        if(i + 1 < num_threads)
            offset_end = read_offset_index[i + 1];	   
        uint32_t line_num = 0;
        uint32_t read_length = 0;
        uint32_t read_ns = 0;
        
        for(uint32_t j = offset_begin; j < offset_end; j ++){
            if(read_content[j] == '\n'){
                line_num += 1;
                if((line_num & 3) == 1){
                    //if(read_num[i] * 4 < ttt)
                    read_offset_tmp[(read_num[i] << 2) ] = j + 1;  
                    read_length = 0;
                    read_ns = 0;
                }else if((line_num & 3) == 2){	
                    read_offset_tmp[(read_num[i] << 2) + 2] = read_length;  
                    read_offset_tmp[(read_num[i] << 2) + 3] = read_ns;
                }else if((line_num & 3) == 3){	
                    read_offset_tmp[(read_num[i] << 2) + 1] = j + 1;  
                    j += read_offset_tmp[(read_num[i] << 2) + 2];
                    read_num[i] += 1;	
                    if(read_num[i] > read_offset_thread_size){
                        printf("read is too many\n");
                        break;

                    }
                }

            } else if((line_num & 3) == 1){
                read_length += 1;
                read_content[j] = toupper(read_content[j]);
                if(read_content[j] != 'A' && read_content[j] != 'C' && read_content[j] != 'G' && read_content[j] != 'T'){
                    read_ns += 1;
                    read_content[j] = 'A';
                }
            }
        }
        
        if(read_num[i] * 4 > read_offset_thread_size)
            printf("read is too many!! %d %d\n", read_offset_thread_size, read_num[i]);
        //printf("i iter done %d\n", i);
    }



    uint32_t total_read_num = read_num[0];
    for(uint32_t i = 1; i < num_threads; i ++){
        memcpy(read_offset + (total_read_num * 4), read_offset + i * read_offset_thread_size,  (read_num[i] * 4) * sizeof(uint32_t));
        total_read_num += read_num[i];
    }

    func_mem funcm[num_threads];
    for (int i = 0; i < num_threads; i++) {
        funcm[i].ubuffer = (uint32_t*)_mm_malloc((READ_MAX_LEN * 4 * 2 + 76 * READ_MAX_LEN * 3 + 2 * READ_MAX_LEN + 64 * 75 * READ_MAX_LEN) * sizeof(uint32_t), 64);
        funcm[i].cbuffer = (char*)_mm_malloc(6 * READ_MAX_LEN + 76 * READ_MAX_LEN * 3  + 2 * READ_MAX_LEN + 64 * 75 * READ_MAX_LEN, 64);
    }
    uint32_t num_trimmed_bases_tmp = 0;
    uint32_t num_corrected_errors_local = 0;
    uint32_t trim_5_end = 0;
    uint32_t trim_3_end = 0;
    uint32_t max_trimmed_bases = 0;
    uint32_t min_check_length;
    uint32_t max_allowed_ns;
    bool too_many_errors;
    uint32_t num_corrected_errors_tmp = 0;
    uint32_t num_corrected_reads_tmp = 0;

#pragma omp parallel for schedule(dynamic, 10)  private( min_check_length, max_allowed_ns,\
        max_trimmed_bases, num_corrected_errors_local, trim_5_end, trim_3_end, too_many_errors) \
        reduction(+: num_corrected_errors_tmp, num_corrected_reads_tmp, num_trimmed_bases_tmp)  num_threads(num_threads)
    for (uint32_t i = 0; i < total_read_num; i++) {
        uint32_t read_length = read_offset[(i << 2) + 2];
        min_check_length = read_length * CHECK_RANGE_RATIO;
        max_allowed_ns = read_length * MAX_N_RATIO;
        // forward
        // too short read: no trimming
        if (read_length <= MIN_BASES_AFTER_TRIMMING) {
            max_trimmed_bases = 0;
        } else {
            uint32_t rmax = read_length * MAX_TRIMMING_RATE;
            uint32_t lmax = read_length - MIN_BASES_AFTER_TRIMMING;
            max_trimmed_bases =  rmax  <lmax?rmax:lmax;
        }
        // initialize variables
        num_corrected_errors_local = 0;
        trim_5_end = 0;
        trim_3_end = 0;
        int ind = omp_get_thread_num();
        char read_modification[256];
        memset(read_modification, '0', 256);
        funcm[ind].kmer_ind = 0;
        funcm[ind].sequence_ind = 0;
        funcm[ind].ubuffer_offset = 0;
        funcm[ind].cbuffer_offset = 0;

        //char *kmer = funcm[ind].cbuffer;
        if(read_offset[(i << 2) + 3] <=  max_allowed_ns && read_offset[(i << 2) + 2] >= arg.kmer_length){
             
            mic_correct_errors(&read_content[read_offset[(i << 2)]],
                    read_content + read_offset[(i << 2) + 1],
                    read_modification , hash_seed, bit_vector,
                    arg, num_corrected_errors_local, trim_5_end, trim_3_end,
                    max_trimmed_bases, min_check_length, read_length, i,
                    funcm[ind]);
            
            
            // no trim

            if (arg.notrim == true) {
                trim_5_end = 0;
                trim_3_end = 0;
            }		
            // adjust the number of trimmed bases
            else {
                if ((trim_5_end + trim_3_end) > max_trimmed_bases) {
                    if (trim_3_end <= max_trimmed_bases) {
                        trim_5_end = 0;
                    }
                    else if (trim_5_end <= max_trimmed_bases) {
                        trim_3_end = 0;
                    }
                    else {
                        trim_5_end = 0;
                        trim_3_end = max_trimmed_bases;
                    }
                }
            }

            num_trimmed_bases_tmp += (trim_5_end + trim_3_end);

            // update num_corrected_reads
            too_many_errors = false;
            if (num_corrected_errors_local > (read_length * MAX_ERROR_RATE)) {
                too_many_errors = true;
            }
            else if (num_corrected_errors_local > 0) {
                num_corrected_errors_tmp += num_corrected_errors_local;

                num_corrected_reads_tmp++;
            }
            else if (arg.notrim == false) {
                if ((trim_5_end > 0) || (trim_3_end > 0)) {
                    num_corrected_reads_tmp++;
                }
            }

            // make a corrected read
            if (too_many_errors == false) {
                // apply modifications to the read
                //if(read_length  - trim_3_end> 101 )
                //    printf("this is wrong\n");
                for (unsigned int it_base = trim_5_end; it_base < (read_length - trim_3_end); it_base++) {
                    if (read_modification[it_base] != '0') {
                        read_content[read_offset[i << 2] + it_base] = read_modification[it_base];
                    }
                }
                for(unsigned int it_base = 0; it_base < trim_5_end; it_base ++){
                    read_content[read_offset[i << 2] + it_base] = '0';
                    read_content[read_offset[(i << 2) + 1] + it_base] = '0';
                    
                }
                for(unsigned int it_base = read_length - trim_3_end; it_base < read_length; it_base ++){
                    read_content[read_offset[i << 2] + it_base] = '0';
                    read_content[read_offset[(i << 2) + 1] + it_base] = '0';
                }
            }
            

        }
        
    }
    
   //printf("num_corrected errors tmp %d num corrected_reads %d, num_trimmed bases == %d\n", num_corrected_errors_tmp, num_corrected_reads_tmp, num_trimmed_bases_tmp); 
   


#pragma omp parallel for  
    for(uint32_t i = 0; i < num_threads; i ++){
        uint32_t begin = i * (data_size / num_threads);
        uint32_t end = (i + 1) * (data_size / num_threads);
        if(i == num_threads - 1)
            end = data_size;
        uint32_t r_index = begin;
        for(uint32_t j = begin; j < end; j ++){
            if(read_content[j] != 0){
                if(r_index != j)
                    read_content[r_index] = read_content[j];
                r_index ++;
            }
        }
        read_offset_index[i] = r_index;
    }
    
    uint32_t total_offset = read_offset_index[0]; 
    for(uint32_t i = 1; i < num_threads; i ++){
        uint32_t begin = i * (data_size / num_threads);
        uint32_t block_length = read_offset_index[i] - begin;
        if(total_offset != begin)
            memcpy(read_content + total_offset, read_content + begin, block_length * sizeof(char));
        total_offset += block_length;
    } 

    *out_size = total_offset;

    
    statistics.num_corrected_errors += num_corrected_errors_tmp;
    statistics.num_corrected_reads += num_corrected_reads_tmp;
    statistics.num_trimmed_bases += num_trimmed_bases_tmp;    

  //printf("num_corrected errors tmp %d num corrected_reads %d, num_trimmed bases == %d\n", num_corrected_errors_tmp, num_corrected_reads_tmp, num_trimmed_bases_tmp); 
    for (int i = 0; i < num_threads; i++) {
        _mm_free(funcm[i].ubuffer);
        _mm_free(funcm[i].cbuffer);
    }
    
}

  uint64_t total_reads_num1 = 0; 
  uint64_t total_reads_num2 = 0; 
  uint64_t total_reads_num3 = 0; 

void cpu_correct(char* read_content , uint32_t data_size, uint32_t *out_data_size, uint32_t* read_offset, uint32_t read_offset_size, init_args& arg, uint8_t* bit_vector, stats_info& statistics){

   int num_threads_exe = omp_get_num_procs();
   num_threads_exe = omp_get_max_threads();
   if(num_threads_exe > 8)
       num_threads_exe -= 8;
   uint32_t read_offset_index[300];
   uint32_t read_num[300];
   if(data_size < num_threads_exe * 2048)
   	num_threads_exe = 1; 
   uint32_t ii = 0;
   for(uint32_t i = 0; i < data_size;i ++ )
       if(read_content[i] == '\n')
            ii ++;
    printf("total number == %d\n", ii); 

#pragma omp parallel for num_threads(num_threads_exe)	
    for(uint32_t i = 0; i < num_threads_exe; i ++){
        read_offset_index[i] = i * (data_size / num_threads_exe);
        if(i != 0){
            uint32_t n_num = 0;
            uint32_t offset_pos[8];
            uint32_t ii = read_offset_index[i];
            uint32_t j;
            uint32_t plus_num = 0;
            //printf("begin compute\n");
           // printf("data size %d \n", data_size); 
            
            for(j = 0; j < 2048; j ++){
                if(read_content[ii - j] == '\n'){
                     offset_pos[n_num] = j;
                     n_num += 1;
                     if(read_content[ii - j + 1] == '+'){
                        plus_num = n_num;
                     }
                     if(n_num - plus_num == 2 && read_content[ii - j + 1] == '@'){ 
                        if(n_num >= 5)
                           break;
                     }
                 }   
             }
            read_offset_index[i] = ii - offset_pos[n_num - 5] + 1;
        }
        
    }
    
   //if(my_rank == 3)
    //   printf("middle\n");
    
    uint32_t read_offset_thread_size = read_offset_size / num_threads_exe;
#pragma omp parallel for num_threads(num_threads_exe)
    for(uint32_t i = 0; i < num_threads_exe; i ++){
        read_num[i] = 0;
        uint32_t offset_begin = read_offset_index[i];
        uint32_t offset_end = data_size;
        uint32_t *read_offset_tmp = read_offset +  i * read_offset_thread_size;;
        if(i + 1 < num_threads_exe)
            offset_end = read_offset_index[i + 1];	   
        uint32_t line_num = 0;
        uint32_t read_length = 0;
        uint32_t read_ns = 0;
        
        for(uint32_t j = offset_begin; j < offset_end; j ++){
            if(read_content[j] == '\n'){
                line_num += 1;
                if((line_num & 3) == 1){
                    //if(read_num[i] * 4 < ttt)
                    read_offset_tmp[(read_num[i] << 2) ] = j + 1;  
                    read_length = 0;
                    read_ns = 0;
                }else if((line_num & 3) == 2){	
                    read_offset_tmp[(read_num[i] << 2) + 2] = read_length;  
                    read_offset_tmp[(read_num[i] << 2) + 3] = read_ns;
                }else if((line_num & 3) == 3){	
                    read_offset_tmp[(read_num[i] << 2) + 1] = j + 1;  
                    j += read_offset_tmp[(read_num[i] << 2) + 2];
                    read_num[i] += 1;	
                    if(read_num[i] > read_offset_thread_size){
                        printf("read is too many\n");
                        break;

                    }
                }

            } else if((line_num & 3) == 1){
                read_length += 1;
                read_content[j] = toupper(read_content[j]);
                if(read_content[j] != 'A' && read_content[j] != 'C' && read_content[j] != 'G' && read_content[j] != 'T'){
                    read_ns += 1;
                    read_content[j] = 'A';
                }
            }
        }
        
        if(read_num[i] * 4 > read_offset_thread_size)
            printf("read is too many!! %d %d\n", read_offset_thread_size, read_num[i]);
    }
    
   //  if(my_rank == 3)
     //  printf("end\n");
    uint32_t total_read_num = read_num[0];
    for(uint32_t i = 1; i < num_threads_exe; i ++){
        memcpy(read_offset + (total_read_num << 2), read_offset + i * read_offset_thread_size, (read_num[i]  << 2) * sizeof(uint32_t));
        total_read_num += read_num[i];
    }
    //printf("total read num == %d\n", total_read_num);
    
    num_threads_exe = omp_get_max_threads();
    if(num_threads_exe > 8)
        num_threads_exe -= 8;
    uint32_t* hash_seed = arg.hash_seed;
    func_mem funcm[300];
    for (int i = 0; i < num_threads_exe; i++) {
        funcm[i].ubuffer = (uint32_t*)_mm_malloc((READ_MAX_LEN * 4 * 2 + 76 * READ_MAX_LEN * 3 + 2 * READ_MAX_LEN + 64 * 75 * READ_MAX_LEN) * sizeof(uint32_t), 64);
        funcm[i].cbuffer = (char*)malloc(6 * READ_MAX_LEN + 76 * READ_MAX_LEN * 3  + 2 * READ_MAX_LEN + 64 * 75 * READ_MAX_LEN);
    }
    uint32_t num_trimmed_bases_tmp = 0;
    uint32_t num_corrected_errors_local = 0;
    uint32_t trim_5_end = 0;
    uint32_t trim_3_end = 0;
    uint32_t max_trimmed_bases = 0;
    uint32_t min_check_length;
    uint32_t max_allowed_ns;
    bool too_many_errors;
    uint32_t num_corrected_errors_tmp = 0;
    uint32_t num_corrected_reads_tmp = 0;
    
    //printf("%d %d %d %d %d\n", read_num[0], read_offset[t << 2], read_offset[(t << 2) + 1], read_offset[(t << 2) + 2], read_offset[(t << 2) + 3]);

    total_reads_num1 += total_read_num; 
#pragma omp parallel for schedule(dynamic, 10)  private( min_check_length, max_allowed_ns,\
        max_trimmed_bases, num_corrected_errors_local, trim_5_end, trim_3_end, too_many_errors) \
        reduction(+: num_corrected_errors_tmp, num_corrected_reads_tmp, num_trimmed_bases_tmp, total_reads_num2, total_reads_num3) num_threads(1)
    for (uint32_t i = 0; i < total_read_num; i++) {
        int th = omp_get_num_threads();
        total_reads_num2 += 1;
        //printf("threaf num == %d \n", th);
        uint32_t read_length = read_offset[(i << 2) + 2];
        min_check_length = read_length * CHECK_RANGE_RATIO;
        max_allowed_ns = read_length * MAX_N_RATIO;
        // forward
        // too short read: no trimming
        if (read_length <= MIN_BASES_AFTER_TRIMMING) {
            max_trimmed_bases = 0;
        } else {
            uint32_t rmax = read_length * MAX_TRIMMING_RATE;
            uint32_t lmax = read_length - MIN_BASES_AFTER_TRIMMING;
            max_trimmed_bases =  rmax  <lmax?rmax:lmax;
        }
        // initialize variables
        num_corrected_errors_local = 0;
        trim_5_end = 0;
        trim_3_end = 0;
        int ind = omp_get_thread_num();
        char read_modification[256];
        memset(read_modification, '0', 256);
        funcm[ind].kmer_ind = 0;
        funcm[ind].sequence_ind = 0;
        funcm[ind].ubuffer_offset = 0;
        funcm[ind].cbuffer_offset = 0;

        //char *kmer = funcm[ind].cbuffer;
        if(read_offset[(i << 2) + 3] <=  max_allowed_ns && read_offset[(i << 2) + 2] >= arg.kmer_length){
            //if(i == 58) 
             
            total_reads_num3 += 1;
            
            cpu_correct_errors(read_content + read_offset[(i << 2)],
                    read_content + read_offset[(i << 2) + 1],
                    read_modification , hash_seed, bit_vector,
                    arg, num_corrected_errors_local, trim_5_end, trim_3_end,
                    max_trimmed_bases, min_check_length, read_length, i,
                    funcm[ind]);
             
            // no trim
            if (arg.notrim == true) {
                trim_5_end = 0;
                trim_3_end = 0;
            }		
            // adjust the number of trimmed bases
            else {
                if ((trim_5_end + trim_3_end) > max_trimmed_bases) {
                    if (trim_3_end <= max_trimmed_bases) {
                        trim_5_end = 0;
                    }
                    else if (trim_5_end <= max_trimmed_bases) {
                        trim_3_end = 0;
                    }
                    else {
                        trim_5_end = 0;
                        trim_3_end = max_trimmed_bases;
                    }
                }
            }

            num_trimmed_bases_tmp += (trim_5_end + trim_3_end);

            // update num_corrected_reads
            too_many_errors = false;
            if (num_corrected_errors_local > (read_length * MAX_ERROR_RATE)) {
                too_many_errors = true;
                //printf("too many errors\n");
            }
            else if (num_corrected_errors_local > 0) {
                num_corrected_errors_tmp += num_corrected_errors_local;

                num_corrected_reads_tmp++;
            }
            else if (arg.notrim == false) {
                if ((trim_5_end > 0) || (trim_3_end > 0)) {
                    num_corrected_reads_tmp++;
                }
            }

            // make a corrected read
            if (too_many_errors == false) {
                // apply modifications to the read
               for (unsigned int it_base = trim_5_end; it_base < (read_length - trim_3_end); it_base++) {
                    if (read_modification[it_base] != '0') {
                        read_content[read_offset[i << 2] + it_base] = read_modification[it_base];
                    }
                }
                for(unsigned int it_base = 0; it_base < trim_5_end; it_base ++){
                    read_content[read_offset[i << 2] + it_base] = '0';
                    read_content[read_offset[(i << 2) + 1] + it_base] = '0';
                }
                for(unsigned int it_base = read_length - trim_3_end; it_base < read_length; it_base ++){
                    read_content[read_offset[i << 2] + it_base] = '0';
                    read_content[read_offset[(i << 2) + 1] + it_base] = '0';
                }

            }

        }

    }

    printf("total read num is == %ld %ld %ld\n", total_reads_num1, total_reads_num2, total_reads_num3);
    ii = 0;
 for(uint32_t i = 0; i < data_size;i ++ )
       if(read_content[i] == '\n')
            ii ++;
    printf("total number == %d\n", ii); 

    ii = 0;
#pragma omp parallel for num_threads(num_threads_exe) reduction(+: ii)	
    for(uint32_t i = 0; i < num_threads_exe; i ++){
        uint32_t begin = i * (data_size / num_threads_exe);
        uint32_t end = (i + 1) * (data_size / num_threads_exe);
        if(i == num_threads_exe - 1)
            end = data_size;
        uint32_t r_index = begin;
        for(uint32_t j = begin; j < end; j ++){
            if(read_content[j] != '0'){
                if(r_index != j)
                    read_content[r_index] = read_content[j];
                r_index ++;
            }
        }
        read_offset_index[i] = r_index;
        for(uint32_t j = begin;j < r_index;j ++ )
            if(read_content[j] == '\n')
                ii ++;
    }
    printf("total number == %d\n", ii); 


    uint32_t total_offset = read_offset_index[0]; 
    for(uint32_t i = 1; i < num_threads_exe; i ++){
        uint32_t begin = i * (data_size / num_threads_exe);
        uint32_t block_length = read_offset_index[i] - begin;
        for(uint32_t j = 0; j < block_length; j ++)
            read_content[total_offset ++] = read_content[begin + j];  
        //if(total_offset != begin)
         //   memcpy(read_content + total_offset, read_content + begin, block_length * sizeof(char));
        //total_offset += block_length;
    } 
    ii = 0;
    for(uint32_t i = 0; i < total_offset;i ++ )
       if(read_content[i] == '\n')
            ii ++;
    printf("total number == %d\n", ii); 
    *out_data_size = total_offset;
   
    //printf("total read num is == %d\n", num_corrected_reads_tmp);
  //  *out_data_size = data_size;
    statistics.num_corrected_errors += num_corrected_errors_tmp;
    statistics.num_corrected_reads += num_corrected_reads_tmp;
    statistics.num_trimmed_bases += num_trimmed_bases_tmp;    
  //printf("num_corrected errors tmp %d num corrected_reads %d, num_trimmed bases == %d\n", num_corrected_errors_tmp, num_corrected_reads_tmp, num_trimmed_bases_tmp); 
    for (int i = 0; i < num_threads_exe; i++) {
        _mm_free(funcm[i].ubuffer);
        free(funcm[i].cbuffer);
    }
    
}
void args_2_array(uint32_t* hash_array, init_args& arg){
    hash_array[0] = arg.bit_vector_width;
    hash_array[1] = arg.bit_vector_width_byte;
    hash_array[2] = arg.num_hash_func;
    hash_array[3] = arg.kmer_length;
    hash_array[4] = arg.kmer_occurrence_threshold;
    hash_array[5] = arg.random_seed;
    hash_array[6] = arg.num_unique_solid_kmers;
    hash_array[7] = arg.read_num;
    hash_array[8] = arg.read_len;
    hash_array[9] = arg.quality_score_offset;
    hash_array[10] = arg.quality_score_cutoff;
    hash_array[11] = arg.extremely_low_quality_score;
    hash_array[12] = arg.max_extension;
}
__ONMIC__ void array_2_args(uint32_t* hash_array, init_args& arg){
    arg.bit_vector_width = hash_array[0];
    arg.bit_vector_width_byte = hash_array[1];
    arg.num_hash_func = hash_array[2];
    arg.kmer_length = hash_array[3];
    arg.kmer_occurrence_threshold = hash_array[4];
    arg.random_seed = hash_array[5];
    arg.num_unique_solid_kmers = hash_array[6];
    arg.read_num = hash_array[7];
    arg.read_len = hash_array[8];
    arg.quality_score_offset = hash_array[9];
    arg.quality_score_cutoff = hash_array[10];
    arg.extremely_low_quality_score = hash_array[11];
    arg.max_extension = hash_array[12];
}
void init_statistics(stats_info& statistics){ 
        statistics.num_corrected_errors = 0;
        statistics.num_corrected_reads = 0;
        statistics.num_trimmed_bases = 0;
}
void merge_statistics(stats_info& statistics_result, stats_info& statistics){
        statistics_result.num_corrected_errors += statistics.num_corrected_errors;
        statistics_result.num_corrected_reads += statistics.num_corrected_reads;
        statistics_result.num_trimmed_bases +=statistics.num_trimmed_bases;
}
void init_control_args(omp_lock_t &read_lock, omp_lock_t &write_lock, sem_t &read_sem_count, sem_t &write_sem_count, sem_t *compute_sem_count, uint64_t *device_data_size, uint32_t *read_flag, uint64_t *device_buffer_size, uint32_t mic_buffer_size, uint32_t cpu_buffer_size, uint32_t mic_num, uint32_t device_num){
    omp_init_lock(&read_lock);
    omp_init_lock(&write_lock);
    sem_init(&read_sem_count, 0, device_num);
    sem_init(&write_sem_count, 0, 0);
    for(uint32_t i = 0; i < device_num; i ++){
        sem_init(&compute_sem_count[i], 0, 0);
        device_data_size[i] = 0;
        read_flag[i] = i; 
        if(i < mic_num){
            device_buffer_size[i] = mic_buffer_size;
        }else{
            device_buffer_size[i] = cpu_buffer_size;
        }
    }

}


void process_data_to_compute(sem_t &read_sem_count, sem_t& rbuffer_count_sem, sem_t *compute_sem_count, char** device_compute_buffer, uint64_t *device_buffer_size, uint64_t *device_data_size, uint32_t *read_flag, char** rbuffer_list, uint64_t *rbuffer_data_size, uint64_t pfile_len, uint32_t mic_num, uint32_t device_num, double *read_time, double &read_total_time){
            int my_rank = 0;
            MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
            uint32_t device_id = 0;
            uint32_t iter_count = 0;
            uint64_t total_read_size = 0;
            uint64_t rbuffer_index = 0;
            uint64_t total_compute_size = 0;
            struct timeval s_time, e_time;
            struct timeval read_start_time, read_done_time;
            GETTIME(read_start_time);
            sem_wait(&rbuffer_count_sem);
            while(1){
                sem_wait(&read_sem_count);
                device_id = read_flag[iter_count];
                uint64_t remain_buffer_size = device_buffer_size[device_id] - device_data_size[device_id];

                GETTIME(s_time); 
                uint64_t device_buffer_size_t = rbuffer_data_size[rbuffer_index] - total_compute_size;
                offset_data(rbuffer_list[rbuffer_index] + total_compute_size, remain_buffer_size, device_buffer_size_t);
                
                memcpy(device_compute_buffer[device_id] + device_data_size[device_id] , rbuffer_list[rbuffer_index] + total_compute_size, sizeof(char) * device_buffer_size_t);

                total_read_size += device_buffer_size_t;
                total_compute_size += device_buffer_size_t;
                device_data_size[device_id] += device_buffer_size_t;

                GETTIME(e_time);
                read_time[device_id] += DIFFTIME(s_time, e_time); 
                if(total_read_size == pfile_len){
                    sem_post(&compute_sem_count[device_id]);
                    for(uint32_t inc = 0;  inc < device_num; inc ++){
                        uint32_t device_id_tmp = inc;
                        sem_post(&compute_sem_count[device_id_tmp]);
                    } 
                    free(rbuffer_list[rbuffer_index]); 
                    break;
                }
                if(total_compute_size < rbuffer_data_size[rbuffer_index]){
                    sem_post(&compute_sem_count[device_id]);
                    iter_count += 1;
                    
                }else{
                    sem_post(&read_sem_count);
                    total_compute_size = 0;
                    free(rbuffer_list[rbuffer_index]); 
                    sem_wait(&rbuffer_count_sem);
                    rbuffer_index += 1; 
                }

            }

            GETTIME(read_done_time);
            read_total_time = DIFFTIME(read_start_time, read_done_time);

}
void process_computed_data_to_buffer(sem_t& write_sem_count, uint64_t *trimmed_data_size, uint32_t *write_flag, char** wbuffer_list, uint64_t *computed_data_size, uint64_t  pfile_len, double* write_time, double &write_total_time, FILE* f_corrected_read){

            int my_rank = 0;
            MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
            uint32_t device_id = 0;
            uint32_t iter_count = 0;
            uint32_t wbuffer_index = 0;
            uint32_t total_out_data_size = 0;
            uint64_t readed_size = 0;
            struct timeval s_time, e_time;
            struct timeval write_start_time, write_done_time;
            GETTIME(write_start_time);
            while(1){
                sem_wait(&write_sem_count);
                GETTIME(s_time);
                device_id = write_flag[iter_count];
                
                fwrite(wbuffer_list[iter_count], sizeof(char), trimmed_data_size[iter_count], f_corrected_read);
                free(wbuffer_list[iter_count]);
                readed_size += computed_data_size[iter_count];
                iter_count ++;
                GETTIME(e_time);
                write_time[device_id] += DIFFTIME(s_time, e_time);
                if(readed_size == pfile_len)
                    break;

            }
            GETTIME(write_done_time);
            write_total_time = DIFFTIME(write_start_time, write_done_time);

}
void manage_data(init_args arg, C_time& c_inst_time, C_arg& c_inst_arg, uint8_t* bit_vector, int my_rank, struct stat& pfile_stat, FILE* pfile, FILE* f_corrected_read, double& mpi_total_time){



    uint32_t mic_buffer_size = 200 * 1024 * 1024;
    uint32_t cpu_buffer_size = 200 * 1024 * 1024;
    uint32_t mic_num = 0; 
    uint32_t read_offset_size = 20000000;
    char* cpu_buffer = (char*)malloc(cpu_buffer_size * sizeof(char));
    uint32_t *read_offset = (uint32_t*)malloc(read_offset_size * sizeof(uint32_t));
    char* mmic_buffer=(char*)malloc(mic_buffer_size * sizeof(char));
    uint32_t* hash_seed = arg.hash_seed; 

//#ifdef __INTEL_OFFLOAD
//    mic_num = _Offload_number_of_devices();
//#endif
   // mic_num = 1;
    //mic_num = 1;
    mic_num = 0;
    uint32_t device_num = mic_num + 1;
    //device_num = mic_num;
    device_num = 1;
    
    //printf("mic num == %d %d\n", mic_num, my_rank);
    char* device_compute_buffer[10];
    for(int i = 0; i < mic_num; i ++){

        device_compute_buffer[i] = (char*)malloc(mic_buffer_size * sizeof(char));
#ifdef __MIC__
#pragma offload target(mic:(i))\
        nocopy(read_offset:length(read_offset_size) alloc_if(1) free_if(0))\
        in(bit_vector:length(arg.bit_vector_width_byte) alloc_if(1) free_if(0))\
        in(hash_seed:length(arg.num_hash_func) alloc_if(1) free_if(0))
        {
        }
#pragma offload target(mic:(i))\
        nocopy(mmic_buffer:length(mic_buffer_size) alloc_if(1) free_if(0))
        {
        }
#endif
    }
    device_compute_buffer[mic_num] = (char*)malloc(cpu_buffer_size * sizeof(char));

    sem_t read_sem_count, compute_sem_count[device_num], write_sem_count;
    omp_lock_t read_lock, write_lock;



    uint64_t device_buffer_size[device_num], device_data_size[device_num];
    uint32_t read_flag[600];
    uint32_t write_flag[600];
    uint32_t read_iteration = device_num;
    uint32_t write_iteration = 0;
    init_control_args(read_lock, write_lock, read_sem_count, write_sem_count, compute_sem_count, device_data_size, read_flag, device_buffer_size, mic_buffer_size, cpu_buffer_size, mic_num, device_num);


    stats_info statistics_result;
    init_statistics(statistics_result);

    uint64_t pfile_len = 0; 
    fseek(pfile, 0, SEEK_END);
    pfile_len = ftell(pfile);
    fseek(pfile, 0, SEEK_SET);

    sem_t rbuffer_count_sem;
    char *rbuffer_list[400];
    uint64_t  trimmed_data_size[1200], computed_data_size[1200];
    char *wbuffer_list[1200];
    uint64_t rbuffer_data_size[1200];
    uint64_t rbuffer_compute_size[1200];
    uint64_t rbuffer_size = 2L * 1024L * 1024L * 1024L;
    bool notrim = false;

    sem_init(&rbuffer_count_sem, 0, 1);
    pthread_mutex_init(&malloc_lock, NULL);
    rbuffer_list[0] = (char*)malloc(rbuffer_size * sizeof(char));
    //wbuffer_list[0] = (char*)malloc(wbuffer_size * sizeof(char));
    rbuffer_data_size[0] = fread(rbuffer_list[0], sizeof(char), rbuffer_size, pfile);

    uint32_t hash_array[16];    
    args_2_array(hash_array, arg);


    double all_time[device_num], read_time[device_num], offload_time[device_num], compute_time[device_num], offload_out_time[device_num], write_time[device_num];
    double read_total_time, write_total_time;
    double total_time;
    for(int i = 0; i < device_num; i ++){
        all_time[i] = 0;
        read_time[i] = 0;
        offload_time[i] = 0;
        compute_time[i] = 0;
        offload_out_time[i] = 0;
        write_time[i] = 0;
    }
    struct timeval start_time, end_time;
    GETTIME(start_time);

    time_t rawtime;
    time(&rawtime);
    c_inst_time.start_correct_errors_in_reads = asctime(localtime(&rawtime));
    //printf("start compute %lld \n", rbuffer_data_size[0]);
    arg.notrim = false; 
    omp_set_nested(1);
#pragma omp parallel for num_threads(8) 
    for(int i = 0; i < device_num + 3; i ++){
        if(i  < device_num){
            uint32_t out_data_size;
            stats_info statistics;
            init_statistics(statistics);
            uint32_t compute_data_size;
            struct timeval s_time, e_time;
            char* compute_buffer = device_compute_buffer[i];
            char* dwrite_buffer;
            int iter_count = 0; 
            while(1){
                GETTIME(s_time);
                sem_wait(&compute_sem_count[i]);

                compute_data_size = device_data_size[i];	
                out_data_size = device_data_size[i];
                device_data_size[i] = 0;

                if(compute_data_size == 0)
                    break;
                if(i == mic_num){
                    //memcpy(cpu_buffer, compute_buffer, compute_data_size * sizeof(char));
                    char *tmp_buffer = cpu_buffer;
                    cpu_buffer = device_compute_buffer[i];    
                    device_compute_buffer[i] = tmp_buffer;
                
                }
                else{
#ifdef __MIC__
                    #pragma offload target(mic:i)\
                        in(compute_buffer[0: mic_buffer_size] : into(mmic_buffer[0:mic_buffer_size]) alloc_if(0) free_if(0))
                        {
                        }
#endif
                }
                     
                omp_set_lock(&read_lock);
                read_flag[read_iteration ++] = i;
                omp_unset_lock(&read_lock);

                sem_post(&read_sem_count);

                GETTIME(e_time);
                offload_time[i] += DIFFTIME(s_time, e_time);

                GETTIME(s_time);
                //printf("iter count %d\n", iter_count ++);
                if(i == mic_num){
                    //cpu_correct(cpu_buffer, compute_data_size, &out_data_size, read_offset, read_offset_size, arg, bit_vector, statistics);
                    mic_correct(cpu_buffer, compute_data_size, &out_data_size, read_offset, read_offset_size, arg, bit_vector, hash_seed, statistics);
                    //printf("compute done\n");
                }
                else{

#ifdef __MIC__
                    #pragma offload target(mic:i)\
                        nocopy(mmic_buffer:length(mic_buffer_size) alloc_if(0) free_if(0))\
                        in(read_offset:length(0) alloc_if(0) free_if(0))\
                        in(bit_vector:length(0) alloc_if(0) free_if(0))\
                        in(hash_seed:length(0) alloc_if(0) free_if(0))\
                        in(hash_array, notrim)\
                        in(compute_data_size, read_offset_size)\
                        inout(out_data_size, statistics)
                        {
                            init_args arg;
                            arg.notrim = notrim;
                            array_2_args(hash_array, arg);
                            mic_correct(mmic_buffer,compute_data_size, &out_data_size, read_offset, read_offset_size, arg, bit_vector, hash_seed, statistics);
                        }
#endif
                }
                GETTIME(e_time);
                compute_time[i] += DIFFTIME(s_time, e_time);

                GETTIME(s_time);
                
                pthread_mutex_lock(&malloc_lock);
                {
                   dwrite_buffer = (char*)malloc(device_buffer_size[i] * sizeof(char));
                } 
                pthread_mutex_unlock(&malloc_lock);
                if(i ==  mic_num){     
                    memcpy(dwrite_buffer, cpu_buffer, compute_data_size);
                }else{

#ifdef __MIC__
    #pragma offload target(mic:i)\
                    out(mmic_buffer[0:mic_buffer_size]:into(dwrite_buffer[0 :mic_buffer_size])  alloc_if(0) free_if(0))
                    {
                    }
#endif
                }

                omp_set_lock(&write_lock);
                trimmed_data_size[write_iteration] = out_data_size;
                computed_data_size[write_iteration] = compute_data_size;
                wbuffer_list[write_iteration] = dwrite_buffer;
                write_flag[write_iteration ++] = i;
                omp_unset_lock(&write_lock);

                sem_post(&write_sem_count);
                
                GETTIME(e_time);
                offload_out_time[i] += DIFFTIME(s_time, e_time);
            }
            #pragma omp critical
            {
                merge_statistics(statistics_result, statistics);
            }
            //printf("%d thread %d done\n", my_rank, i);
        } if(i == device_num){

            process_data_to_compute(read_sem_count, rbuffer_count_sem, compute_sem_count, device_compute_buffer, device_buffer_size, device_data_size, read_flag, rbuffer_list, rbuffer_data_size, pfile_len, mic_num, device_num, read_time, read_total_time);

            //printf("%d thread %d done\n", my_rank, i);
        }else if(i == device_num + 1){

            process_computed_data_to_buffer(write_sem_count, trimmed_data_size, write_flag, wbuffer_list, computed_data_size, pfile_len, write_time, write_total_time, f_corrected_read);

            //printf("%d thread %d done\n", my_rank, i);
        }else if(i == device_num + 2){ 
             
            uint32_t iter_count = 0;
            uint64_t total_readed_size = rbuffer_data_size[0];
            if(total_readed_size < pfile_len){
                while(1){
                    iter_count += 1;
                    pthread_mutex_lock(&malloc_lock);
                    {
                        rbuffer_list[iter_count] = (char*)malloc(rbuffer_size * sizeof(char));
                    }
                    pthread_mutex_unlock(&malloc_lock);
                
                    rbuffer_data_size[iter_count] = fread(rbuffer_list[iter_count], sizeof(char), rbuffer_size,  pfile);
                //printf("read data %d\n", my_rank);
                    total_readed_size += rbuffer_data_size[iter_count];
                    sem_post(&rbuffer_count_sem);
                    if(total_readed_size == pfile_len)
                        break;
                }
                //printf("%d thread %d done\n", my_rank, i);
                //printf("%d thread %d done total readed size %lld\n", my_rank, i, total_readed_size);
            }
        }

    }
    omp_destroy_lock(&read_lock);
    omp_destroy_lock(&write_lock);
    sem_destroy(&read_sem_count);
    sem_destroy(&write_sem_count);
    sem_destroy(&rbuffer_count_sem);
    pthread_mutex_destroy(&malloc_lock);
    for(int i = 0; i < device_num; i ++){
        sem_destroy(&compute_sem_count[i]);
    }

    GETTIME(end_time);
    //time(&rawtime);
    //c_inst_time.end_correct_errors_in_reads = asctime(localtime(&rawtime));
    //total_time = DIFFTIME(start_time, end_time);
    //printf("total correced errors %ld\ntotal trimmed bases %ld\ntotal corrected reads %ld\n", statistics_result.num_corrected_errors, statistics_result.num_trimmed_bases, statistics_result.num_corrected_reads);
    //printf("total time is %lf\n", total_time);
    //printf("read total time is %lf\n", read_total_time);
    //printf("write total time is %lf\n", write_total_time);
    //printf("read time is ");
    //for(int i = 0; i < device_num; i ++)
    //    printf("%lf ", read_time[i]);
    //printf("\n");
    //printf("offload time is ");
    //for(int i = 0; i < device_num; i ++)
    //    printf("%lf ", offload_time[i]);
    //printf("\n");
    //printf("compute time is ");
    //for(int i = 0; i < device_num; i ++)
    //    printf("%lf ", compute_time[i]);
    //printf("\n");
    //printf("offload out  time is ");
    //for(int i = 0; i < device_num; i ++)
    //    printf("%lf ", offload_out_time[i]);
    //printf("\n");
    //printf("write time is ");
    //for(int i = 0; i < device_num; i ++)
    //    printf("%lf ", write_time[i]);
    //printf("\n");

    for(int i = 0; i < mic_num; i ++){
#ifdef __MIC__
#pragma offload target(mic:i)\
        in(read_offset:length(read_offset_size) alloc_if(0) free_if(1))\
        in(bit_vector:length(arg.bit_vector_width_byte) alloc_if(0) free_if(1))\
        in(hash_seed:length(arg.num_hash_func) alloc_if(0) free_if(1))
        {
        }
#pragma offload target(mic:i)\
        in(mmic_buffer:length(mic_buffer_size) alloc_if(0) free_if(1))
        {
        }
#endif
        free(device_compute_buffer[i]);
    }


    free(cpu_buffer);
    free(read_offset);
    free(mmic_buffer);
    free(device_compute_buffer[mic_num]);
    mpi_total_time += total_time;


}

