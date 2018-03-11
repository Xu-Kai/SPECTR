/*
 * correct_errors.c
 *
 *  Created on: 2015年12月14日
 *      Author: xk
 */


#include "correct_errors.h"
//#include "query_text.h"
//#include "cpu_vec_query_text.h"
#include "new_query_text.h"
#define BLOCK_READ_NUM 2000000
#define READ_MAX_LEN 256
//#define READ_LEN 101
#include "struct.h"
#include "read_inputs.h"

 void cpu_extend_a_kmer(const char* kmer, const char* sequence, uint32_t index_kmer, const uint32_t& index_last_mod, 
        uint32_t* pos_path_vec, char* c_path_vec, uint32_t& c_path_size, const uint32_t org_boundary_left, const uint32_t  org_boundary_right, 
        const char* quality_score, const unsigned char* bit_vector, const uint32_t* hash_seed, bool& run_exploration, init_args& arg, func_mem& funcm){
    
	char *path_content = funcm.cbuffer + (MAX_CANDIDATE_PATHS + 1 + 3) * READ_MAX_LEN;
    int index_pos = 0;
   	char *kmer_ext = funcm.cbuffer + (MAX_CANDIDATE_PATHS + 1 + 7) * READ_MAX_LEN;
    char *kmer_new = funcm.cbuffer + (MAX_CANDIDATE_PATHS + 1 + 8) * READ_MAX_LEN;
    uint32_t *rec_ind = funcm.ubuffer + (MAX_CANDIDATE_PATHS + 1 + 4) * READ_MAX_LEN;
    uint32_t flag = 1;
    uint32_t index_start = index_kmer;
    memcpy(kmer_ext, kmer + 1, arg.kmer_length - 1);
    while(index_pos > 0 || flag == 1){
        flag = 0;
        //printf("index_pos %d \n", index_pos);
        memcpy(kmer_new, kmer_ext + index_pos, arg.kmer_length - 1);
        //kmer_new[arg.kmer_length - 1] = sequence[index_kmer + arg.kmer_length -1];
        //path_query_text(kmer_new, bit_vector, hash_seed, arg.kmer_length, arg.num_hash_func, arg.bit_vector_width, 1);
        uint32_t k_ind = 0;
        if ((index_kmer == (org_boundary_left - 1)) || (index_kmer == (org_boundary_right - 1)) || 
                (((uint32_t)quality_score[index_kmer + arg.kmer_length] - arg.quality_score_offset) <= arg.extremely_low_quality_score)){

            for(unsigned short int it_alter = A; it_alter <= T; it_alter ++){ 
                kmer_new[arg.kmer_length - 1] = NEOCLEOTIDE[it_alter];
                //if(kmer_new[funcm.kmer_length + it_alter]){
                if(query_text(kmer_new, bit_vector, &hash_seed[0], arg.kmer_length, arg.num_hash_func, arg.bit_vector_width) == true){
                    path_content[(index_pos << 2) + k_ind ++] = NEOCLEOTIDE[it_alter];      
                }
            }                
        }else{
            //printf("index_kmer %d %d \n", index_kmer, arg.kmer_length);
            kmer_new[arg.kmer_length - 1] = sequence[index_kmer + arg.kmer_length];
            if(query_text(kmer_new, bit_vector, &hash_seed[0], arg.kmer_length, arg.num_hash_func, arg.bit_vector_width) == true){ 
                path_content[(index_pos << 2) + k_ind ++] = kmer_new[arg.kmer_length - 1];      
            }else{
                 for(unsigned short int it_alter = A; it_alter <= T; it_alter ++){ 
                    kmer_new[arg.kmer_length - 1] = NEOCLEOTIDE[it_alter];
                    //if(kmer_new[funcm.kmer_length + it_alter]){
                    if(sequence[index_kmer + arg.kmer_length] != NEOCLEOTIDE[it_alter])
                    if(query_text(kmer_new, bit_vector, &hash_seed[0], arg.kmer_length, arg.num_hash_func, arg.bit_vector_width) == true){
                        path_content[(index_pos << 2) + k_ind ++] = NEOCLEOTIDE[it_alter];      
                    }
                }      
            }
        
        }

        //printf("index %d %d \n",k_ind, index_kmer);
        if(k_ind < 4){
            path_content[(index_pos << 2) + k_ind] = 0;
        }
        if((index_kmer + 1 == index_last_mod) && (k_ind > 0)){
            for(uint32_t it = 0; it < k_ind; it ++){
                kmer_ext[index_pos + arg.kmer_length - 1] = path_content[(index_pos << 2) + it];
            	pos_path_vec[c_path_size * READ_MAX_LEN] = 0;
                for(uint32_t ct = arg.kmer_length - 2; ct < index_pos + arg.kmer_length; ct ++){
                    if(kmer_ext[ct] != sequence[index_start + ct + 1]){
                        pos_path_vec[c_path_size * READ_MAX_LEN + pos_path_vec[c_path_size * READ_MAX_LEN] + 1] = index_start + ct + 1;
                        c_path_vec[c_path_size * READ_MAX_LEN + pos_path_vec[c_path_size * READ_MAX_LEN] + 1] = kmer_ext[ct];
                        pos_path_vec[c_path_size * READ_MAX_LEN] += 1;
                    }
                }
                c_path_size += 1;        
                //printf("%d\n", c_path_size);
                if(c_path_size > MAX_CANDIDATE_PATHS){
                    //退出寻找路径
                    index_pos = -1;        
                    run_exploration = false;
                    break;
                }else{
                    if(kmer_ext[arg.kmer_length - 2] != sequence[index_start + arg.kmer_length - 1]){
                        pos_path_vec[c_path_size * READ_MAX_LEN + 1] = index_start + arg.kmer_length - 1;
                        c_path_vec[c_path_size * READ_MAX_LEN + 1] = kmer_ext[arg.kmer_length - 2];
                        pos_path_vec[c_path_size * READ_MAX_LEN] = 1;
                    }else{
                        pos_path_vec[c_path_size * READ_MAX_LEN] = 0; 
                    }       
                }  

            }
            
            while(index_pos - 1 >= 0){
                    if(rec_ind[index_pos - 1] < 4 && (path_content[((index_pos - 1) << 2 ) + rec_ind[index_pos - 1]] > 0)){
                    kmer_ext[arg.kmer_length - 1 + index_pos - 1] = path_content[((index_pos - 1)<< 2) + rec_ind[index_pos - 1]];
                    rec_ind[index_pos - 1] += 1;
                    break;
                }else{
                    index_kmer -= 1;
                    index_pos -= 1;
                }
            }
            /*
            index_pos -= 1;
            index_kmer -= 1;
            while(index_pos >= 0){
                if(rec_ind[index_pos] < 4 && (path_content[(index_pos << 2 ) + rec_ind[index_pos]] > 0)) {
                    kmer_ext[arg.kmer_length - 1 + index_pos] = path_content[(index_pos << 2) + rec_ind[index_pos]];
                    rec_ind[index_pos] += 1;
                    break;
                }else{
                    index_kmer -= 1;
                    index_pos -=1;
                }
            }
            */
        }
        else if (k_ind > 0){
            index_kmer += 1;
            rec_ind[index_pos] = 0;
            kmer_ext[arg.kmer_length  + index_pos - 1] = path_content[(index_pos << 2) + rec_ind[index_pos]];
            rec_ind[index_pos] += 1;
            index_pos += 1;             
            //printf("p 1 \n"); 
        }else{
            //printf("%d %d %x\n", index_pos, rec_ind[index_pos], path_content[(index_pos << 2) + rec_ind[index_pos]]);
            while(index_pos - 1 >= 0){
                    if(rec_ind[index_pos - 1] < 4 && (path_content[((index_pos - 1) << 2 ) + rec_ind[index_pos - 1]] > 0)){
                    kmer_ext[arg.kmer_length - 1 + index_pos - 1] = path_content[((index_pos - 1)<< 2) + rec_ind[index_pos - 1]];
                    rec_ind[index_pos - 1] += 1;
                    break;
                }else{
                    index_kmer -= 1;
                    index_pos -= 1;
                }
            }
            //printf("p 2 \n");
        }
        /*
        for(int ii = 0; ii < index_pos; ii ++){
            for(uint32_t jj = 0; jj < 4; jj ++){
                if(path_content[(ii << 2) + jj] > 0){
                    printf("%c", path_content[(ii << 2) + jj]);
                }else{
                    break;
                }
            }
            printf(" ");
        }
        printf("\n");
        for(uint32_t ii = 0; ii < index_pos + arg.kmer_length - 1; ii ++){
            printf("%c", kmer_ext[ii]);
        }
        printf("\n");
        for(int ii = 0; ii < index_pos; ii ++){
            printf("%d ", rec_ind[ii]);
        }
        printf("*****\n");
        char cc, ccc;
        scanf("%c%c", &cc, &ccc);
        printf("cc %c \n", cc); 
        */
    }

}
 void cpu_modify_sequence(const char* org_sequence, char* sequence_modified, const char* quality_score, char* sequence_modification, 
        uint32_t& trim_5_end, uint32_t& trim_3_end, const uint32_t read_length, uint32_t& num_corrected_errors_local, uint32_t* pos_path_vec, 
        char* c_path_vec, uint32_t c_path_size, init_args& arg, func_mem& funcm, uint32_t it_path_1st, uint32_t qs_1st, uint32_t qs_2nd, 
        bool too_many_correction, uint32_t max_trimmed_bases, bool trim_flag){
     
        if(too_many_correction == false && qs_1st + MIN_QS_DIFF <= qs_2nd){
            for(uint32_t it_base = 0; it_base < pos_path_vec[it_path_1st * READ_MAX_LEN]; it_base ++){
                sequence_modification[pos_path_vec[it_path_1st * READ_MAX_LEN + it_base + 1]] = c_path_vec[it_path_1st * READ_MAX_LEN + it_base + 1];
                sequence_modified[pos_path_vec[it_path_1st * READ_MAX_LEN + it_base + 1]] = c_path_vec[it_path_1st * READ_MAX_LEN + it_base + 1];
                num_corrected_errors_local ++;
            }
        }else{
            /*
            printf("KKKKKKKKK\n");
            for(uint32_t ii = 0; ii < c_path_size; ii ++){
                for(uint32_t jj = 0; jj < pos_path_vec[ii * READ_MAX_LEN]; jj ++){
                    printf("%d %c | ", pos_path_vec[ii * READ_MAX_LEN + jj + 1], c_path_vec[ii * READ_MAX_LEN + jj + 1]);
                }
                printf("****\n");
            }
            */
            //uint32_t pos_base_vec_intersection[READ_MAX_LEN];
            //char c_base_vec_intersection[READ_MAX_LEN];
			
			uint32_t *pos_base_vec_intersection = funcm.ubuffer +  (MAX_CANDIDATE_PATHS + 1 + 4) * READ_MAX_LEN;
			char *c_base_vec_intersection = funcm.cbuffer + (MAX_CANDIDATE_PATHS + 1 + 2) * READ_MAX_LEN;
			
			
            //uint32_t pos_base_vec_intersection_prev[READ_MAX_LEN];
            //char c_base_vec_intersection_prev[READ_MAX_LEN];
			uint32_t *pos_base_vec_intersection_prev = funcm.ubuffer + (MAX_CANDIDATE_PATHS + 1 + 5) * READ_MAX_LEN;
			char *c_base_vec_intersection_prev = funcm.cbuffer + (MAX_CANDIDATE_PATHS + 1 + 3) * READ_MAX_LEN;
			
            //uint32_t pos_base_vec_union[75*READ_MAX_LEN];
            //char c_base_vec_union[75*READ_MAX_LEN];
            uint32_t *pos_base_vec_union = funcm.ubuffer + (MAX_CANDIDATE_PATHS + 1 + 6) * READ_MAX_LEN;
			char *c_base_vec_union = funcm.cbuffer + (MAX_CANDIDATE_PATHS + 1 + 4) * READ_MAX_LEN;
			

            //uint32_t pos_base_vec_union_prev[75*READ_MAX_LEN];
            //char c_base_vec_union_prev[75*READ_MAX_LEN];
			uint32_t *pos_base_vec_union_prev = funcm.ubuffer + (MAX_CANDIDATE_PATHS * 2 + 2 + 6) * READ_MAX_LEN;
			char *c_base_vec_union_prev = funcm.cbuffer + (MAX_CANDIDATE_PATHS * 2 + 2  + 4) * READ_MAX_LEN;
            memcpy(pos_base_vec_intersection_prev, pos_path_vec, (pos_path_vec[0] + 1) * sizeof(uint32_t));
            memcpy(c_base_vec_intersection_prev, c_path_vec, pos_path_vec[0] + 1);
 
            memcpy(pos_base_vec_union_prev, pos_path_vec, (pos_path_vec[0] + 1) * sizeof(uint32_t));
            memcpy(c_base_vec_union_prev, c_path_vec, pos_path_vec[0] + 1);
            //printf("c_path_size == %d \n", c_path_size);
            for(uint32_t it_p = 1; it_p < c_path_size; it_p ++){
                uint32_t in1_begin = 0, in2_begin = 0, in1_end = pos_base_vec_intersection_prev[0], in2_end = pos_path_vec[it_p * READ_MAX_LEN];  
                uint32_t offset = it_p * READ_MAX_LEN;
                pos_base_vec_intersection[0] = 0;
                pos_base_vec_union[0] = 0;
                c_base_vec_union[0] = 0;
                //printf("ppppp\n");
                while(in1_begin != in1_end && in2_begin != in2_end){
                    if(pos_base_vec_intersection_prev[in1_begin + 1] < pos_path_vec[offset + in2_begin + 1]){
                        in1_begin += 1;
                    }else if(pos_base_vec_intersection_prev[in1_begin + 1] > pos_path_vec[offset + in2_begin + 1]){
                        in2_begin += 1;
                    }else{
                        if(c_base_vec_intersection_prev[in1_begin + 1] ==  c_path_vec[offset + in2_begin + 1]){
                            pos_base_vec_intersection[pos_base_vec_intersection[0] + 1] = pos_base_vec_intersection_prev[in1_begin +1];
                            c_base_vec_intersection[pos_base_vec_intersection[0] + 1] = c_base_vec_intersection_prev[in1_begin +1];
                            pos_base_vec_intersection[0] += 1;
                        }
                        in1_begin += 1;
                        in2_begin += 1;
                    }
                } 

                //printf("HHHHHHHHH\n");
                in1_begin = 0, in2_begin = 0, in1_end = pos_base_vec_union_prev[0], in2_end = pos_path_vec[it_p * READ_MAX_LEN];  
                while(1){
                    //printf("in1 in2 %d %d \n", in1_begin, in2_begin);
                    if(in1_begin == in1_end){
                        memcpy(pos_base_vec_union + pos_base_vec_union[0] + 1, pos_path_vec + offset + in2_begin + 1, 
                                (pos_path_vec[it_p * READ_MAX_LEN] - in2_begin)* sizeof(uint32_t));
                        memcpy(c_base_vec_union + pos_base_vec_union[0] + 1, c_path_vec + offset + in2_begin + 1, 
                                (pos_path_vec[it_p * READ_MAX_LEN] - in2_begin));
                        pos_base_vec_union[0] += pos_path_vec[offset] - in2_begin;

                    //printf("uion1 %d \n", pos_base_vec_union[0]);
                        break;
                    }
                    else if(in2_begin == in2_end){ 
                        memcpy(pos_base_vec_union + pos_base_vec_union[0] + 1, pos_base_vec_union_prev +  in1_begin + 1, 
                                (pos_base_vec_union_prev[0] - in1_begin)* sizeof(uint32_t));
                        memcpy(c_base_vec_union + pos_base_vec_union[0] + 1, c_base_vec_union_prev +  in1_begin + 1, 
                                (pos_base_vec_union_prev[0] - in1_begin));
                        pos_base_vec_union[0] += pos_base_vec_union_prev[0] - in1_begin;
                    //printf("uion2 %d \n", pos_base_vec_union[0]);
                        break;
                    } else {
                        if(pos_base_vec_union_prev[in1_begin + 1] < pos_path_vec[offset + in2_begin + 1]){
                            pos_base_vec_union[pos_base_vec_union[0] + 1] = pos_base_vec_union_prev[in1_begin + 1];
                            c_base_vec_union[pos_base_vec_union[0] + 1] = c_base_vec_union_prev[in1_begin + 1];
                            pos_base_vec_union[0] += 1;
                            in1_begin += 1;
                        } else if(pos_base_vec_union_prev[in1_begin + 1] > pos_path_vec[offset + in2_begin + 1]){ 
                            pos_base_vec_union[pos_base_vec_union[0] + 1] = pos_path_vec[offset + in2_begin + 1];
                            c_base_vec_union[pos_base_vec_union[0] + 1] = c_path_vec[offset + in2_begin + 1];
                            pos_base_vec_union[0] += 1;
                            in2_begin += 1;
                        } else{ 
                            pos_base_vec_union[pos_base_vec_union[0] + 1] = pos_base_vec_union_prev[in1_begin + 1];
                            c_base_vec_union[pos_base_vec_union[0] + 1] = c_base_vec_union_prev[in1_begin + 1];
                            pos_base_vec_union[0] += 1;
                            in1_begin += 1;
                            in2_begin += 1;
                        }
                    }
                    //printf("uion %d \n", pos_base_vec_union[0]);
                }
                //printf("MMMMMMMMMMM\n");
                //printf("pos_base_vec_intersection == %d %d \n", pos_base_vec_intersection[0], pos_base_vec_union[0]);
                memcpy(pos_base_vec_intersection_prev, pos_base_vec_intersection, (pos_base_vec_intersection[0] + 1) * sizeof(uint32_t));
                memcpy(c_base_vec_intersection_prev, c_base_vec_intersection, pos_base_vec_intersection[0] + 1);
 
                memcpy(pos_base_vec_union_prev, pos_base_vec_union, (pos_base_vec_union[0] + 1) * sizeof(uint32_t));
                memcpy(c_base_vec_union_prev, c_base_vec_union, pos_base_vec_union[0] + 1);
                //printf("nnnnnnnnnn\n");
            }
           	/* 
			for(uint32_t ii = 0; ii < pos_base_vec_union[0]; ii ++){
				printf("%d %c\n", pos_base_vec_union[ii + 1], c_base_vec_union[ii + 1]);
			}
			printf("union\n");
			for(uint32_t ii = 0; ii < pos_base_vec_intersection[0]; ii ++){
				printf("%d %c\n", pos_base_vec_intersection[ii + 1], c_base_vec_intersection[ii + 1]);
			}
			printf("intersection\n");
			*/
            uint32_t *pos_base_vec_difference = pos_base_vec_union_prev;
            char *c_base_vec_difference = c_base_vec_union_prev;
            pos_base_vec_difference[0] = 0;
            uint32_t in1_begin = 0, in2_begin = 0, in1_end = pos_base_vec_union[0], in2_end = pos_base_vec_intersection[0];
            while(in1_begin != in1_end && in2_begin != in2_end){
                if(pos_base_vec_union[in1_begin + 1] < pos_base_vec_intersection[in2_begin + 1]){
                    pos_base_vec_difference[pos_base_vec_difference[0] + 1] = pos_base_vec_union[in1_begin + 1];
                    c_base_vec_difference[pos_base_vec_difference[0] + 1] = c_base_vec_union[in1_begin + 1];
                    pos_base_vec_difference[0] += 1;
                    in1_begin += 1;
                }else if(pos_base_vec_intersection[in2_begin + 1] < pos_base_vec_union[in1_begin + 1]){
                    in2_begin += 1;
                }else{
                    in1_begin += 1;
                    in2_begin += 1;
                }
            }
            //printf("in1 in2 %d %d %d %d\n", in1_begin, in1_end, in2_begin, in2_end);
            memcpy(pos_base_vec_difference + pos_base_vec_difference[0] + 1, pos_base_vec_union + in1_begin + 1, (in1_end - in1_begin) * sizeof(uint32_t));
            memcpy(c_base_vec_difference + pos_base_vec_difference[0] + 1, c_base_vec_union + in1_begin + 1, in1_end - in1_begin);
            pos_base_vec_difference[0] += in1_end - in1_begin;
			/*	
			for(uint32_t ii = 0; ii < pos_base_vec_difference[0]; ii ++){
				printf("%d %c\n", pos_base_vec_difference[ii + 1], c_base_vec_difference[ii + 1]);
			}
			printf("difference\n");
			*/
            
            //printf("%d\n", pos_base_vec_difference[0]);
			/*
            for(uint32_t ii = 0; ii < pos_base_vec_intersection[0]; ii ++){
                printf("%d %c\n", pos_base_vec_intersection[ii + 1], c_base_vec_intersection[ii + 1]);
            }
            printf("hhhhhh\n");

            for(uint32_t ii = 0; ii < pos_base_vec_union[0]; ii ++){
                printf("%d %c\n", pos_base_vec_union[ii + 1], c_base_vec_union[ii + 1]);
            }
            printf("hhhhhh\n");
            for(uint32_t ii = 0; ii < pos_base_vec_difference[0]; ii ++){
                printf("%d %c\n", pos_base_vec_difference[ii + 1], c_base_vec_difference[ii + 1]);
            }
            printf("hhhhhh\n");
            */
            if(pos_base_vec_difference[0] > 0){
              uint32_t vector_index_leftmost = 0;
              uint32_t vector_index_rightmost = pos_base_vec_difference[0] - 1;
              bool keep_going = true;
              while(keep_going){
                if(pos_base_vec_difference[vector_index_leftmost + 1] + 1 < read_length - pos_base_vec_difference[vector_index_rightmost + 1]){
                    if(pos_base_vec_difference[vector_index_leftmost + 1] + 1 + trim_3_end <= max_trimmed_bases){
                        if(pos_base_vec_difference[vector_index_leftmost + 1] + 1 > trim_5_end || trim_flag){
                            trim_5_end =pos_base_vec_difference[vector_index_leftmost + 1] + 1;
                        }

                        if(vector_index_leftmost == vector_index_rightmost){
                            keep_going = false;
                        } else{
                            vector_index_leftmost += 1;
                        }
                    } 
                    else{
                        keep_going = false;
                    }
                } else{
                    if(read_length - pos_base_vec_difference[vector_index_rightmost + 1] <= max_trimmed_bases){
                        if(read_length - pos_base_vec_difference[vector_index_rightmost + 1] > trim_3_end || trim_flag){
                            trim_3_end = read_length - pos_base_vec_difference[vector_index_rightmost + 1];
                        }
                        if(vector_index_leftmost == vector_index_rightmost){
                            keep_going = false;
                        }else{
                            vector_index_rightmost -= 1;
                        }
                    }else{
                        keep_going = false;
                    }
                }
              } 
            }
			
			//printf("tirm_5_end == %d\n", trim_5_end);
			//printf("trim_3_end == %d\n", trim_3_end);
			
           /* 
            for(uint32_t ii = 0; ii < pos_base_vec_difference[0]; ii ++){
                printf("%d %c\n", pos_base_vec_difference[ii + 1], c_base_vec_difference[ii + 1]);
            }
            printf("\n");
            printf("trim 3 5 %d %d %d\n", trim_3_end, trim_5_end, read_length);
            */
            for(uint32_t it_inter = 0; it_inter < pos_base_vec_intersection[0]; it_inter ++){
                if(pos_base_vec_intersection[it_inter + 1] < read_length - trim_3_end && pos_base_vec_intersection[it_inter + 1] >= trim_5_end){
                    if(sequence_modified[pos_base_vec_intersection[it_inter + 1]] != c_base_vec_intersection[it_inter + 1]){
                        sequence_modification[pos_base_vec_intersection[it_inter + 1]] = c_base_vec_intersection[it_inter + 1];
                        sequence_modified[pos_base_vec_intersection[it_inter + 1]] = c_base_vec_intersection[it_inter + 1];
                        num_corrected_errors_local += 1; 
                    }
                }           
            }
            
        } 
 
}
    /*
    for(uint32_t ii = 0; ii < read_length; ii ++){
        printf("%c", sequence_modified[ii]);
    } 
    printf("\n");
    */
/*

        if(too_many_correction == false && qs_1st + MIN_QS_DIFF <= qs_2nd){
            for(uint32_t it_base = 0; it_base < pos_path_vec[it_path_1st * READ_MAX_LEN]; it_base ++){
                sequence_modification[pos_path_vec[it_path_1st * READ_MAX_LEN + it_base + 1]] = c_path_vec[it_path_1st * READ_MAX_LEN + it_base + 1];
                sequence[pos_path_vec[it_path_1st * READ_MAX_LEN + it_base + 1]] = c_path_vec[it_path_1st * READ_MAX_LEN + it_base + 1];
                num_corrected_errors_local ++;
            }
        }else{
            uint32_t pos_base_vec_intersection[READ_MAX_LEN];
            char c_base_vec_intersection[READ_MAX_LEN];

            uint32_t pos_base_vec_union[READ_MAX_LEN];
            char c_base_vec_unio[READ_MAX_LEN];
            
            uint32_t pos_base_vec_intersection_prev[READ_MAX_LEN];
            char c_base_vec_intersection_prev[READ_MAX_LEN];

            uint32_t pos_base_vec_union_prev[READ_MAX_LEN];
            char c_base_vec_union_prev[READ_MAX_LEN];

            memcpy(pos_base_vec_intersection_prev, pos_path_vec, pos_path_vec[0] + 1);
            memcpy(c_base_vec_intersection_prev, c_path_vec, pos_path_vec[0] + 1);
 
            memcpy(pos_base_vec_union_prev, pos_path_vec, pos_path_vec[0] + 1);
            memcpy(c_base_vec_union_prev, c_path_vec, pos_path_vec[0] + 1);
            for(uint32_t it_p = 1; it_p < c_path_size; it_p ++){
                pos_base_vec_intersection[0] = 0;
                pos_base_vec_union[0] = 0;
                uint32_t in1_begin = 0, in2_begin = 0, in1_end = pos_base_vec_intersection_prev[0], in2_end = pos_path_vec[it_p * READ_MAX_LEN];  
                uint32_t offset = it_p * READ_MAX_LEN;
                while(in1_begin != in1_end && in2_begin != in2_end){
                    if(pos_base_vec_intersection_prev[in1_begin + 1] < pos_path_vec[offset + in2_begin + 1]){
                        in1_begin += 1;
                    }else if(pos_base_vec_intersection_prev[in1_begin + 1] > pos_path_vec[offset + in2_begin + 1]){
                        in2_begin += 1;
                    }else{
                        if(c_base_vec_intersection_prev[in1_begin + 1] ==  c_path_vec[offset + in2_begin + 1]){
                            pos_base_vec_intersection[pos_base_vec_intersection[0] + 1] = c_base_vec_intersection_prev[in1_begin +1];
                            c_base_vec_intersection[pos_base_vec_intersection[0] + 1] = c_base_vec_intersection_prev[in1_begin +1];
                            pos_base_vec_intersection[0] += 1;
                        }
                        in1_begin += 1;
                        in2_begin += 1;
                    }
                } 
                in1_begin = 0, in2_begin = 0, in1_end = pos_base_vec_union_prev[0], in2_end = pos_path_vec[it_p * READ_MAX_LEN];  
                while(1){
                    if(in1_begin == in1_end){
                        memcpy(pos_base_vec_union + pos_base_vec_union[0] + 1, pos_path_vec + offset + in2_begin + 1, pos_path_vec[it_p * READ_MAX_LEN] - in2_begin);
                        pos_base_vec_union[0] = pos_path_vec[offset] - in2_begin;
                        break;
                    }
                    else if(in2_begin == in2_end){ 
                        memcpy(pos_base_vec_union + pos_base_vec_union[0] + 1, pos_base_vec_union_prev +  in1_begin + 1, pos_base_vec_union_prev[0] - in1_begin);
                        pos_base_vec_union[0] = pos_path_vec[offset] - in2_begin;
                        break;
                    } else {
                        if(pos_base_vec_union_prev[in1_begin + 1] < pos_path_vec[offset + in2_begin + 1]){
                            pos_base_vec_union[pos_base_vec_union[0] + 1] = pos_base_vec_union_prev[in1_begin + 1];
                            c_base_vec_union[pos_base_union[0] + 1] = pos_base_vec_union_prev[in1_begin + 1];
                            pos_base_union[0] += 1;
                            in1_begin += 1;
                        } else if(pos_base_vec_union_prev[in1_begin + 1] > pos_path_vec[offset + in2_begin + 1]){ 
                            pos_base_vec_union[pos_base_vec_union[0] + 1] = pos_path_vec[offset + in2_begin + 1];
                            c_base_vec_union[pos_base_union[0] + 1] = c_path_vec[offset + in2_begin + 1];
                            pos_base_union[0] += 1;
                            in2_begin += 1;
                        } else{ 
                            pos_base_vec_union[pos_base_vec_union[0] + 1] = pos_base_vec_union_prev[in1_begin + 1];
                            c_base_vec_union[pos_base_union[0] + 1] = pos_base_vec_union_prev[in1_begin + 1];
                            pos_base_vec_union[0] += 1;
                            in1_begin += 1;
                            in2_begin += 1;
                        }
                    }
                }
                memcpy(pos_base_vec_intersection_prev, pos_base_vec_intersection, pos_base_vec_intersection[0] + 1);
                memcpy(c_base_vec_intersection_prev, c_base_vec_intersection, pos_base_vec_intersection[0] + 1);
 
                memcpy(pos_base_vec_union_prev, pos_base_vec_union, pos_base_vec_union[0] + 1);
                memcpy(c_base_vec_union_prev, c_base_vec_union, pos_base_vec_union[0] + 1);
            }
            
            uint32_t *pos_base_difference = pos_base_vec_intersection_prev;
            char *c_base_difference = c_base_vec_intersection_prev;
            pos_base_vec_difference[0] = 0;
            uint32_t in1_begin = 0, in2_begin = 0, in1_end = pos_base_vec_union[0], in2_end = pos_base_vec_intersection[0];
            while(in1_begin != in1_end && in2_begin != in2_end){
                if(pos_base_vec_union[in1_begin + 1] < pos_base_vec_intersection[in2_begin + 1]){
                    pos_base_difference[pos_base_vec_difference[0] + 1] = pos_base_vec_union[in1_begin + 1];
                    c_base_difference[pos_base_vec_difference[0] + 1] = c_base_vec_union[in1_begin + 1];
                    pos_base_vec_difference[0] += 1;
                    in1_begin += 1;
                }else if(in2_begin < in1_begin){
                    in2_begin += 1;
                }else{
                    in1_begin += 1;
                    in2_begin += 1;
                }
            }
            memcpy(pos_base_vec_difference, pos_base_union + in1_begin, in1_end - in1_begin);
            memcpy(c_base_vec_difference, c_base_union + in1_begin, in1_end - in1_begin);
            if(pos_base_vec_difference[0])%{
               uint32_t vector_index_leftmost = 0;
                uint32_t vector_index_rightmost = pos_base_vec_difference[0] - 1;
               bool keep_going = true;
              while(keep_going){
                if(pos_base_vec_difference[vector_index_leftmost + 1] + 1 < read_length - pos_base_vec_difference[vec_index_rightmost + 1]){
                    if(pos_base_vec_difference[vector_index_leftmost + 1] + 1 + trim_3_end <= max_trimmed_bases){
                        if(pos_base_vec_difference[vector_index_leftmost + 1] + 1 > trim_5_end){
                            trim_5_end =pos_base_vec_difference[vector_index_leftmost + 1] + 1;
                        }

                        if(vector_index_leftmost == vector_index_rightmost){
                            keep_going = false;
                        } else{
                            vector_index_leftmost += 1;
                        }
                    } 
                    else{
                        keep_going = false;
                    }
                } else{
                    if(read_length - pos_base_vec_difference[vector_index_rightmost + 1] <= max_trimmed_bases){
                        if(read_length - pos_base_vec_difference[vector_index_rightmost] > trim_3_end){
                            trim_3_end = read_length - pos_base_vec_difference[vector_index_rightmost + 1];
                        }
                        if(vector_index_leftmost == vector_index_rightmost){
                            keep_going = false;
                        }else{
                            vector_index_rightmost -= 1;
                        }
                    }else{
                        keep_going = false;
                    }
                }
              } 
            }
            for(uint32_t it_inter = 0; it_inter < pos_base_vec_difference[0]; it_inter ++){
                if(pos_base_difference[it_inter + 1] < read_length - trim_3_end && pos_base_vec_difference[it_inter + 1] >= trim_5_end){
                    if(sequence[pos_base_vec_difference[it_inter + 1]] != c_base_vec_difference[pos_base_vec_difference[it_inter + 1]]){
                        sequence_modification[pos_base_vec_difference[it_inter + 1]] = c_base_vec_difference[it_inter + 1];
                        sequence[pos_base_vec_difference[it_inter + 1]] = c_base_vec_difference[it_inter + 1];
                        num_corrected_errors_local += 1; 
                    }
                }           
            }
       
        }
        */

 void cpu_correct_errors_between_solid_regions(const char* org_sequence, char* sequence_modified, const char* quality_score, const uint32_t left_first, 
        const uint32_t index_start, const uint32_t index_end, const uint32_t right_second, const uint32_t org_boundary_left, const uint32_t org_boundary_right, 
        char* sequence_modification, uint32_t& trim_5_end, uint32_t& trim_3_end, uint32_t& num_solid_islands, const uint32_t read_length, uint32_t max_trimmed_bases, 
        uint32_t min_check_length, uint32_t& num_corrected_errors_local, const uint8_t* bit_vector, const uint32_t* hash_seed, init_args& arg, func_mem& funcm){
    //vector_candidate_path candidate_path_vector_tmp;
    //init_vector_candidate_path(candidate_path_vector_tmp);
    //path_query_text(kmer_initial, bit_vector, hash_seed, arg.kmer_length, arg.num_hash_func, arg.bit_vector_width, 1);

    //for each alternative neocletide 
    
    uint32_t index_last_mod = index_end - arg.kmer_length + 1;
    
  //printf("1 %d  %d\n", index_start, index_last_mod);
    bool run_exploration = true;
	/*
    char c_path_vec[(MAX_CANDIDATE_PATHS + 1) * READ_MAX_LEN];
    uint32_t pos_path_vec[(MAX_CANDIDATE_PATHS + 1)* READ_MAX_LEN];
	*/
	char *c_path_vec = funcm.cbuffer +  READ_MAX_LEN;
	uint32_t *pos_path_vec = funcm.ubuffer + 4 * READ_MAX_LEN;
    uint32_t c_path_size = 0;
    char *kmer_initial = funcm.cbuffer + (MAX_CANDIDATE_PATHS + 2) * READ_MAX_LEN;
    memcpy(kmer_initial, sequence_modified + index_start, arg.kmer_length);
    for(uint16_t it_alter = A; it_alter <= T; it_alter ++){
        kmer_initial[arg.kmer_length - 1] = NEOCLEOTIDE[it_alter];
        //new path 
        if(run_exploration == true) 
        if(query_text(kmer_initial, bit_vector, &hash_seed[0], arg.kmer_length, arg.num_hash_func, arg.bit_vector_width) == true){
			//printf("NE %c\n", NEOCLEOTIDE[it_alter]);
            if(index_start == index_last_mod){
				if(sequence_modified[index_start + arg.kmer_length - 1] != NEOCLEOTIDE[it_alter]){
                	pos_path_vec[c_path_size * READ_MAX_LEN + 0] = 1;
                	pos_path_vec[c_path_size * READ_MAX_LEN + 1] = index_start + arg.kmer_length - 1;
                	c_path_vec[c_path_size * READ_MAX_LEN + 1] = NEOCLEOTIDE[it_alter];
         		}else{	
                	pos_path_vec[c_path_size * READ_MAX_LEN + 0] = 0;
				}
                // this is a new path, should save
                c_path_size += 1;
            }else if((index_start < index_last_mod)){
                cpu_extend_a_kmer(
                        kmer_initial,
                        sequence_modified,
                        index_start,
                        index_last_mod,
                        pos_path_vec,
                        c_path_vec,
                        c_path_size,
                        org_boundary_left,
                        org_boundary_right,
                        quality_score,
                        bit_vector,
                        hash_seed,
                        run_exploration,
                        arg,
                        funcm
                        );
                
            }
        }
    }

  	/* 
    for(uint32_t it = 0; it < c_path_size; it ++){
        for(uint32_t ct = 0; ct < pos_path_vec[it * READ_MAX_LEN]; ct ++){
            printf("%d %c\n", pos_path_vec[it * READ_MAX_LEN + ct + 1], c_path_vec[it * READ_MAX_LEN + ct + 1]);
        }
        printf("****\n");
    }
   printf("********\n"); 
   */
    /*
    for(uint32_t it = 0; it < c_path_size; it ++){
        for(uint32_t ct = 0; ct < pos_path_vec[it * READ_MAX_LEN]; ct ++){
            printf("%d %c\n", pos_path_vec[it * READ_MAX_LEN + ct + 1], c_path_vec[it * READ_MAX_LEN + ct + 1]);
        }
    }
	*/
    if(run_exploration == false){
        c_path_size = 0;
    }
    
    /*
    for(uint32_t it = 0; it < c_path_size; it ++){
        for(uint32_t ct = 0; ct < pos_path_vec[it * READ_MAX_LEN]; ct ++){
            printf("%d %c\n", pos_path_vec[it * READ_MAX_LEN + ct + 1], c_path_vec[it * READ_MAX_LEN + ct + 1]);
        }
        printf("****\n");
    }
    printf("c_path_size %d %d \n", c_path_size, index_last_mod);    
    */
    bool all_solid_wo_modification = false;
    uint32_t c_path_num = 0;
    uint32_t flag = 0;
    for(uint32_t it_path = 0; it_path < c_path_size; it_path ++){
        if(pos_path_vec[it_path * READ_MAX_LEN] == 0){
            all_solid_wo_modification = true;
            break;
        }else{
            uint32_t index_last_modified_base = pos_path_vec[it_path * READ_MAX_LEN + pos_path_vec[it_path * READ_MAX_LEN] ];
            //printf(" ind %d %d\n", index_last_modified_base, index_last_mod);
            if(index_last_modified_base > index_last_mod){
                char *sequence_tmp = funcm.cbuffer + (MAX_CANDIDATE_PATHS + 2) * READ_MAX_LEN;
                memcpy(sequence_tmp, sequence_modified, read_length);

                for(uint32_t it_base = 0; it_base < pos_path_vec[it_path * READ_MAX_LEN]; it_base ++){
                    sequence_tmp[pos_path_vec[it_path * READ_MAX_LEN + it_base + 1]] = c_path_vec[it_path * READ_MAX_LEN + it_base + 1];
                }
                uint32_t num_success = 0;
                uint32_t num_fails = 0;
                for(uint32_t it_check = index_last_mod; it_check<= index_last_modified_base; it_check ++){
                    char kmer_current[READ_MAX_LEN];
                    memcpy(kmer_current, sequence_tmp + it_check, arg.kmer_length);
                    if(query_text(kmer_current, bit_vector, &hash_seed[0], arg.kmer_length, arg.num_hash_func, arg.bit_vector_width) == true){
                        num_success += 1;
                    }else{
                        num_fails += 1;
                        if(num_fails > NUM_ALLOWABLE_FAILS){
                            break;
                        }
                    }
                }
                //printf("success %d \n", num_success);
                if(num_success >= (index_last_modified_base - index_last_mod + 1 - NUM_ALLOWABLE_FAILS)){
                    //printf("%d %d\n", c_path_num, flag);
                    if(c_path_num != it_path){
                        memcpy(pos_path_vec + c_path_num * READ_MAX_LEN, pos_path_vec + it_path * READ_MAX_LEN, (pos_path_vec[it_path * READ_MAX_LEN] + 1) * sizeof(uint32_t));
                        memcpy(c_path_vec + c_path_num * READ_MAX_LEN, c_path_vec + it_path * READ_MAX_LEN, pos_path_vec[it_path * READ_MAX_LEN] + 1);
                    }
                    c_path_num += 1;
                }
            }else{
                if(c_path_num != it_path){
                    memcpy(pos_path_vec + c_path_num * READ_MAX_LEN, pos_path_vec + it_path * READ_MAX_LEN, (pos_path_vec[it_path * READ_MAX_LEN] + 1) * sizeof(uint32_t));
                    memcpy(c_path_vec + c_path_num * READ_MAX_LEN, c_path_vec + it_path * READ_MAX_LEN, pos_path_vec[it_path * READ_MAX_LEN] + 1);
                }
                c_path_num += 1;
            }
        }
    }
    c_path_size = c_path_num;
   /* 
    for(uint32_t it = 0; it < c_path_size; it ++){
        for(uint32_t ct = 0; ct < pos_path_vec[it * READ_MAX_LEN]; ct ++){
            printf("%d %c\n", pos_path_vec[it * READ_MAX_LEN + ct + 1], c_path_vec[it * READ_MAX_LEN + ct + 1]);
        }
        printf("****\n");
    }
    */
    c_path_num = 0; 
    for(uint32_t it_candidate = 0; it_candidate < c_path_size; it_candidate ++){
        bool real_modified = false;
        for(unsigned int it_mod_base = 0; it_mod_base < pos_path_vec[it_candidate * READ_MAX_LEN]; it_mod_base ++){
            if(org_sequence[pos_path_vec[it_candidate * READ_MAX_LEN + it_mod_base + 1]] != c_path_vec[it_candidate * READ_MAX_LEN + it_mod_base + 1]){
                real_modified = true;
            }
        }
        if(real_modified){
            if(c_path_num != it_candidate){
                memcpy(pos_path_vec + c_path_num * READ_MAX_LEN, pos_path_vec + it_candidate * READ_MAX_LEN, (pos_path_vec[it_candidate * READ_MAX_LEN] + 1) * sizeof(uint32_t));
                memcpy(c_path_vec + c_path_num * READ_MAX_LEN, c_path_vec + it_candidate * READ_MAX_LEN, pos_path_vec[it_candidate * READ_MAX_LEN] + 1);
            }
            c_path_num += 1;
        }
    }
    c_path_size = c_path_num;

    
   	/* 
    for(uint32_t it = 0; it < c_path_size; it ++){
        for(uint32_t ct = 0; ct < pos_path_vec[it * READ_MAX_LEN]; ct ++){
            printf("%d %c\n", pos_path_vec[it * READ_MAX_LEN + ct + 1], c_path_vec[it * READ_MAX_LEN + ct + 1]);
        }
        printf("!!!!!!\n");
    }
   printf("!!!!!!!!!!!!!!\n"); 
   */
    /*
    for(uint32_t ii = 0; ii < read_length; ii ++){
        printf("%c", sequence_modified[ii]);
    } 
    printf("\n");
    */
    if(all_solid_wo_modification){
        //do nothing 
    }else if (c_path_size > 1){
        //printf("*******\n");
        uint32_t it_path;
        uint32_t it_path_1st;
        uint32_t it_path_2nd;
        
        uint32_t qs_1st = INIT_MIN_QS;
        uint32_t qs_2nd = INIT_MIN_QS;
        
        for(it_path = 0; it_path < c_path_size; it_path ++){
            uint32_t sum_qs = 0;
            for(uint32_t it_mod = 0; it_mod < pos_path_vec[it_path * READ_MAX_LEN]; it_mod ++){
                sum_qs += quality_score[pos_path_vec[it_path * READ_MAX_LEN + it_mod + 1]] - arg.quality_score_offset;
            }
            if(sum_qs <= qs_1st){
                qs_2nd = qs_1st;
                qs_1st = sum_qs;
                it_path_2nd = it_path_1st;
                it_path_1st = it_path;
            }else if(sum_qs <= qs_2nd){
                qs_2nd = sum_qs;
                it_path_2nd = it_path;
            }
        }

        bool too_many_correction = false;
        
        if(num_solid_islands == 2){
            if(index_start - left_first > MAX_MODIFICATION && right_second - index_end + 2 <= MAX_MODIFICATION){
                if(pos_path_vec[it_path_1st * READ_MAX_LEN + pos_path_vec[it_path_1st * READ_MAX_LEN]] - arg.kmer_length - index_start + 2 >= min_check_length && pos_path_vec[it_path_1st * READ_MAX_LEN] > MAX_MODIFICATION){
                    for(uint32_t it_mod_base = pos_path_vec[it_path_1st * READ_MAX_LEN] - 1; it_mod_base >= MAX_MODIFICATION; it_mod_base --){
                        if(pos_path_vec[it_path_1st * READ_MAX_LEN + it_mod_base + 1] - pos_path_vec[it_path_1st * READ_MAX_LEN + it_mod_base + 1 - MAX_MODIFICATION] < min_check_length){
                            if(read_length - pos_path_vec[it_path_1st * READ_MAX_LEN + it_mod_base + 1]<= max_trimmed_bases){
                                trim_3_end = read_length - pos_path_vec[it_path_1st * READ_MAX_LEN + it_mod_base + 1];
                            }else if (pos_path_vec[it_path_1st * READ_MAX_LEN + it_mod_base + 1] + 1 <= max_trimmed_bases){
                                trim_5_end = pos_path_vec[it_path_1st * READ_MAX_LEN + it_mod_base + 1] + 1;
                            }

                            too_many_correction = true;
                            break;
                        }
                    }
                }
            }else if(index_start - left_first <= MAX_MODIFICATION && right_second - index_end + 2 > MAX_MODIFICATION){
                if(index_end - pos_path_vec[it_path_1st * READ_MAX_LEN + 0 + 1] + 1 >= min_check_length && pos_path_vec[it_path_1st * READ_MAX_LEN] > MAX_MODIFICATION){
                    for(uint32_t it_mod_base = 0; it_mod_base < pos_path_vec[it_path_1st * READ_MAX_LEN] - MAX_MODIFICATION; it_mod_base ++){
                        if(pos_path_vec[it_path_1st * READ_MAX_LEN + it_mod_base + 1 + MAX_MODIFICATION] - pos_path_vec[it_path_1st * READ_MAX_LEN + it_mod_base + 1] < min_check_length){ 
                            if(read_length - pos_path_vec[it_path_1st * READ_MAX_LEN + it_mod_base + 1]<= max_trimmed_bases){
                                trim_3_end = read_length - pos_path_vec[it_path_1st * READ_MAX_LEN + it_mod_base + 1];
                            }else if (pos_path_vec[it_path_1st * READ_MAX_LEN + it_mod_base + 1] + 1 <= max_trimmed_bases){
                                trim_5_end = pos_path_vec[it_path_1st * READ_MAX_LEN + it_mod_base + 1] + 1;
                            }
                            
                            too_many_correction = true;
                            break;
                        }
                    }
                }
            }
        }
    /*
    for(uint32_t ii = 0; ii < read_length; ii ++){
        printf("%c", sequence_modified[ii]);
    } 
    printf("\n");
    printf("%d %d %d %d %d %d %d\n", trim_5_end, trim_3_end, it_path_1st, qs_1st, qs_2nd, too_many_correction, max_trimmed_bases);
    */    
    cpu_modify_sequence(org_sequence, sequence_modified, quality_score, sequence_modification, trim_5_end, trim_3_end, 
         read_length, num_corrected_errors_local, pos_path_vec, c_path_vec, c_path_size, arg, funcm, it_path_1st, qs_1st,
         qs_2nd, too_many_correction, max_trimmed_bases, 0);

    } else if(c_path_size == 1){
        bool too_many_correction = false;
        if(num_solid_islands == 2){
            if(index_start - left_first > MAX_MODIFICATION && right_second - index_end + 2 <= MAX_MODIFICATION){
                if(pos_path_vec[pos_path_vec[0]] - arg.kmer_length - index_start + 2 >= min_check_length && pos_path_vec[0] > MAX_MODIFICATION){
                    for(uint32_t it_mod_base = pos_path_vec[0] - 1; it_mod_base >= MAX_MODIFICATION; it_mod_base --){
                        if(pos_path_vec[it_mod_base + 1] - pos_path_vec[it_mod_base + 1 - MAX_MODIFICATION] < min_check_length){
                            if(read_length - pos_path_vec[it_mod_base + 1] <= max_trimmed_bases){
                                trim_3_end = read_length - pos_path_vec[it_mod_base + 1];
                                if(it_mod_base > 0){
                                    for(uint32_t it_base = 0;it_base < it_mod_base - 1; it_base ++){
                                        sequence_modification[pos_path_vec[it_base + 1]] = c_path_vec[it_base + 1];
                                        sequence_modified[pos_path_vec[it_base +1]] = c_path_vec[it_base + 1];
                                        num_corrected_errors_local += 1;
                                    }
                                }
                            }
                            else if(pos_path_vec[it_mod_base + 1] + 1 <= max_trimmed_bases){
                                trim_5_end = pos_path_vec[it_mod_base + 1] + 1;
                                if(it_mod_base < pos_path_vec[0] - 1){
                                    for(uint32_t it_base = it_mod_base + 1; it_base < pos_path_vec[0]; it_base ++){
                                        sequence_modification[pos_path_vec[it_base + 1]] = c_path_vec[it_base + 1];
                                        sequence_modified[pos_path_vec[it_base + 1]] = c_path_vec[it_base + 1];
                                        num_corrected_errors_local += 1;
                                    }
                                }
                            }
                            too_many_correction = true;
                            break;
                        }
                    }
                }
            }else if(index_start - left_first <=  MAX_MODIFICATION && right_second - index_end + 2 >  MAX_MODIFICATION){
                if(index_end - pos_path_vec[0 + 1] + 1 >= min_check_length && pos_path_vec[0] > MAX_MODIFICATION){
                    for(uint32_t it_mod_base = 0; it_mod_base < pos_path_vec[0] - MAX_MODIFICATION; it_mod_base ++){
                        if(pos_path_vec[it_mod_base + 1 + MAX_MODIFICATION] - pos_path_vec[it_mod_base + 1] < min_check_length){
                            if(read_length - pos_path_vec[it_mod_base + 1] <= max_trimmed_bases){
                                trim_3_end = read_length - pos_path_vec[it_mod_base + 1];
                            }
                            else if(pos_path_vec[it_mod_base + 1] + 1 <= max_trimmed_bases){
                                trim_5_end = pos_path_vec[it_mod_base + 1] + 1;
                            }
                            too_many_correction = true;
                            break;
                        }
                    }
                }
            }
        }
        if(too_many_correction == false){
            for(uint32_t it_base = 0; it_base < pos_path_vec[0]; it_base ++){
                sequence_modification[pos_path_vec[it_base + 1]] = c_path_vec[it_base + 1];
                sequence_modified[pos_path_vec[it_base + 1]] = c_path_vec[it_base + 1];
                num_corrected_errors_local += 1;
            }
        }
    }
    else  if(c_path_size == 0){
       if(read_length - (index_start + arg.kmer_length - 1) <= max_trimmed_bases){
            trim_3_end = read_length - (index_start + arg.kmer_length - 1);
       } 
        else if(index_end + 1 <= max_trimmed_bases){
            trim_5_end = index_end + 1;
        }
        else if(read_length - (org_boundary_left + arg.kmer_length - 1) <= max_trimmed_bases){
            trim_3_end = read_length - (org_boundary_left + arg.kmer_length - 1);
        }else if(org_boundary_right + 1 <=  max_trimmed_bases){
           trim_5_end = org_boundary_right + 1; 
        }
    }
  


}   

 void cpu_extend_a_kmer_5_prime_end(const char* kmer, const char* sequence, uint32_t index_kmer, uint32_t* pos_path_vec,
        char* c_path_vec, uint32_t& c_path_size,  const uint32_t& org_boundary, const char* quality_score, const uint8_t* bit_vector, 
        const uint32_t* hash_seed, bool& run_exploration, init_args& arg, func_mem& funcm){
    
    //char path_content[4 * READ_MAX_LEN];
    char *path_content = funcm.cbuffer + (MAX_CANDIDATE_PATHS + 1 + 2) * READ_MAX_LEN;
    int index_pos = 0;
    //char* kmer_ext = kmer + 1;
    //char  kmer_new[36];
    char *kmer_new = funcm.cbuffer + (MAX_CANDIDATE_PATHS + 1 + 6) * READ_MAX_LEN;
    //uint32_t rec_index[READ_MAX_LEN];
    uint32_t *rec_index = funcm.ubuffer + (MAX_CANDIDATE_PATHS + 1 + 4) * READ_MAX_LEN;
    //char kmer_ext[READ_MAX_LEN];
    char *kmer_ext = funcm.cbuffer + (MAX_CANDIDATE_PATHS + 1 + 7) * READ_MAX_LEN;
    uint32_t index_start = index_kmer;
    //printf("index start %d \n", index_start);
    memcpy(kmer_ext + index_kmer, kmer, arg.kmer_length - 1);
    /*
    for(uint32_t i = 0; i < 101; i ++){
        printf("%c", sequence[i]);
    }
    printf("\n");
    */
    uint32_t flag = 1;
    while(index_pos > 0 || flag == 1){
        flag = 0;
        memcpy(kmer_new + 1, kmer_ext + index_kmer, arg.kmer_length - 1);
        //kmer_new[arg.kmer_length] = sequence[index_kmer - 1];
        //path_query_text(kmer_new, bit_vector, hash_seed, arg.kmer_length, arg.num_hash_func, arg.bit_vector_width, 1);
        /* 
        for(int ii = 0; ii < index_pos; ii ++){
            for(uint32_t jj = 0; jj < 4; jj ++){
                if(path_content[(ii << 2) + jj] > 0){
                    printf("%c", path_content[(ii << 2) + jj]);
                }else{
                    break;
                }
            }
            printf(" ");
        }
        printf("\n");
        for(uint32_t ii = index_kmer; ii < index_start + arg.kmer_length - 1; ii ++){
            printf("%c", kmer_ext[ii]);
        }
        printf("\n");
        for(int ii = 0; ii < index_pos; ii ++){
            printf("%d ", rec_index[ii]);
        }
        printf("\n");
        for(uint32_t ii = 1; ii < arg.kmer_length; ii ++){
            printf("%c", kmer_new[ii]);
        }
        printf("\n"); 
        printf("%d\n", index_kmer);
        printf("*****\n");
        char cc, ccc;
        scanf("%c%c", &cc, &ccc);
        printf("cc %c \n", cc); 
        */
        uint32_t k_ind = 0;

        if ((index_kmer == (org_boundary + 1)) || ((quality_score[index_kmer - 1] - arg.quality_score_offset) <= arg.extremely_low_quality_score)){
            for(unsigned short int it_alter = A; it_alter <= T; it_alter ++){ 
                kmer_new[0] = NEOCLEOTIDE[it_alter];
                //if(kmer_new[funcm.kmer_length + it_alter]){
                if(query_text(kmer_new, bit_vector, &hash_seed[0], arg.kmer_length, arg.num_hash_func, arg.bit_vector_width) == true){
                    path_content[(index_pos << 2) + k_ind ++] = NEOCLEOTIDE[it_alter];      
                }
            }                
        }else{
            kmer_new[0] = sequence[index_kmer - 1];
            if(query_text(kmer_new, bit_vector, &hash_seed[0], arg.kmer_length, arg.num_hash_func, arg.bit_vector_width) == true){ 
                path_content[(index_pos << 2) + k_ind ++] = kmer_new[0];      
            }else{
                 for(unsigned short int it_alter = A; it_alter <= T; it_alter ++){ 
                    kmer_new[0] = NEOCLEOTIDE[it_alter];
                    //if(kmer_new[funcm.kmer_length + it_alter]){
                    if(sequence[index_kmer - 1] != NEOCLEOTIDE[it_alter])
                    if(query_text(kmer_new, bit_vector, &hash_seed[0], arg.kmer_length, arg.num_hash_func, arg.bit_vector_width) == true){
                        path_content[(index_pos << 2) + k_ind ++] = NEOCLEOTIDE[it_alter];      
                    }
                }      
            }
        
        }
        
        if(k_ind < 4)
            path_content[(index_pos << 2) + k_ind] = 0;
        if(index_kmer  == 1 && k_ind > 0 ){
            for(uint32_t it = 0; it < k_ind; it ++){
                //save path 
                pos_path_vec[c_path_size * READ_MAX_LEN] = 0; 
                kmer_ext[index_kmer - 1] = path_content[(index_pos << 2) + it];          
               	/* 
                for(uint32_t ii = 0; ii <= index_start; ii ++){
                    printf("%x ", kmer_ext[ii]);
                }   
                printf("\n");
                */
                for(int ct = index_start; ct >= 0; ct --){
                    if(sequence[ct] != kmer_ext[ct]){
                        pos_path_vec[c_path_size * READ_MAX_LEN + pos_path_vec[c_path_size * READ_MAX_LEN] + 1] = ct;
                        c_path_vec[c_path_size * READ_MAX_LEN + pos_path_vec[c_path_size * READ_MAX_LEN] + 1] = kmer_ext[ct];
                        pos_path_vec[c_path_size * READ_MAX_LEN] += 1;
                    }
                }
                c_path_size += 1;        

               	if(c_path_size > MAX_CANDIDATE_PATHS){
                    //退出寻找路径
                   	run_exploration = false;
                    index_pos = -1;     
                    break;  
                 } 
                 
            }  
            while(index_pos > 0){
                if(rec_index[index_pos - 1] < 4 && path_content[((index_pos - 1) << 2 ) + rec_index[index_pos - 1]] > 0){
                    kmer_ext[index_kmer] = path_content[((index_pos - 1) << 2) + rec_index[index_pos - 1]];
                    rec_index[index_pos - 1] += 1;
                    break;
                }else{
                    index_kmer += 1;
                    index_pos -=1;
                }
            }

        }
        else if (k_ind){
            index_kmer -= 1;
            rec_index[index_pos] = 0;
            kmer_ext[index_kmer] = path_content[(index_pos << 2) + rec_index[index_pos]];
            rec_index[index_pos] += 1;
            index_pos += 1;              
        }else{
            while(index_pos > 0){
                if(rec_index[index_pos - 1] < 4 && path_content[((index_pos - 1) << 2 ) + rec_index[index_pos - 1]]){
                    kmer_ext[index_kmer] = path_content[((index_pos - 1) << 2) + rec_index[index_pos - 1]];
                    rec_index[index_pos - 1] += 1;
                    break;
                }else{
                    index_kmer += 1;
                    index_pos -=1;
                }
            }
            //printf("index_pos %d \n", index_pos);
        }
               

    }

}

 void cpu_correct_errors_5_prime_end(const char* org_sequence, char* sequence, const char* quality_score, const uint32_t& index_start, char* sequence_modification, uint32_t& trim_5_end,
       uint32_t& trim_3_end, const uint32_t& org_boundary, const uint32_t& read_length, const uint32_t& max_trimmed_bases, const uint32_t& min_check_length, uint32_t& num_corrected_errors_local,
    const uint8_t* bit_vector, const uint32_t* hash_seed, init_args& arg, func_mem& funcm){
    //vector_candidate_path candidate_path_vector_tmp;
    //init_vector_candidate_path(candidate_path_vector_tmp);
    
    bool run_exploration = true;
    /*
    char c_path_vec[(MAX_CANDIDATE_PATHS  + 1) * READ_MAX_LEN];
    uint32_t pos_path_vec[(MAX_CANDIDATE_PATHS + 1) * READ_MAX_LEN];
    */
    char *c_path_vec = funcm.cbuffer + 1 * READ_MAX_LEN;
    uint32_t *pos_path_vec = funcm.ubuffer + 4 * READ_MAX_LEN; 
    uint32_t c_path_size = 0;
    char* kmer_initial = funcm.cbuffer + (MAX_CANDIDATE_PATHS + 1 + 1) * READ_MAX_LEN;
    memcpy(kmer_initial, sequence + index_start, arg.kmer_length);
    //for(uint32_t ii = 0; ii < arg.kmer_length; ii ++){
    //    printf("%c", kmer_initial[ii]);
    //}
    //printf("******\n");
    //path_query_text(kmer_initial, bit_vector, hash_seed, arg.kmer_length, arg.num_hash_func, arg.bit_vector_width, 0);
	//printf("%d \n", index_start);
    for(unsigned short int it_alter = A; it_alter <= T; it_alter ++){
        kmer_initial[0] = NEOCLEOTIDE[it_alter];
        //if(kmer_initial[funcm.kmer_length + it_alter]){
        if(run_exploration == true)
        if(query_text(kmer_initial, bit_vector, &hash_seed[0], arg.kmer_length, arg.num_hash_func, arg.bit_vector_width) == true){
              //generate new path
            
            if(index_start == 0){
                //new path  
				if(sequence[index_start] != NEOCLEOTIDE[it_alter]){
                	pos_path_vec[c_path_size * READ_MAX_LEN] = 1;
                	pos_path_vec[c_path_size * READ_MAX_LEN + 1] = index_start;
                	c_path_vec[c_path_size * READ_MAX_LEN + 1] = NEOCLEOTIDE[it_alter];              
            	}else{
                	pos_path_vec[c_path_size * READ_MAX_LEN] = 0;
				}
                c_path_size += 1;              
           } else if(index_start > 0){
                cpu_extend_a_kmer_5_prime_end(
                    kmer_initial,
                    sequence,
                    index_start,
                    pos_path_vec,
                    c_path_vec,
                    c_path_size,
                    org_boundary,
                    quality_score,
                    bit_vector,
                    hash_seed,
                    run_exploration,
                    arg,
                    funcm
                    );
           }          
        }
    }    
    
    //printf("hhh %d \n", c_path_size);
   	/* 
    for(uint32_t ii = 0; ii < c_path_size; ii ++){
        for(uint32_t jj = 0; jj < pos_path_vec[ii * READ_MAX_LEN]; jj ++){
            printf("%d %c\n", pos_path_vec[ii * READ_MAX_LEN + jj + 1], c_path_vec[ii * READ_MAX_LEN + jj + 1]);
        }
        printf("****\n");
    }
    
  	for(uint32_t ii = 0; ii < 101; ii ++){
		printf("%c", sequence[ii]);
	}	
	printf("\n");
   */
    if(run_exploration == false){
        c_path_size = 0;
    }
	    
    uint32_t c_path_num = 0;
    bool all_solid_wo_modification = false;
    for(uint32_t it_path = 0; it_path < c_path_size; it_path ++){
        if(pos_path_vec[it_path * READ_MAX_LEN] == 0){
            all_solid_wo_modification = true;
            break;
        }
        else{
            char *sequence_tmp = funcm.cbuffer + (MAX_CANDIDATE_PATHS + 1 + 1) * READ_MAX_LEN;
            memcpy(sequence_tmp, sequence, read_length);
            uint32_t index_smallest_modified = pos_path_vec[it_path * READ_MAX_LEN + pos_path_vec[it_path * READ_MAX_LEN]];
            uint32_t extend_amount;
            if(index_smallest_modified >= arg.kmer_length - 1){
                if(c_path_num != it_path){
                    memcpy(pos_path_vec + c_path_num * READ_MAX_LEN, pos_path_vec + it_path * READ_MAX_LEN, (pos_path_vec[it_path * READ_MAX_LEN] + 1) * sizeof(uint32_t));
                    memcpy(c_path_vec + c_path_num * READ_MAX_LEN, c_path_vec + it_path * READ_MAX_LEN, pos_path_vec[it_path * READ_MAX_LEN] + 1);
                }
                c_path_num += 1;
            }
            else{
                for(uint32_t it_base = 0; it_base < pos_path_vec[it_path * READ_MAX_LEN]; it_base ++){
                     sequence_tmp[pos_path_vec[it_path * READ_MAX_LEN + it_base + 1]] = c_path_vec[it_path * READ_MAX_LEN + it_base + 1];
                }
                if(index_smallest_modified >= arg.kmer_length - arg.max_extension -1){
                    extend_amount = arg.kmer_length - index_smallest_modified - 1;
                }
                else{
                    extend_amount = arg.max_extension;
                }
                bool extension = false;
               char *kmer_initial = funcm.cbuffer + (MAX_CANDIDATE_PATHS + 1 + 2) * READ_MAX_LEN;
                char *kmer_ext = funcm.cbuffer + (MAX_CANDIDATE_PATHS + 1 + 3) * READ_MAX_LEN;
                memcpy(kmer_ext + extend_amount, sequence_tmp, arg.kmer_length - 1);
                uint32_t *rec_index = funcm.ubuffer + (MAX_CANDIDATE_PATHS + 1 + 4) * READ_MAX_LEN;
                char *path_content = funcm.cbuffer + (MAX_CANDIDATE_PATHS  +  1 + 4) * READ_MAX_LEN;
                int index_pos = 0;
                uint32_t num_extend = 0;
                uint32_t flag = 1;
                //printf("%d \n", extend_amount);
                while(index_pos > 0 || flag == 1){
                    flag = 0;
                    memcpy(kmer_initial + 1, kmer_ext + extend_amount - num_extend, arg.kmer_length - 1);
                   //path_query_text(kmer_initial, bit_vector, hash_seed, arg.kmer_length, arg.num_hash_func, arg.bit_vector_width, 1);
                    uint32_t k_ind = 0;
                    for(uint32_t it_alter = A; it_alter <= T; it_alter ++){
                        //if(kmer_initial[funcm.kmer_length + it_alter])
                        kmer_initial[0] = NEOCLEOTIDE[it_alter];
                        if(query_text(kmer_initial, bit_vector, hash_seed, arg.kmer_length, arg.num_hash_func, arg.bit_vector_width) == true)
                            path_content[(index_pos << 2) + k_ind ++] = NEOCLEOTIDE[it_alter];
                    }
                    if(k_ind < 4)
                        path_content[(index_pos << 2) + k_ind] = 0;
                    if(num_extend + 1 == extend_amount && k_ind > 0){
                        extension = true;
                        break;
                    }else if(k_ind > 0){
                        num_extend += 1;
                        rec_index[index_pos] = 0;
                        kmer_ext[extend_amount - index_pos - 1] = path_content[index_pos << 2];
                        rec_index[index_pos] = 1;
                        index_pos += 1;
                    } else{
                        while(index_pos > 0){
                            if(rec_index[index_pos - 1] < 4  && path_content[((index_pos - 1 ) << 2) + rec_index[index_pos - 1]]){
                                kmer_ext[extend_amount - index_pos] = path_content[((index_pos - 1) << 2) + rec_index[index_pos - 1]];
                                rec_index[index_pos - 1] += 1;
                                break;
                            }else{
                                num_extend -= 1;
                                index_pos -=1;
                            }
                        }
                    }
                }
                if(extension == true){
                    if(c_path_num != it_path){
                        memcpy(pos_path_vec + c_path_num * READ_MAX_LEN, pos_path_vec + it_path * READ_MAX_LEN, (pos_path_vec[it_path * READ_MAX_LEN] + 1) * sizeof(uint32_t));
                        memcpy(c_path_vec + c_path_num * READ_MAX_LEN, c_path_vec + it_path * READ_MAX_LEN, pos_path_vec[it_path * READ_MAX_LEN] + 1);
                    }
                    c_path_num += 1;
                }
            }
        }
    }
     
    c_path_size = c_path_num;
    c_path_num = 0;
    for(uint32_t it_candidate = 0; it_candidate < c_path_size; it_candidate ++){
        bool real_modified = false;
        for(uint32_t it_mod_base = 0; it_mod_base < pos_path_vec[it_candidate * READ_MAX_LEN]; it_mod_base ++){
            if(org_sequence[pos_path_vec[it_candidate * READ_MAX_LEN + it_mod_base + 1]] != c_path_vec[it_candidate * READ_MAX_LEN + it_mod_base + 1]){
                real_modified = true;
            }
        }
        if(real_modified){
            if(it_candidate != c_path_num){
                memcpy(pos_path_vec + c_path_num * READ_MAX_LEN, pos_path_vec + it_candidate* READ_MAX_LEN, (pos_path_vec[it_candidate * READ_MAX_LEN] + 1) * sizeof(uint32_t));
                memcpy(c_path_vec + c_path_num * READ_MAX_LEN, c_path_vec + it_candidate * READ_MAX_LEN, pos_path_vec[it_candidate * READ_MAX_LEN] + 1);
            }
            c_path_num += 1; 
        }
    }

	//printf("cpath size %d cpath %d\n", c_path_size, c_path_num);
    c_path_size = c_path_num;
   
    

    if(all_solid_wo_modification == true){
    
    } else if(c_path_size > 1){
        uint32_t it_path;
        uint32_t it_path_1st;
        uint32_t it_path_2nd;
        
        uint32_t qs_1st = INIT_MIN_QS;
        uint32_t qs_2nd = INIT_MIN_QS;
        
        for(it_path = 0; it_path < c_path_size; it_path ++){
            uint32_t sum_qs = 0;
            for(uint32_t it_mod = 0; it_mod < pos_path_vec[it_path * READ_MAX_LEN]; it_mod ++){
                if(sequence[pos_path_vec[it_path * READ_MAX_LEN + it_mod + 1]] != c_path_vec[it_path * READ_MAX_LEN + it_mod + 1]){ 
                    sum_qs += (uint32_t)quality_score[pos_path_vec[it_path * READ_MAX_LEN + it_mod + 1]] - arg.quality_score_offset;
                }
            }
            if(sum_qs <= qs_1st){
                qs_2nd = qs_1st;
                qs_1st = sum_qs;
                it_path_2nd = it_path_1st;
                it_path_1st = it_path;
            }else if(sum_qs <= qs_2nd){
                qs_2nd = sum_qs;
                it_path_2nd = it_path;
            }
        }
        //printf("it_path_1st %d \n", it_path_1st);
        bool too_many_correction = false;
       

        if(pos_path_vec[pos_path_vec[it_path_1st * READ_MAX_LEN]] + 1 >= min_check_length && pos_path_vec[it_path_1st * READ_MAX_LEN] > MAX_MODIFICATION){
            for(uint32_t it_mod_base = pos_path_vec[it_path_1st * READ_MAX_LEN] - 1; it_mod_base >= MAX_MODIFICATION; it_mod_base --){
                if(pos_path_vec[it_path_1st * READ_MAX_LEN + it_mod_base + 1] - pos_path_vec[it_path_1st * READ_MAX_LEN + it_mod_base + 1 - MAX_MODIFICATION] < min_check_length){
                    if(read_length - pos_path_vec[it_path_1st * READ_MAX_LEN + it_mod_base + 1] <= max_trimmed_bases){
                        trim_3_end = read_length - pos_path_vec[it_path_1st * READ_MAX_LEN + it_mod_base + 1];
                    }else if(pos_path_vec[it_path_1st * READ_MAX_LEN + it_mod_base + 1] + 1 <= max_trimmed_bases){
                        trim_5_end = pos_path_vec[it_path_1st * READ_MAX_LEN + it_mod_base + 1];
                    }
                    too_many_correction = true;
                    break;
                }
            }
        }
        
    
        cpu_modify_sequence(org_sequence, sequence, quality_score, sequence_modification, trim_5_end, trim_3_end, 
         read_length, num_corrected_errors_local, pos_path_vec, c_path_vec, c_path_size, arg, funcm, it_path_1st, 
         qs_1st, qs_2nd, too_many_correction, max_trimmed_bases, 0); 
     } else if(c_path_size == 1){
        bool too_many_correction = false;
        if(pos_path_vec[0 + 1] + 1 >= min_check_length && pos_path_vec[0] > MAX_MODIFICATION){
           for(uint32_t it_mod_base = pos_path_vec[0] - 1; it_mod_base >= MAX_MODIFICATION; it_mod_base --){
                if(pos_path_vec[it_mod_base + 1] - pos_path_vec[it_mod_base + 1 - MAX_MODIFICATION] < min_check_length){
                    if(read_length - pos_path_vec[it_mod_base + 1] <= max_trimmed_bases){
                        trim_3_end = read_length - pos_path_vec[it_mod_base + 1];
                        
                        if(it_mod_base > 0){
                            for(uint32_t it_base = 0; it_base < (it_mod_base - 1); it_base ++){
                                sequence_modification[pos_path_vec[it_base + 1]] = c_path_vec[it_base + 1];
                                sequence[pos_path_vec[it_base + 1]] = c_path_vec[it_base + 1];
                                num_corrected_errors_local ++;
                            }
                        }
                    }else if(pos_path_vec[it_mod_base + 1]  + 1 <= max_trimmed_bases){
                        trim_5_end = pos_path_vec[it_mod_base + 1] + 1;
                        if(it_mod_base < pos_path_vec[0] -1 ){
                            for(uint32_t it_base = it_mod_base + 1; it_base < pos_path_vec[0]; it_base ++){
                                sequence_modification[pos_path_vec[it_base + 1]] = c_path_vec[it_base + 1];
                                sequence[pos_path_vec[it_base + 1]] = c_path_vec[it_base + 1];
                                num_corrected_errors_local ++;
                            }
                        }
                    }
                    too_many_correction = true;
                    break;
                }
           }
        }
        if(too_many_correction == false){
            for(uint32_t it_base = 0; it_base < pos_path_vec[0]; it_base ++){
                sequence_modification[pos_path_vec[it_base + 1]] = c_path_vec[it_base + 1];
                sequence[pos_path_vec[it_base + 1]] = c_path_vec[it_base + 1];
                num_corrected_errors_local ++;
            }
        } 
    }
    
}
 inline void cpu_extend_a_kmer_3_prime_end(const char* kmer, const char* sequence, uint32_t index_kmer, uint32_t* pos_path_vec, 
        char* c_path_vec, uint32_t& c_path_size, const uint32_t& org_boundary, const char* quality_score, const uint8_t* bit_vector, 
        uint32_t* hash_seed, bool& run_exploration, uint32_t read_length, init_args& arg, func_mem& funcm){

    //char path_content[4 * READ_MAX_LEN];
    int index_pos = 0;
    //char* kmer_ext = kmer + 1;
    //char  kmer_new[36];
    //uint32_t rec_index[READ_MAX_LEN];
    //char kmer_ext[READ_MAX_LEN];

	char *path_content = funcm.cbuffer + (MAX_CANDIDATE_PATHS + 1 + 2) * READ_MAX_LEN;
	char *kmer_ext = funcm.cbuffer + (MAX_CANDIDATE_PATHS + 1 + 6) * READ_MAX_LEN;
	uint32_t *rec_index = funcm.ubuffer + (MAX_CANDIDATE_PATHS + 1 + 4) * READ_MAX_LEN;
    memcpy(kmer_ext + index_pos, kmer + 1, arg.kmer_length - 1);
    char *kmer_new = funcm.cbuffer + (MAX_CANDIDATE_PATHS + 1 + 7) * READ_MAX_LEN;
    uint32_t flag = 1;
    uint32_t index_start = index_kmer;
	/*
    for(uint32_t ii = 0; ii < arg.kmer_length + index_pos - 1; ii ++){
         printf("%c", kmer_ext[ii]);
     }
     printf("\n");
	 */
    //printf("hhhhhh%d\n", index_kmer);
    while(index_pos > 0 || flag == 1){
        flag = 0;
        memcpy(kmer_new, kmer_ext + index_pos, arg.kmer_length - 1);
       /* 
        for(uint32_t ii = 0; ii < arg.kmer_length + index_pos - 1; ii ++){
            printf("%c", kmer_ext[ii]);
        }
        printf("\n");
        
        
        for(uint32_t ii = 0; ii < arg.kmer_length - 1; ii ++){
            printf("%c", kmer_new[ii]);
        }
        printf("\n");
        */
        /*
        kmer_new[arg.kmer_length] = sequence[index_kmer + arg.kmer_length];
        path_query_text(kmer_new, bit_vector, hash_seed, arg.kmer_length, arg.num_hash_func, arg.bit_vector_width, 1);
        uint32_t k_ind = 0;
        for(unsigned short int it_alter = A; it_alter <= T; it_alter ++){ 
            if(kmer_new[funcm.kmer_length + it_alter]){
                path_content[(index_pos << 2) + k_ind ++] = NEOCLEOTIDE[it_alter];
            }
                       
        }
        */
        /*
        for(int ii = 0; ii < index_pos; ii ++){
            for(uint32_t jj = 0; jj < 4; jj ++){
                if(path_content[(ii << 2) + jj] > 0){
                    printf("%c", path_content[(ii << 2) + jj]);
                }else{
                    break;
                }
            }
            printf(" ");
        }
        printf("\n");
        for(uint32_t ii = 0; ii < index_pos + arg.kmer_length - 1; ii ++){
            printf("%c", kmer_ext[ii]);
        }
        printf("\n");
        for(int ii = 0; ii < index_pos; ii ++){
            printf("%d ", rec_index[ii]);
        }
        printf("\n");
        for(uint32_t ii = 1; ii < arg.kmer_length; ii ++){
            printf("%c", kmer_new[ii]);
        }
        printf("\n"); 
        printf("%d\n", index_kmer);
        printf("*****\n");
        char cc, ccc;
        scanf("%c%c", &cc, &ccc);
        printf("cc %c \n", cc); 
        */
        uint32_t k_ind = 0;
        if ((index_kmer == (org_boundary - 1)) || (((uint32_t)quality_score[index_kmer + arg.kmer_length] - arg.quality_score_offset) <= arg.extremely_low_quality_score)){
            for(unsigned short int it_alter = A; it_alter <= T; it_alter ++){ 
                kmer_new[arg.kmer_length - 1] = NEOCLEOTIDE[it_alter];
                //if(kmer_new[funcm.kmer_length + it_alter]){
                if(query_text(kmer_new, bit_vector, &hash_seed[0], arg.kmer_length, arg.num_hash_func, arg.bit_vector_width) == true){
                    path_content[(index_pos << 2) + k_ind ++] = NEOCLEOTIDE[it_alter];      
                }
            }                
        }else{
            kmer_new[arg.kmer_length - 1] = sequence[index_kmer + arg.kmer_length];
            if(query_text(kmer_new, bit_vector, &hash_seed[0], arg.kmer_length, arg.num_hash_func, arg.bit_vector_width) == true){ 
                path_content[(index_pos << 2) + k_ind ++] = kmer_new[arg.kmer_length - 1];      
            }else{
                 for(unsigned short int it_alter = A; it_alter <= T; it_alter ++){ 
                    kmer_new[arg.kmer_length - 1] = NEOCLEOTIDE[it_alter];
                    //if(kmer_new[funcm.kmer_length + it_alter]){
                    if(sequence[index_kmer + arg.kmer_length] != NEOCLEOTIDE[it_alter])
                    if(query_text(kmer_new, bit_vector, &hash_seed[0], arg.kmer_length, arg.num_hash_func, arg.bit_vector_width) == true){
                        path_content[(index_pos << 2) + k_ind ++] = NEOCLEOTIDE[it_alter];      
                    }
                }      
            } 
        }

        if(k_ind < 4)
            path_content[(index_pos << 2) + k_ind] = 0;
        if(index_kmer + 1 == read_length - arg.kmer_length && k_ind > 0){
            for(uint32_t it = 0; it < k_ind; it ++){
                //save path
                kmer_ext[arg.kmer_length - 1 + index_pos] = path_content[(index_pos << 2) + it];
                pos_path_vec[c_path_size * READ_MAX_LEN] = 0;
                /*
                for(uint32_t ii = 0; ii <= index_pos + 1; ii ++){
                    printf("%c", kmer_ext[arg.kmer_length -2 + ii]);
                }
                printf("\n");
                for(uint32_t ii = 0; ii <= index_pos + 1; ii ++){
                    printf("%c", sequence[arg.kmer_length - 1 + index_start + ii]);
                }
                printf("\n");
                //printf("*********\n");
                */
                for(uint32_t ct = 0; ct <= index_pos + 1; ct ++){
                    if(kmer_ext[arg.kmer_length - 2 + ct] != sequence[arg.kmer_length - 1 + index_start + ct]){
                        pos_path_vec[c_path_size * READ_MAX_LEN + pos_path_vec[c_path_size * READ_MAX_LEN] + 1] = index_start + ct + arg.kmer_length - 1 ;
                        c_path_vec[c_path_size * READ_MAX_LEN + pos_path_vec[c_path_size * READ_MAX_LEN] + 1] = kmer_ext[arg.kmer_length - 2 + ct];
                        pos_path_vec[c_path_size * READ_MAX_LEN] += 1;
                    }
                }

   /*             for(uint32_t ii = 0; ii < pos_path_vec[c_path_size * READ_MAX_LEN]; ii ++){
                    printf("%c", c_path_vec[c_path_size * READ_MAX_LEN + ii  + 1]);
                }
                printf("\n******\n");
    */            c_path_size += 1;        
                if(c_path_size > MAX_CANDIDATE_PATHS){
                    //退出寻找路径
                    index_pos = -1;       
                    run_exploration = false;
                    break;
                } 
            }  
            
            while(index_pos > 0){
                if(rec_index[index_pos - 1] < 4 && path_content[((index_pos - 1) << 2 ) + rec_index[index_pos - 1]] > 0){
                    kmer_ext[arg.kmer_length - 1 + index_pos - 1] = path_content[((index_pos - 1) << 2) + rec_index[index_pos - 1]];
                    rec_index[index_pos - 1] += 1;
                    break;
                }else{
                    index_kmer -= 1;
                    index_pos -=1;
                }
            }

        }
        else if (k_ind){
            index_kmer += 1;
            rec_index[index_pos] = 0;
            kmer_ext[arg.kmer_length - 1 + index_pos] = path_content[(index_pos << 2) + rec_index[index_pos]];
            rec_index[index_pos] += 1;
            index_pos += 1;              
        }else{
            while(index_pos > 0){
                if(rec_index[index_pos - 1] < 4 && path_content[((index_pos - 1) << 2) + rec_index[index_pos - 1]] > 0){
                    kmer_ext[arg.kmer_length - 1 + index_pos - 1] = path_content[((index_pos - 1) << 2) + rec_index[index_pos - 1]];
                    rec_index[index_pos - 1] += 1;
                    break;
                }else{
                    index_kmer -= 1;
                    index_pos -=1;
                }
            }
        }
    }
   
}

 void cpu_correct_errors_3_prime_end(char* org_sequence, char* sequence, char* quality_score, const uint32_t index_start, char* sequence_modification, uint32_t& trim_5_end, uint32_t& trim_3_end, const uint32_t org_boundary, const uint32_t read_length, const uint32_t max_trimmed_bases, const uint32_t min_check_length, uint32_t& num_corrected_errors_local, uint8_t* bit_vector,
        uint32_t* hash_seed, init_args& arg, func_mem& funcm){
	//printf("correct 3 prime end\n");    
   	/* 
    char c_path_vec[(MAX_CANDIDATE_PATHS + 1) * READ_MAX_LEN];
    uint32_t pos_path_vec[(MAX_CANDIDATE_PATHS + 1)* READ_MAX_LEN];\
	*/
	char *c_path_vec = funcm.cbuffer + READ_MAX_LEN;
	uint32_t *pos_path_vec = funcm.ubuffer + 4 * READ_MAX_LEN;
    uint32_t c_path_size = 0;
    //char kmer_initial[36];
	char *kmer_initial = funcm.cbuffer + (MAX_CANDIDATE_PATHS + 1 + 1) * READ_MAX_LEN;
    memcpy(kmer_initial, sequence + index_start, arg.kmer_length); 
    /*
    for(uint32_t ii = 0; ii < read_length; ii ++){
        printf("%c", sequence[ii]);
    }
    printf("\n");
    */
    bool run_exploration = true;
    //path_query_text(kmer_initial, bit_vector, hash_seed, arg.kmer_length, arg.num_hash_func, arg.bit_vector_width, 0);
    for(uint32_t it_alter = A; it_alter <= T; it_alter ++){
        kmer_initial[arg.kmer_length - 1] = NEOCLEOTIDE[it_alter];
        /* 
                    for(uint32_t ii = 0; ii < arg.kmer_length; ii ++){
                        printf("%c", kmer_initial[ii]);
                    }
                    printf("\n");
        */
        //if(kmer_initial[funcm.kmer_length + it_alter]){
        if(run_exploration == true)
        if(query_text(kmer_initial, bit_vector, hash_seed, arg.kmer_length, arg.num_hash_func, arg.bit_vector_width) == true){
            //printf("hhhhhh\n"); 
			//printf("mmmmm %c\n", NEOCLEOTIDE[it_alter]);
            if(index_start == read_length - arg.kmer_length){
                //generate new path
				if(sequence[index_start + arg.kmer_length - 1] != NEOCLEOTIDE[it_alter]){
                    pos_path_vec[c_path_size * READ_MAX_LEN] = 1;
                    pos_path_vec[c_path_size * READ_MAX_LEN + 1] = index_start + arg.kmer_length - 1;
                    c_path_vec[c_path_size * READ_MAX_LEN + 1 ] = NEOCLEOTIDE[it_alter];
				}else{	
                    pos_path_vec[c_path_size * READ_MAX_LEN] = 0;
				}
                //save new path        
                c_path_size += 1;
            }else if(index_start < read_length - arg.kmer_length){
                cpu_extend_a_kmer_3_prime_end(
                        kmer_initial,
                        sequence,
                        index_start,
                        pos_path_vec,
                        c_path_vec,
                        c_path_size,
                        org_boundary,
                        quality_score,
                        bit_vector,
                        hash_seed,
                        run_exploration,
                        read_length,
                        arg,
                        funcm
                        );
            }
        }
    } 
  /* 
   for(uint32_t ii = 0; ii < c_path_size; ii ++){
        for(uint32_t jj = 0; jj < pos_path_vec[ii * READ_MAX_LEN]; jj ++){
            printf("%d %c\n", pos_path_vec[ii * READ_MAX_LEN + jj + 1], c_path_vec[ii * READ_MAX_LEN + jj + 1]);
        }
        printf("****\n");
   }
	*/
    //printf("******\n");     
    
    if(run_exploration == false){
        c_path_size = 0;
    } 
    uint32_t c_path_num = 0;
    bool all_solid_wo_modification = false;

    for(uint32_t it_path = 0; it_path < c_path_size; it_path ++){
        if(pos_path_vec[it_path * READ_MAX_LEN] == 0){
            all_solid_wo_modification = true;
            break;
        }
        else{
            char *sequence_tmp = funcm.cbuffer + (MAX_CANDIDATE_PATHS + 2) * READ_MAX_LEN;
            memcpy(sequence_tmp, sequence, read_length);
            uint32_t index_largest_modified = pos_path_vec[it_path * READ_MAX_LEN + pos_path_vec[it_path * READ_MAX_LEN]];
            uint32_t extend_amount;
            if(index_largest_modified <= read_length - arg.kmer_length){
                if(c_path_num != it_path){
                    memcpy(pos_path_vec + c_path_num * READ_MAX_LEN, pos_path_vec + it_path * READ_MAX_LEN, (pos_path_vec[it_path * READ_MAX_LEN] + 1) * sizeof(uint32_t));
                    memcpy(c_path_vec + c_path_num * READ_MAX_LEN, c_path_vec + it_path * READ_MAX_LEN, pos_path_vec[it_path * READ_MAX_LEN] + 1);
                }
                c_path_num += 1;
            }
            else{
                for(uint32_t it_base = 0; it_base < pos_path_vec[it_path * READ_MAX_LEN]; it_base ++){
                    sequence_tmp[pos_path_vec[it_path * READ_MAX_LEN + it_base + 1]] = c_path_vec[it_path * READ_MAX_LEN + it_base + 1];
                }
                if(index_largest_modified <= read_length + arg.max_extension - arg.kmer_length){
                    extend_amount = arg.kmer_length - (read_length - index_largest_modified);
                }
                else{
                    extend_amount = arg.max_extension;
                }
                bool extension = false;
				/*
 				char kmer_initial[36];
				char kmer_ext[READ_MAX_LEN]
				*/
                char *kmer_initial = funcm.cbuffer + (MAX_CANDIDATE_PATHS + 3) * READ_MAX_LEN;
                char *kmer_ext = funcm.cbuffer + (MAX_CANDIDATE_PATHS + 4) * READ_MAX_LEN;
                memcpy(kmer_ext, sequence_tmp + read_length -arg.kmer_length + 1, arg.kmer_length - 1);
                
                //uint32_t rec_index[extend_amount];
                uint32_t *rec_index = funcm.ubuffer + (MAX_CANDIDATE_PATHS + 5) * READ_MAX_LEN;
                char *path_content = funcm.cbuffer + (MAX_CANDIDATE_PATHS + 5) * READ_MAX_LEN;
                //char path_content[extend_amount * 4];
                uint32_t index_pos = 0;
                uint32_t num_extend = 0;
                uint32_t flag = 1;
                while(index_pos > 0 || flag == 1){
                    flag = 0;
                    memcpy(kmer_initial, kmer_ext + index_pos, arg.kmer_length - 1);
                   
                    //kmer_initial[arg.kmer_length - 1] = '0'; 
                    //path_query_text(kmer_initial, bit_vector, hash_seed, arg.kmer_length, arg.num_hash_func, arg.bit_vector_width, 1);
                    uint32_t k_ind = 0;
                    for(uint32_t it_alter = A; it_alter <= T; it_alter ++){
                        kmer_initial[arg.kmer_length - 1] = NEOCLEOTIDE[it_alter];
                       // if(kmer_initial[funcm.kmer_length + it_alter])
                        if(query_text(kmer_initial, bit_vector, hash_seed, arg.kmer_length, arg.num_hash_func, arg.bit_vector_width))
                            path_content[(index_pos << 2) + k_ind ++] = NEOCLEOTIDE[it_alter];
                    }
                    if(k_ind < 4)
                        path_content[(index_pos << 2) + k_ind] = 0;
                    if(num_extend + 1 == extend_amount && k_ind > 0){
                        extension = true;
                        break;
                    }else if(k_ind){
                        num_extend += 1;
                        rec_index[index_pos] = 0;
                        kmer_ext[arg.kmer_length - 1 + index_pos] = path_content[(index_pos << 2)];
                        rec_index[index_pos] = 1;
                        index_pos += 1;
                    } else{
                        while(index_pos > 0){
                            if(rec_index[index_pos - 1] < 4  && path_content[((index_pos - 1) << 2 ) + rec_index[index_pos - 1]]){
                                kmer_ext[arg.kmer_length - 1 + index_pos - 1] = path_content[((index_pos - 1) << 2) + rec_index[index_pos - 1]];
                                rec_index[index_pos - 1] += 1;
                                break;
                            }else{
                                num_extend -= 1;
                                index_pos -=1;
                            }
                        }

                    }
                }
                if(extension == true){
                    if(c_path_num != it_path){
                        memcpy(pos_path_vec + c_path_num * READ_MAX_LEN, pos_path_vec + it_path * READ_MAX_LEN, (pos_path_vec[it_path * READ_MAX_LEN] + 1) * sizeof(uint32_t));
                        memcpy(c_path_vec + c_path_num * READ_MAX_LEN, c_path_vec + it_path * READ_MAX_LEN, pos_path_vec[it_path * READ_MAX_LEN] + 1);
                    }
                    c_path_num += 1;
                }
            }
        }
    }   
    c_path_size = c_path_num;
   
    c_path_num = 0;
    for(uint32_t it_candidate = 0; it_candidate < c_path_size; it_candidate ++){
        bool real_modified = false;
        for(uint32_t it_mod_base = 0; it_mod_base < pos_path_vec[it_candidate * READ_MAX_LEN]; it_mod_base ++){
            if(org_sequence[pos_path_vec[it_candidate * READ_MAX_LEN + it_mod_base + 1]] != c_path_vec[it_candidate * READ_MAX_LEN + it_mod_base + 1]){
                real_modified = true;
            }
        }
        if(real_modified){
            if(it_candidate != c_path_num){
                memcpy(pos_path_vec + c_path_num * READ_MAX_LEN, pos_path_vec + it_candidate* READ_MAX_LEN, (pos_path_vec[it_candidate * READ_MAX_LEN] + 1) * sizeof(uint32_t));
                memcpy(c_path_vec + c_path_num * READ_MAX_LEN, c_path_vec + it_candidate * READ_MAX_LEN, pos_path_vec[it_candidate * READ_MAX_LEN] + 1);
            }
            c_path_num += 1; 
        }
    }

	//printf("cpath size %d cpath %d\n", c_path_size, c_path_num);
    c_path_size = c_path_num;

  /* 
   for(uint32_t ii = 0; ii < c_path_size; ii ++){
        for(uint32_t jj = 0; jj < pos_path_vec[ii * READ_MAX_LEN]; jj ++){
            printf("%d %c\n", pos_path_vec[ii * READ_MAX_LEN + jj + 1], c_path_vec[ii * READ_MAX_LEN + jj + 1]);
        }
        printf("****\n");
   }
   */
   // printf("******\n"); 
    
    
    if(all_solid_wo_modification == true){
    
    } else if(c_path_size > 1){
        uint32_t it_path;
        uint32_t it_path_1st;
        uint32_t it_path_2nd;
        
        uint32_t qs_1st = INIT_MIN_QS;
        uint32_t qs_2nd = INIT_MIN_QS;
        
        for(it_path = 0; it_path < c_path_size; it_path ++){
            uint32_t sum_qs = 0;
            for(uint32_t it_mod = 0; it_mod < pos_path_vec[it_path * READ_MAX_LEN]; it_mod ++){
                if(sequence[pos_path_vec[it_path * READ_MAX_LEN + it_mod + 1]] != c_path_vec[it_path * READ_MAX_LEN + it_mod + 1]){ 
                    sum_qs += (uint32_t)quality_score[pos_path_vec[it_path * READ_MAX_LEN + it_mod + 1]] - arg.quality_score_offset;
                }
            }
            if(sum_qs <= qs_1st){
                qs_2nd = qs_1st;
                qs_1st = sum_qs;
                it_path_2nd = it_path_1st;
                it_path_1st = it_path;
            }else if(sum_qs <= qs_2nd){
                qs_2nd = sum_qs;
                it_path_2nd = it_path;
            }
        }

        bool too_many_correction = false;
       
		
        if(read_length - pos_path_vec[it_path_1st * READ_MAX_LEN + 1] >= min_check_length && pos_path_vec[it_path_1st * READ_MAX_LEN] > MAX_MODIFICATION){
            for(uint32_t it_mod_base = 0; it_mod_base < pos_path_vec[it_path_1st * READ_MAX_LEN] - MAX_MODIFICATION; it_mod_base ++){
                if(pos_path_vec[it_path_1st * READ_MAX_LEN + it_mod_base + 1 + MAX_MODIFICATION] - pos_path_vec[it_path_1st * READ_MAX_LEN + it_mod_base + 1] < min_check_length){
                    if(read_length - pos_path_vec[it_path_1st * READ_MAX_LEN + it_mod_base + 1] <= max_trimmed_bases){
                        trim_3_end = read_length - pos_path_vec[it_path_1st * READ_MAX_LEN + it_mod_base + 1];
                    }else if(pos_path_vec[it_path_1st * READ_MAX_LEN + it_mod_base + 1] <= max_trimmed_bases){ 
                        trim_5_end = pos_path_vec[it_path_1st * READ_MAX_LEN + it_mod_base + 1] + 1;
                    }
                    too_many_correction = true;
                    break;
                }
            }
        }
        
       	/* 
        printf("trim 3 end %d \n", trim_3_end);
        printf("trim_5_end %d \n", trim_5_end);
        printf("too many corrections %d \n", too_many_correction);
        printf("qs 1st %d \n", qs_1st);
        printf("qs_2nd %d \n", qs_2nd);
        for(uint32_t ii = 0; ii < pos_path_vec[it_path_1st * READ_MAX_LEN]; ii ++){
            printf("%d %c\n", pos_path_vec[it_path_1st * READ_MAX_LEN + ii + 1], c_path_vec[it_path_1st * READ_MAX_LEN + ii + 1]);
        }
        printf("*****\n");
       */
        cpu_modify_sequence(org_sequence, sequence, quality_score, sequence_modification, trim_5_end, trim_3_end, 
         read_length, num_corrected_errors_local, pos_path_vec, c_path_vec, c_path_size, arg, funcm, it_path_1st, 
         qs_1st, qs_2nd, too_many_correction, max_trimmed_bases, 0); 
       	/*	
        for(uint32_t ii = 0; ii < read_length; ii ++){
            printf("%c", sequence[ii]);
        }
        printf("\n");
        */
    }
    else if(c_path_size == 1){
        bool too_many_correction = false;
        if(read_length - pos_path_vec[0 + 1] >= min_check_length && pos_path_vec[0] > MAX_MODIFICATION){
           for(uint32_t it_mod_base = 0; it_mod_base < pos_path_vec[0] - MAX_MODIFICATION; it_mod_base ++){
                if(pos_path_vec[it_mod_base + 1 + MAX_MODIFICATION] - pos_path_vec[it_mod_base + 1] < min_check_length){
                    if(read_length - pos_path_vec[it_mod_base + 1] <= max_trimmed_bases){
                        trim_3_end = read_length - pos_path_vec[it_mod_base + 1];
                        
                        if(it_mod_base > 0){
                            for(uint32_t it_base = 0; it_base < (it_mod_base - 1); it_base ++){
                                sequence_modification[pos_path_vec[it_base + 1]] = c_path_vec[it_base + 1];
                                sequence[pos_path_vec[it_base + 1]] = c_path_vec[it_base + 1];
                                num_corrected_errors_local ++;
                            }
                        }
                    }else if(pos_path_vec[it_mod_base + 1] + 1 <= max_trimmed_bases){
                        trim_5_end = pos_path_vec[it_mod_base + 1] + 1;
                        if(it_mod_base < pos_path_vec[0] -1 ){
                            for(uint32_t it_base = it_mod_base + 1; it_base < pos_path_vec[0]; it_base ++){
                                sequence_modification[pos_path_vec[it_base + 1]] = c_path_vec[it_base + 1];
                                sequence[pos_path_vec[it_base + 1]] = c_path_vec[it_base + 1];
                                num_corrected_errors_local ++;
                            }
                        }
                    }
                    too_many_correction = true;
                    break;
                }
           }
        }
        
    /*
        for(uint32_t ii = 0; ii < pos_path_vec[0]; ii ++){
            printf("%d %c\n", pos_path_vec[0 + ii + 1], c_path_vec[0 + ii + 1]);
        }
        printf("*****\n");
        for(uint32_t ii = 0; ii < read_length; ii ++){
            printf("%c", sequence[ii]);
        }
        printf("\n");
    */
        if(too_many_correction == false){
            for(uint32_t it_base = 0; it_base < pos_path_vec[0]; it_base ++){
                sequence_modification[pos_path_vec[it_base + 1]] = c_path_vec[it_base + 1];
                sequence[pos_path_vec[it_base + 1]] = c_path_vec[it_base + 1];
                num_corrected_errors_local ++;
            }
        } 
        
    }
    
}

 void cpu_correct_errors_first_kmer(const char* sequence, const char* quality_score, char* sequence_modification, uint32_t* pos_path_vec, char* c_path_vec, 
        uint32_t& c_path_size, uint8_t* bit_vector, uint32_t* hash_seed, init_args& arg, func_mem& funcm){
        //char first_kmer[36];
        //uint32_t low_qs_indexes[36];
        char *first_kmer = funcm.cbuffer + (MAX_CANDIDATE_PATHS + 1 + 1) * READ_MAX_LEN;
        uint32_t *low_qs_indexes = funcm.ubuffer + (MAX_CANDIDATE_PATHS + 1 + 4) * READ_MAX_LEN;
        uint32_t low_qs_indexes_ind = 0;

        //char new_kmer[36];
        //uint32_t c_path_size = 0;
        
        uint32_t rec_index[3];
        char path_content[3 * 4];
        int index_pos = 0;
        uint32_t num_extend = 0;
        uint32_t flag = 1;
        memcpy(first_kmer, sequence, arg.kmer_length);
        for(uint32_t it_bases = 0; it_bases < arg.kmer_length; it_bases ++){
            if((uint32_t)quality_score[it_bases] - arg.quality_score_offset < arg.quality_score_cutoff){
                low_qs_indexes[low_qs_indexes_ind ++] = it_bases;   
            }
        }
        //printf("low_qs_indexes == %d \n", low_qs_indexes_ind);
        if(low_qs_indexes_ind == 0){
            //printf("low qs index == %d \n", low_qs_indexes_ind);
            if(query_text(first_kmer, bit_vector, &hash_seed[0], arg.kmer_length, arg.num_hash_func, arg.bit_vector_width) == true){
                //save a path
                pos_path_vec[c_path_size * 4] = 0;
                c_path_size += 1;
            }
            else{
                for(uint32_t i = 0; i < arg.kmer_length; i ++){
                    low_qs_indexes[i] = i;
                }
                for(int i = 0; i < arg.kmer_length; i ++){
                    for(int j = i; j > 0; j --){
                        if(quality_score[low_qs_indexes[j]] < quality_score[low_qs_indexes[j - 1]]){
                            uint32_t t = low_qs_indexes[j];
                            low_qs_indexes[j] = low_qs_indexes[j - 1];
                            low_qs_indexes[j - 1] = t;
                        }else{
                            break;
                        }
                    }
                }
               /* 
                for(uint32_t ii = 0; ii < arg.kmer_length; ii ++){
                    printf("%d %x\n", low_qs_indexes[ii], quality_score[low_qs_indexes[ii]]);
                }
                printf("*********\n"); 
            */
               for(uint32_t it_alter1 = A; it_alter1 <= T; it_alter1 ++){
                    first_kmer[low_qs_indexes[0]] = NEOCLEOTIDE[it_alter1];
                    for(uint32_t it_alter2 = A; it_alter2 <= T; it_alter2 ++){
                        first_kmer[low_qs_indexes[1]] = NEOCLEOTIDE[it_alter2];
                        for(uint32_t it_alter3 = A; it_alter3 <= T; it_alter3 ++ ){
                            first_kmer[low_qs_indexes[2]] = NEOCLEOTIDE[it_alter3];
                            if(query_text(first_kmer, bit_vector, &hash_seed[0], arg.kmer_length, arg.num_hash_func, arg.bit_vector_width) == true){
                                //save path
                                pos_path_vec[c_path_size * 4] = 0;
                                if(sequence[low_qs_indexes[0]] != NEOCLEOTIDE[it_alter1]){
                                    pos_path_vec[c_path_size * 4 + pos_path_vec[c_path_size * 4] + 1] = low_qs_indexes[0]; 
                                    c_path_vec[c_path_size * 4 + pos_path_vec[c_path_size * 4] + 1] = NEOCLEOTIDE[it_alter1];
                                    pos_path_vec[c_path_size * 4] += 1;
                                }
                                
                                if(sequence[low_qs_indexes[1]] != NEOCLEOTIDE[it_alter2]){
                                    pos_path_vec[c_path_size * 4 + pos_path_vec[c_path_size * 4] + 1] = low_qs_indexes[1]; 
                                    c_path_vec[c_path_size * 4 + pos_path_vec[c_path_size * 4] + 1] = NEOCLEOTIDE[it_alter2];
                                    pos_path_vec[c_path_size * 4] += 1;
                                }
                                
                                if(sequence[low_qs_indexes[2]] != NEOCLEOTIDE[it_alter3]){
                                    pos_path_vec[c_path_size * 4 + pos_path_vec[c_path_size * 4] + 1] = low_qs_indexes[2]; 
                                    c_path_vec[c_path_size * 4 + pos_path_vec[c_path_size * 4] + 1] = NEOCLEOTIDE[it_alter3];
                                    pos_path_vec[c_path_size * 4] += 1;
                                }
                                c_path_size += 1;
                            }
                        }
                    }
                } 
                first_kmer[low_qs_indexes[0]] = sequence[low_qs_indexes[0]];
                first_kmer[low_qs_indexes[1]] = sequence[low_qs_indexes[1]];
                first_kmer[low_qs_indexes[2]] = sequence[low_qs_indexes[2]];
                for(uint32_t it_bases = 0; it_bases < arg.kmer_length; it_bases ++){
                    char c = first_kmer[it_bases];
                    for(uint32_t it_alter = A; it_alter <= T; it_alter ++){
                        if(first_kmer[it_bases] != NEOCLEOTIDE[it_alter]){
                            first_kmer[it_bases] = NEOCLEOTIDE[it_alter];
                            if(query_text(first_kmer, bit_vector, &hash_seed[0], arg.kmer_length, arg.num_hash_func, arg.bit_vector_width) == true){
                                //save_path
                                pos_path_vec[c_path_size * 4] = 1;
                                pos_path_vec[c_path_size * 4 + 1] = it_bases;
                                c_path_vec[c_path_size * 4 + 1] = NEOCLEOTIDE[it_alter];
                                c_path_size += 1;
                            }
                            first_kmer[it_bases] = c;
                        }
                    }
                }                    
            }
        } else if(low_qs_indexes_ind > 0){
            if(low_qs_indexes_ind <= MAX_LOW_QS_BASES){
                switch(low_qs_indexes_ind){
                    case 3:
                    for(uint32_t it_alter1 = A; it_alter1 <= T; it_alter1 ++){
                        first_kmer[low_qs_indexes[0]] = NEOCLEOTIDE[it_alter1];
                        for(uint32_t it_alter2 = A; it_alter2 <= T; it_alter2 ++){
                            first_kmer[low_qs_indexes[1]] = NEOCLEOTIDE[it_alter2];
                            for(uint32_t it_alter3 = A; it_alter3 <= T; it_alter3 ++) {
                                first_kmer[low_qs_indexes[2]] = NEOCLEOTIDE[it_alter3];
                                if(query_text(first_kmer, bit_vector, &hash_seed[0], arg.kmer_length, arg.num_hash_func, arg.bit_vector_width) == true){
                                    //save path                                
                                    pos_path_vec[c_path_size * 4] = 0;
                                    if(sequence[low_qs_indexes[0]] != NEOCLEOTIDE[it_alter1]){
                                        pos_path_vec[c_path_size * 4 + pos_path_vec[c_path_size * 4] + 1] = low_qs_indexes[0]; 
                                        c_path_vec[c_path_size * 4 + pos_path_vec[c_path_size * 4] + 1] = NEOCLEOTIDE[it_alter1];
                                        pos_path_vec[c_path_size * 4] += 1;
                                    }
                                
                                    if(sequence[low_qs_indexes[1]] != NEOCLEOTIDE[it_alter2]){
                                        pos_path_vec[c_path_size * 4 + pos_path_vec[c_path_size * 4] + 1] = low_qs_indexes[1]; 
                                        c_path_vec[c_path_size * 4 + pos_path_vec[c_path_size * 4] + 1] = NEOCLEOTIDE[it_alter2];
                                        pos_path_vec[c_path_size * 4] += 1;
                                    }
                                
                                    if(sequence[low_qs_indexes[2]] != NEOCLEOTIDE[it_alter3]){
                                        pos_path_vec[c_path_size * 4 + pos_path_vec[c_path_size * 4] + 1] = low_qs_indexes[2]; 
                                        c_path_vec[c_path_size * 4 + pos_path_vec[c_path_size * 4] + 1] = NEOCLEOTIDE[it_alter3];
                                        pos_path_vec[c_path_size * 4] += 1;
                                    }
                                    c_path_size += 1;
                                }
                            }
                        }
                    }
                    break;
                    case 2:
                    for(uint32_t it_alter1 = A; it_alter1 <= T; it_alter1 ++){
                        first_kmer[low_qs_indexes[0]] = NEOCLEOTIDE[it_alter1];
                        for(uint32_t it_alter2 = A; it_alter2 <= T; it_alter2 ++){
                            first_kmer[low_qs_indexes[1]] = NEOCLEOTIDE[it_alter2];
                            if(query_text(first_kmer, bit_vector, &hash_seed[0], arg.kmer_length, arg.num_hash_func, arg.bit_vector_width) == true){
                                //save path                                
                                pos_path_vec[c_path_size * 4] = 0;
                                if(sequence[low_qs_indexes[0]] != NEOCLEOTIDE[it_alter1]){
                                    pos_path_vec[c_path_size * 4 + pos_path_vec[c_path_size * 4] + 1] = low_qs_indexes[0]; 
                                    c_path_vec[c_path_size * 4 + pos_path_vec[c_path_size * 4] + 1] = NEOCLEOTIDE[it_alter1];
                                    pos_path_vec[c_path_size * 4] += 1;
                                }
                                
                                if(sequence[low_qs_indexes[1]] != NEOCLEOTIDE[it_alter2]){
                                    pos_path_vec[c_path_size * 4 + pos_path_vec[c_path_size * 4] + 1] = low_qs_indexes[1]; 
                                    c_path_vec[c_path_size * 4 + pos_path_vec[c_path_size * 4] + 1] = NEOCLEOTIDE[it_alter2];
                                    pos_path_vec[c_path_size * 4] += 1;
                                }
                                
                                c_path_size += 1;

                            } 
                        }
                    }
                    break;
                    case 1:
                    for(uint32_t it_alter1 = A; it_alter1 <= T; it_alter1 ++){
                        first_kmer[low_qs_indexes[0]] = NEOCLEOTIDE[it_alter1];
                        if(query_text(first_kmer, bit_vector, &hash_seed[0], arg.kmer_length, arg.num_hash_func, arg.bit_vector_width) == true){
                                //save path
                                pos_path_vec[c_path_size * 4] = 0;
                                if(sequence[low_qs_indexes[0]] != NEOCLEOTIDE[it_alter1]){
                                    pos_path_vec[c_path_size * 4 + pos_path_vec[c_path_size * 4] + 1] = low_qs_indexes[0]; 
                                    c_path_vec[c_path_size * 4 + pos_path_vec[c_path_size * 4] + 1] = NEOCLEOTIDE[it_alter1];
                                    pos_path_vec[c_path_size * 4] += 1;
                                }
                                c_path_size += 1;
                        }
                    }
                    break;
                }
               
                for(uint32_t i = 0; i < low_qs_indexes_ind; i ++){ 
                    first_kmer[low_qs_indexes[i]] = sequence[low_qs_indexes[i]];
                } 

                if(c_path_size == 0){
                    if(query_text(first_kmer, bit_vector, &hash_seed[0], arg.kmer_length, arg.num_hash_func, arg.bit_vector_width) == true){
                        // save nothing path
                        pos_path_vec[c_path_size * 4] = 0;
                        c_path_size += 1;
                    }
                    else{
                        for(uint32_t it_bases = 0; it_bases < arg.kmer_length; it_bases ++){
                            char c = first_kmer[it_bases];
                            for(uint32_t it_alter = A; it_alter <= T; it_alter ++){
                                if(first_kmer[it_bases] != NEOCLEOTIDE[it_alter]){
                                    first_kmer[it_bases] = NEOCLEOTIDE[it_alter];
                                    if(query_text(first_kmer, bit_vector, &hash_seed[0], arg.kmer_length, arg.num_hash_func, arg.bit_vector_width) == true){
                                        pos_path_vec[c_path_size * 4] = 1;
                                        pos_path_vec[c_path_size * 4 + 1] = it_bases;
                                        c_path_vec[c_path_size * 4 + 1] = NEOCLEOTIDE[it_alter];
                                        c_path_size += 1;
                                    }
                                    first_kmer[it_bases] = c;
                                }
                            }
                        }                           
                    }
                }
            }else{
                  for(int i = 0; i < low_qs_indexes_ind; i ++){
                    for(int j = i; j > 0; j --){
                        if(quality_score[low_qs_indexes[j]] < quality_score[low_qs_indexes[j - 1]]){
                            uint32_t t = low_qs_indexes[j];
                            low_qs_indexes[j] = low_qs_indexes[j - 1];
                            low_qs_indexes[j - 1] = t;
                        }else{
                            break;
                        }
                    }
                } 
                for(uint32_t it_alter1 = A; it_alter1 <= T; it_alter1 ++){
                    first_kmer[low_qs_indexes[0]] = NEOCLEOTIDE[it_alter1];
                    for(uint32_t it_alter2 = A; it_alter2 <= T; it_alter2 ++){
                        first_kmer[low_qs_indexes[1]] = NEOCLEOTIDE[it_alter2];
                        for(uint32_t it_alter3 = A; it_alter3 <= T; it_alter3 ++){
                            first_kmer[low_qs_indexes[2]] = NEOCLEOTIDE[it_alter3];
                            if(query_text(first_kmer, bit_vector, &hash_seed[0], arg.kmer_length, arg.num_hash_func, arg.bit_vector_width) == true){
                                //save path                                    
                                    pos_path_vec[c_path_size * 4] = 0;
                                    if(sequence[low_qs_indexes[0]] != NEOCLEOTIDE[it_alter1]){
                                        pos_path_vec[c_path_size * 4 + pos_path_vec[c_path_size * 4] + 1] = low_qs_indexes[0]; 
                                        c_path_vec[c_path_size * 4 + pos_path_vec[c_path_size * 4] + 1] = NEOCLEOTIDE[it_alter1];
                                        pos_path_vec[c_path_size * 4] += 1;
                                    }
                                
                                    if(sequence[low_qs_indexes[1]] != NEOCLEOTIDE[it_alter2]){
                                        pos_path_vec[c_path_size * 4 + pos_path_vec[c_path_size * 4] + 1] = low_qs_indexes[1]; 
                                        c_path_vec[c_path_size * 4 + pos_path_vec[c_path_size * 4] + 1] = NEOCLEOTIDE[it_alter2];
                                        pos_path_vec[c_path_size * 4] += 1;
                                    }
                                
                                    if(sequence[low_qs_indexes[2]] != NEOCLEOTIDE[it_alter3]){
                                        pos_path_vec[c_path_size * 4 + pos_path_vec[c_path_size * 4] + 1] = low_qs_indexes[2]; 
                                        c_path_vec[c_path_size * 4 + pos_path_vec[c_path_size * 4] + 1] = NEOCLEOTIDE[it_alter3];
                                        pos_path_vec[c_path_size * 4] += 1;
                                    }
                                    c_path_size += 1;
                            }
                        }
                    }
                }

                first_kmer[low_qs_indexes[0]] = sequence[low_qs_indexes[0]];
                first_kmer[low_qs_indexes[1]] = sequence[low_qs_indexes[1]];
                first_kmer[low_qs_indexes[2]] = sequence[low_qs_indexes[2]];
                for(uint32_t it_bases = 0; it_bases < arg.kmer_length; it_bases ++){
                    char c = first_kmer[it_bases];
                    for(uint32_t it_alter = A; it_alter <= T; it_alter ++){
                        if(first_kmer[it_bases] != NEOCLEOTIDE[it_alter]){
                            first_kmer[it_bases] = NEOCLEOTIDE[it_alter];
                            if(query_text(first_kmer, bit_vector, &hash_seed[0], arg.kmer_length, arg.num_hash_func, arg.bit_vector_width) == true){
                                //save_path
                                pos_path_vec[c_path_size * 4] = 1;
                                pos_path_vec[c_path_size * 4 + 1] = it_bases;
                                c_path_vec[c_path_size * 4 + 1] = NEOCLEOTIDE[it_alter];
                                c_path_size += 1;
                            }
                            first_kmer[it_bases] = c;
                        }
                    }
                } 
                if(query_text(first_kmer, bit_vector, &hash_seed[0], arg.kmer_length, arg.num_hash_func, arg.bit_vector_width) == true){
                    //save nothing path
                    pos_path_vec[c_path_size * 4] = 0;
                    c_path_size += 1;;
                }                
            }
        }
    
}

 void cpu_solid_first_kmer(uint32_t* pos_path_vec, char* c_path_vec, const char* sequence, bool& extension_success, uint8_t* bit_vector, uint32_t* hash_seed, 
       init_args& arg, func_mem& funcm){

    uint32_t index_smallest_modified = pos_path_vec[0 + 1];
    //char first_kmer[36];
    uint32_t extend_amount;
    if(index_smallest_modified >= arg.kmer_length - arg.max_extension - 1){
        extend_amount = arg.kmer_length - index_smallest_modified - 1;
    }
    else{
        extend_amount = arg.max_extension;
    }
    
    //char kmer_ext[READ_MAX_LEN];
    char *kmer_ext = funcm.cbuffer + (MAX_CANDIDATE_PATHS + 1 + 1) * READ_MAX_LEN;
    memcpy(kmer_ext  + extend_amount, sequence, arg.kmer_length);
    for(uint32_t it_bases = 0; it_bases < pos_path_vec[0]; it_bases ++){
        kmer_ext[extend_amount + pos_path_vec[it_bases + 1]] = c_path_vec[it_bases + 1];
    }
    bool extension = false;
    //char kmer_initial[36];
    char *kmer_initial = funcm.cbuffer + (MAX_CANDIDATE_PATHS + 1 + 2) * READ_MAX_LEN;
    //memcpy(kmer_ext + extend_amount, first_kmer, arg.kmer_length - 1);
                
    //uint32_t rec_index[extend_amount];
    //char path_content[extend_amount * 4];
    
    uint32_t *rec_index = funcm.ubuffer + (MAX_CANDIDATE_PATHS + 1 + 4) * READ_MAX_LEN;
    char *path_content = funcm.cbuffer + (MAX_CANDIDATE_PATHS + 1 + 3) * READ_MAX_LEN;

    uint32_t index_pos = 0;
    uint32_t num_extend = 0;
    uint32_t flag = 1;
    while(index_pos > 0 || flag == 1){
        flag = 0;
        memcpy(kmer_initial + 1, kmer_ext + extend_amount - num_extend, arg.kmer_length - 1);
        /*
        for(uint32_t ii = 0; ii < index_pos + arg.kmer_length; ii ++){
            printf("%c", kmer_ext[ii + extend_amount - num_extend]);
        }
        printf("******\n");
        */
        //path_query_text(kmer_initial, bit_vector, hash_seed, arg.kmer_length, arg.num_hash_func, arg.bit_vector_width, 1);
        uint32_t k_ind = 0;
        for(uint32_t it_alter = A; it_alter <= T; it_alter ++){
            kmer_initial[0] = NEOCLEOTIDE[it_alter];
            //if(kmer_initial[funcm.kmer_length + it_alter])
            if(query_text(kmer_initial, bit_vector, hash_seed, arg.kmer_length, arg.num_hash_func, arg.bit_vector_width) == true)
                path_content[(index_pos << 2) + k_ind ++] = NEOCLEOTIDE[it_alter];
        }
        if(k_ind < 4)
            path_content[(index_pos << 2) + k_ind] = 0;
        if(num_extend + 1 == extend_amount && k_ind > 0){
            extension_success = true;
            break;
        }else if(k_ind){
            num_extend += 1;
            //rec_index[index_pos] = 0;
            kmer_ext[extend_amount - num_extend] = path_content[(index_pos << 2)];
            rec_index[index_pos] = 1;
            index_pos += 1;
        } else{
            while(index_pos > 0){
                if(rec_index[index_pos - 1] < 4  && path_content[((index_pos - 1)<< 2 ) + rec_index[index_pos - 1]]){
                    kmer_ext[extend_amount - index_pos] = path_content[((index_pos - 1) << 2) + rec_index[index_pos - 1]];
                    rec_index[index_pos - 1] += 1;
                    break;
                }else{
                    num_extend -= 1;
                    index_pos -=1;
                }
            }

        }
    }
                
    
}
 inline void cpu_extend_a_kmer_3_prime_end(const char* kmer, const char* sequence, uint32_t index_kmer, uint32_t* pos_path_vec, 
        char* c_path_vec, uint32_t& c_path_size, uint32_t* pos_path_vec_first, char* c_path_vec_first,const uint32_t& org_boundary, const char* quality_score, const uint8_t* bit_vector, 
        uint32_t* hash_seed, bool& run_exploration, uint32_t read_length, init_args& arg, func_mem& funcm){

    //char path_content[4 * READ_MAX_LEN];
    int index_pos = 0;
    //char* kmer_ext = kmer + 1;
    //char  kmer_new[36];
    //uint32_t rec_index[READ_MAX_LEN];
    //char kmer_ext[READ_MAX_LEN];

	char *path_content = funcm.cbuffer + (MAX_CANDIDATE_PATHS + 1 + 2) * READ_MAX_LEN;
	char *kmer_ext = funcm.cbuffer + (MAX_CANDIDATE_PATHS + 1 + 6) * READ_MAX_LEN;
	uint32_t *rec_index = funcm.ubuffer + (MAX_CANDIDATE_PATHS + 1 + 4) * READ_MAX_LEN;
    memcpy(kmer_ext + index_pos, kmer + 1, arg.kmer_length - 1);
    char *kmer_new = funcm.cbuffer + (MAX_CANDIDATE_PATHS + 1 + 7) * READ_MAX_LEN;
    uint32_t flag = 1;
    uint32_t index_start = index_kmer;
    //printf("hhhhhh%d\n", index_kmer);
    while((index_pos > 0) || (flag == 1)){
        flag = 0;
        memcpy(kmer_new, kmer_ext + index_pos, arg.kmer_length - 1);
       	/* 
        for(uint32_t ii = 0; ii < arg.kmer_length + index_pos - 1; ii ++){
            printf("%c", kmer_ext[ii]);
        }
        printf("\n");
        
        
        for(uint32_t ii = 0; ii < arg.kmer_length - 1; ii ++){
            printf("%c", kmer_new[ii]);
        }
        printf("\n");
        */
        /*
        kmer_new[arg.kmer_length] = sequence[index_kmer + arg.kmer_length];
        path_query_text(kmer_new, bit_vector, hash_seed, arg.kmer_length, arg.num_hash_func, arg.bit_vector_width, 1);
        uint32_t k_ind = 0;
        for(unsigned short int it_alter = A; it_alter <= T; it_alter ++){ 
            if(kmer_new[funcm.kmer_length + it_alter]){
                path_content[(index_pos << 2) + k_ind ++] = NEOCLEOTIDE[it_alter];
            }
                       
        }
        */
        /*
        for(int ii = 0; ii < index_pos; ii ++){
            for(uint32_t jj = 0; jj < 4; jj ++){
                if(path_content[(ii << 2) + jj] > 0){
                    printf("%c", path_content[(ii << 2) + jj]);
                }else{
                    break;
                }
            }
            printf(" ");
        }
        printf("\n");
        for(uint32_t ii = 0; ii < index_pos + arg.kmer_length - 1; ii ++){
            printf("%c", kmer_ext[ii]);
        }
        printf("\n");
        for(int ii = 0; ii < index_pos; ii ++){
            printf("%d ", rec_index[ii]);
        }
        printf("\n");
		*/
	
       	/*
	   	char cc, ccc;
        scanf("%c%c", &cc, &ccc);
        printf("cc %c \n", cc); 
        */
        uint32_t k_ind = 0;
        if ((index_kmer == (org_boundary - 1)) || (((uint32_t)quality_score[index_kmer + arg.kmer_length] - arg.quality_score_offset) <= arg.extremely_low_quality_score)){
            for(unsigned short int it_alter = A; it_alter <= T; it_alter ++){ 
                kmer_new[arg.kmer_length - 1] = NEOCLEOTIDE[it_alter];
                //if(kmer_new[funcm.kmer_length + it_alter]){
                if(query_text(kmer_new, bit_vector, &hash_seed[0], arg.kmer_length, arg.num_hash_func, arg.bit_vector_width) == true){
                    path_content[(index_pos << 2) + k_ind ++] = NEOCLEOTIDE[it_alter];      
                }
            }                
        }else{
            kmer_new[arg.kmer_length - 1] = sequence[index_kmer + arg.kmer_length];
            if(query_text(kmer_new, bit_vector, &hash_seed[0], arg.kmer_length, arg.num_hash_func, arg.bit_vector_width) == true){ 
                path_content[(index_pos << 2) + k_ind ++] = kmer_new[arg.kmer_length - 1];      
            }else{
                 for(unsigned short int it_alter = A; it_alter <= T; it_alter ++){ 
                    kmer_new[arg.kmer_length - 1] = NEOCLEOTIDE[it_alter];
                    //if(kmer_new[funcm.kmer_length + it_alter]){
                    if(sequence[index_kmer + arg.kmer_length] != NEOCLEOTIDE[it_alter])
                    if(query_text(kmer_new, bit_vector, &hash_seed[0], arg.kmer_length, arg.num_hash_func, arg.bit_vector_width) == true){
                        path_content[(index_pos << 2) + k_ind ++] = NEOCLEOTIDE[it_alter];      
                    }
                }      
            } 
        }

        if(k_ind < 4)
            path_content[(index_pos << 2) + k_ind] = 0;
        if(index_kmer + 1 == read_length - arg.kmer_length && k_ind > 0){
            for(uint32_t it = 0; it < k_ind; it ++){
                //save path
                kmer_ext[arg.kmer_length - 1 + index_pos] = path_content[(index_pos << 2) + it];
                pos_path_vec[c_path_size * READ_MAX_LEN] = 0;
				for(uint32_t nt = 0; nt < pos_path_vec_first[0]; nt ++){
                	pos_path_vec[c_path_size * READ_MAX_LEN + nt + 1] = pos_path_vec_first[nt + 1];
                 	c_path_vec[c_path_size * READ_MAX_LEN + nt + 1] = c_path_vec_first[nt + 1];
                    pos_path_vec[c_path_size * READ_MAX_LEN] += 1;
				}
                /*
                for(uint32_t ii = 0; ii <= index_pos + 1; ii ++){
                    printf("%c", kmer_ext[arg.kmer_length -2 + ii]);
                }
                printf("\n");
                for(uint32_t ii = 0; ii <= index_pos + 1; ii ++){
                    printf("%c", sequence[arg.kmer_length - 1 + index_start + ii]);
                }
                printf("\n");
                //printf("*********\n");
                */
                for(uint32_t ct = 0; ct <= index_pos + 1; ct ++){
                    if(kmer_ext[arg.kmer_length - 2 + ct] != sequence[arg.kmer_length - 1 + index_start + ct]){
                        pos_path_vec[c_path_size * READ_MAX_LEN + pos_path_vec[c_path_size * READ_MAX_LEN] + 1] = index_start + ct + arg.kmer_length - 1 ;
                        c_path_vec[c_path_size * READ_MAX_LEN + pos_path_vec[c_path_size * READ_MAX_LEN] + 1] = kmer_ext[arg.kmer_length - 2 + ct];
                        pos_path_vec[c_path_size * READ_MAX_LEN] += 1;
                    }
                }

   /*             for(uint32_t ii = 0; ii < pos_path_vec[c_path_size * READ_MAX_LEN]; ii ++){
                    printf("%c", c_path_vec[c_path_size * READ_MAX_LEN + ii  + 1]);
                }
                printf("\n******\n");
    */            c_path_size += 1;        
                if(c_path_size > MAX_CANDIDATE_PATHS){
                    //退出寻找路径
                    index_pos = -1;       
                    run_exploration = false;
                    break;
                } 
            }  
            
            while(index_pos > 0){
                if(rec_index[index_pos - 1] < 4 && path_content[((index_pos - 1) << 2 ) + rec_index[index_pos - 1]] > 0){
                    kmer_ext[arg.kmer_length - 1 + index_pos - 1] = path_content[((index_pos - 1) << 2) + rec_index[index_pos - 1]];
                    rec_index[index_pos - 1] += 1;
                    break;
                }else{
                    index_kmer -= 1;
                    index_pos -=1;
                }
            }

        }
        else if (k_ind){
            index_kmer += 1;
            rec_index[index_pos] = 0;
            kmer_ext[arg.kmer_length - 1 + index_pos] = path_content[(index_pos << 2) + rec_index[index_pos]];
            rec_index[index_pos] += 1;
            index_pos += 1;              
        }else{
            while(index_pos > 0){
                if(rec_index[index_pos - 1] < 4 && path_content[((index_pos - 1) << 2) + rec_index[index_pos - 1]] > 0){
                    kmer_ext[arg.kmer_length - 1 + index_pos - 1] = path_content[((index_pos - 1) << 2) + rec_index[index_pos - 1]];
                    rec_index[index_pos - 1] += 1;
                    break;
                }else{
                    index_kmer -= 1;
                    index_pos -=1;
                }
            }
		}
		/*
        for(uint32_t ii = 0; ii < arg.kmer_length; ii ++){
            printf("%x ", kmer_new[ii]);
        }
        printf("\n"); 
        for(uint32_t ii = 0; ii < arg.kmer_length; ii ++){
            printf("%x ", kmer[ii]);
        }
        printf("\n"); 
        printf("%d\n", index_kmer);
        printf("%d\n", index_pos);
        printf("*****\n");
		*/
    }
   
}
 void cpu_extend_first_kmer_to_right(const char* sequence, const char* quality_score, uint32_t* pos_path_vec_first, char* c_path_vec_first,
        uint32_t* pos_path_vec, char* c_path_vec, uint32_t& c_path_size, uint32_t read_length, uint8_t* bit_vector, uint32_t* hash_seed, 
        init_args& arg, func_mem& funcm, bool& run_exploration){
    
    //char first_kmer[36];
    /*
    char *first_kmer = funcm.cbuffer + (MAX_CANDIDATE_PATHS * 2 + 2 + 1) * READ_MAX_LEN;
    memcpy(first_kmer, sequence, arg.kmer_length + 1);
    for(uint32_t it_base = 0; it_base < pos_path_vec_first[0]; it_base ++){
        first_kmer[pos_path_vec_first[it_base + 1]] = c_path_vec_first[it_base + 1];
    }
    */

    //char second_kmer[36];
	//printf("hhhhhhhhh\n");	
	
	//uint32_t *pos_path_vec_tmp = funcm.ubuffer + 4 * READ_MAX_LEN;
    //char *c_path_vec_tmp = funcm.cbuffer + READ_MAX_LEN;
    //uint32_t c_path_size_tmp = 0;  
    char *second_kmer = funcm.cbuffer + (MAX_CANDIDATE_PATHS  + 1 + 1) * READ_MAX_LEN;
    memcpy(second_kmer, sequence + 1, arg.kmer_length);
    for(uint32_t it_base = 0; it_base < pos_path_vec_first[0]; it_base ++){
        if(pos_path_vec_first[it_base + 1]  >= 1){
            //printf(" mm %d  %d\n", pos_path_vec_first[it_base + 1] - 1, it_base);
            second_kmer[pos_path_vec_first[it_base + 1] - 1] = c_path_vec_first[it_base + 1];
        }
    }
   // memcpy(second_kmer, first_kmer + 1, arg.kmer_length);

    if(query_text(second_kmer, bit_vector, &hash_seed[0], arg.kmer_length, arg.num_hash_func, arg.bit_vector_width) == true){
        if(read_length - arg.kmer_length == 1){
         //save a path
         memcpy(pos_path_vec + c_path_size * READ_MAX_LEN, pos_path_vec_first, 4 * sizeof(uint32_t));
         memcpy(c_path_vec + c_path_size * READ_MAX_LEN, c_path_vec_first, 4);  
         c_path_size += 1;
        } else if(read_length - arg.kmer_length > 1){  
			//printf("nnnnnnnnn111\n");
            cpu_extend_a_kmer_3_prime_end(
                        second_kmer,
                        sequence,
                        1,
                        pos_path_vec,
                        c_path_vec,
                        c_path_size,
						pos_path_vec_first,
						c_path_vec_first,
                        read_length,
                        quality_score,
                        bit_vector,
                        hash_seed,
                        run_exploration,
                        read_length,
                        arg,
                        funcm
                        );
			//printf("ggggggggggg111\n");
             
        }
    }
    else{
        //path_query_text(second_kmer, bit_vector, hash_seed, arg.kmer_length, arg.num_hash_func, arg.bit_vector_width, 1);
        for(unsigned short int it_alter = A; it_alter <= T; it_alter ++){
            if(sequence[arg.kmer_length] != NEOCLEOTIDE[it_alter]){
                second_kmer[arg.kmer_length - 1] = NEOCLEOTIDE[it_alter];
				/*
				printf("33333\n");
				for(uint32_t ii = 0; ii < 31; ii ++){
					printf("%x", second_kmer[ii]);
				}
				printf("\n");
    			for(uint32_t it_base = 0; it_base < pos_path_vec_first[0]; it_base ++){
        			printf("%d %x\n", pos_path_vec_first[it_base + 1], c_path_vec_first[it_base + 1]);
    			}
				*/
				if(run_exploration == true)
                if(query_text(second_kmer, bit_vector, hash_seed, arg.kmer_length, arg.num_hash_func, arg.bit_vector_width) == true){
                    //generate new path
					
         			memcpy(pos_path_vec + c_path_size * READ_MAX_LEN, pos_path_vec_first, 4 * sizeof(uint32_t));
         			memcpy(c_path_vec + c_path_size * READ_MAX_LEN, c_path_vec_first, 4);  
                    //pos_path_vec_tmp[c_path_size_tmp * READ_MAX_LEN] = 1;
                    pos_path_vec[c_path_size * READ_MAX_LEN + pos_path_vec[c_path_size * READ_MAX_LEN] + 1] = arg.kmer_length;
                    c_path_vec[c_path_size * READ_MAX_LEN + pos_path_vec[c_path_size * READ_MAX_LEN] + 1] = NEOCLEOTIDE[it_alter];
					pos_path_vec[c_path_size * READ_MAX_LEN] += 1;
                    if(read_length - arg.kmer_length == 1){
                        //save new path
                        c_path_size += 1;
                    }else if(read_length - arg.kmer_length > 1){
                        //printf("nnnnnnnnnnnnnnnnn222\n");
                        cpu_extend_a_kmer_3_prime_end(
                                    second_kmer,
                                    sequence,
                                    1,
                                    pos_path_vec,
                                    c_path_vec,
                                    c_path_size,
									pos_path_vec_first,
									c_path_vec_first,
                                    read_length,
                                    quality_score,
                                    bit_vector,
                                    hash_seed,
                                    run_exploration,
                                    read_length,
                                    arg,
                                    funcm
                                    );
                        //printf("ggggggggggg222\n");
                    }
                }
				/*
				printf("4444444\n");
   				for(uint32_t ii = 0; ii < 31; ii ++){
					printf("%x", second_kmer[ii]);
				}
				printf("\n");
				*/
            }
        }
    }
	/*	
	for(uint32_t ii = 0; ii < c_path_size; ii ++){
		for(uint32_t jj = 0; jj < pos_path_vec[ii * READ_MAX_LEN]; jj ++){
			printf("%d %c\n", pos_path_vec[ii * READ_MAX_LEN + jj + 1], c_path_vec[ii * READ_MAX_LEN + jj + 1]);
		}
		printf("****\n");
	}
	printf("mmmmmmm\n");
	*/
    if(run_exploration == false){
        //clear candidate_path_vector
        c_path_size = 0;
    }
	//printf("kkkkkkkkkk\n");
	uint32_t c_path_num = 0;
    for(uint32_t it_path = 0; it_path < c_path_size; it_path ++){
        uint32_t index_largest_modified = pos_path_vec[it_path * READ_MAX_LEN + pos_path_vec[it_path * READ_MAX_LEN]];
        uint32_t extend_amount;

       if(index_largest_modified <= read_length - arg.kmer_length){
            //save path
			if(c_path_num != it_path){
            	memcpy(pos_path_vec + c_path_num * READ_MAX_LEN, pos_path_vec + it_path * READ_MAX_LEN, (pos_path_vec[it_path * READ_MAX_LEN] + 1) * sizeof(uint32_t));
            	memcpy(c_path_vec + c_path_num * READ_MAX_LEN, c_path_vec + it_path * READ_MAX_LEN, pos_path_vec[it_path * READ_MAX_LEN] + 1);
			}
           c_path_num += 1;
       } 
       else {
            //char sequence_tmp[READ_MAX_LEN];
            char *sequence_tmp = funcm.cbuffer + (MAX_CANDIDATE_PATHS  + 1 + 1) * READ_MAX_LEN;
            memcpy(sequence_tmp, sequence, read_length);
			/*
			for(uint32_t it_base = 0; it_base < pos_path_vec[it_path * READ_MAX_LEN]; it_base ++){
				sequence_tmp[pos_path_vec[it_path * READ_MAX_LEN + it_base + 1]] = c_path_vec[it_path * READ_MAX_LEN + it_base + 1];
			}

            uint32_t index_largest_modified = pos_path_vec[it_path * READ_MAX_LEN + pos_path_vec[it_path * READ_MAX_LEN]];
            uint32_t extend_amount;
            if(index_largest_modified <= read_length - arg.kmer_length){
                if(c_path_num != it_path){
                    memcpy(pos_path_vec + c_path_num * READ_MAX_LEN, pos_path_vec + it_path * READ_MAX_LEN, (pos_path_vec[it_path * READ_MAX_LEN] + 1) * sizeof(uint32_t));
                    memcpy(c_path_vec + c_path_num * READ_MAX_LEN, c_path_vec + it_path * READ_MAX_LEN, pos_path_vec[it_path * READ_MAX_LEN] + 1);
                }
                c_path_num += 1;
            }
            else{
			*/
                for(uint32_t it_base = 0; it_base < pos_path_vec[it_path * READ_MAX_LEN]; it_base ++){
                    sequence_tmp[pos_path_vec[it_path * READ_MAX_LEN + it_base + 1]] = c_path_vec[it_path * READ_MAX_LEN + it_base + 1];
                }
                if(index_largest_modified <= read_length + arg.max_extension - arg.kmer_length){
                    extend_amount = arg.kmer_length - (read_length - index_largest_modified);
                }
                else{
                    extend_amount = arg.max_extension;
                }
                bool extension = false;
                //char kmer_initial[36];
                //char kmer_ext[READ_MAX_LEN];
                char *kmer_initial = funcm.cbuffer + (MAX_CANDIDATE_PATHS  + 1 + 2) * READ_MAX_LEN;
                char *kmer_ext = funcm.cbuffer + (MAX_CANDIDATE_PATHS  + 1 + 3) * READ_MAX_LEN;
                memcpy(kmer_ext, sequence_tmp + read_length -arg.kmer_length + 1, arg.kmer_length - 1);
                /* 
                uint32_t rec_index[extend_amount];
                char path_content[extend_amount * 4];
                */
                uint32_t *rec_index = funcm.ubuffer + (MAX_CANDIDATE_PATHS  + 1 + 4) * READ_MAX_LEN;
                char *path_content = funcm.cbuffer + (MAX_CANDIDATE_PATHS  + 1 + 4) * READ_MAX_LEN;
                uint32_t index_pos = 0;
                uint32_t num_extend = 0;
                uint32_t flag = 1;
                while(index_pos > 0 || flag == 1){
                    flag = 0;
                    memcpy(kmer_initial, kmer_ext + index_pos, arg.kmer_length - 1);
                   
                    //kmer_initial[arg.kmer_length - 1] = '0'; 
                    //path_query_text(kmer_initial, bit_vector, hash_seed, arg.kmer_length, arg.num_hash_func, arg.bit_vector_width, 1);
                    uint32_t k_ind = 0;
                    for(uint32_t it_alter = A; it_alter <= T; it_alter ++){
                        kmer_initial[arg.kmer_length - 1] = NEOCLEOTIDE[it_alter];
                       // if(kmer_initial[funcm.kmer_length + it_alter])
                        if(query_text(kmer_initial, bit_vector, hash_seed, arg.kmer_length, arg.num_hash_func, arg.bit_vector_width))
                            path_content[(index_pos << 2) + k_ind ++] = NEOCLEOTIDE[it_alter];
                    }
                    if(k_ind < 4)
                        path_content[(index_pos << 2) + k_ind] = 0;
                    if(num_extend + 1 == extend_amount && k_ind > 0){
                        extension = true;
                        break;
                    }else if(k_ind){
                        num_extend += 1;
                        rec_index[index_pos] = 0;
                        kmer_ext[arg.kmer_length - 1 + index_pos] = path_content[(index_pos << 2)];
                        rec_index[index_pos] = 1;
                        index_pos += 1;
                    } else{
                        while(index_pos > 0){
                            if(rec_index[index_pos - 1] < 4  && path_content[((index_pos - 1) << 2 ) + rec_index[index_pos - 1]]){
                                kmer_ext[arg.kmer_length - 1 + index_pos - 1] = path_content[((index_pos - 1) << 2) + rec_index[index_pos - 1]];
                                rec_index[index_pos - 1] += 1;
                                break;
                            }else{
                                num_extend -= 1;
                                index_pos -=1;
                            }
                        }

                    }
                }
                if(extension == true){
					if(c_path_num != it_path){
                   		memcpy(pos_path_vec + c_path_num * READ_MAX_LEN, pos_path_vec + it_path * READ_MAX_LEN, (pos_path_vec[it_path * READ_MAX_LEN] + 1) * sizeof(uint32_t));
                    	memcpy(c_path_vec + c_path_num * READ_MAX_LEN, c_path_vec + it_path * READ_MAX_LEN, pos_path_vec[it_path * READ_MAX_LEN] + 1);
					}
                    c_path_num += 1;
                }
           // }

            }
    }
	c_path_size = c_path_num;
	/*
	for(uint32_t ii = 0; ii < c_path_size; ii ++){
		for(uint32_t jj = 0; jj < pos_path_vec[ii * READ_MAX_LEN]; jj ++){
			printf("%d %c\n", pos_path_vec[ii * READ_MAX_LEN + jj + 1], c_path_vec[ii * READ_MAX_LEN + jj + 1]);
		}
		printf("****\n");
	}
	printf("kkkkkkkkk\n");
	*/
	/*
	for(uint32_t ii = 0; ii < c_path_size; ii ++){
		for(uint32_t jj = 0; jj < pos_path_vec[ii * READ_MAX_LEN]; jj ++){
			printf("%d %c\n", pos_path_vec[ii * READ_MAX_LEN + jj + 1], c_path_vec[ii * READ_MAX_LEN + jj + 1]);
		}
		printf("****\n");
	}
	printf("nnnnnnnnn\n");
	*/
	
}
 void cpu_correct_errors(char* read_content, char* read_quality, char* sequence_modification,
        uint32_t* hash_seed, uint8_t* bit_vector, init_args& arg,
        uint32_t& num_corrected_errors_local,uint32_t& trim_5_end, uint32_t& trim_3_end,
        uint32_t max_trimmed_bases,uint32_t min_check_length, uint32_t read_length, int ssss, func_mem& funcm) {
    	
	uint32_t *solid_regions1 = funcm.ubuffer;
	uint32_t *solid_regions2 = funcm.ubuffer + READ_MAX_LEN;
	uint32_t *solid_regions1_org = funcm.ubuffer + READ_MAX_LEN * 2;
	uint32_t *solid_regions2_org = funcm.ubuffer + READ_MAX_LEN * 3;
    uint32_t solid_regions_org_ind = 0;
	uint32_t num_kmers =read_length- arg.kmer_length + 1;
	uint32_t solid_regions_ind = 0;
	int regions_ind_prev = 0;
    bool is_solid_kmer_prev = false;
	uint32_t *regions_text = funcm.ubuffer + READ_MAX_LEN * 4;
    char *sequence_modified = funcm.cbuffer;
    memcpy(sequence_modified, read_content, read_length);   

        
    char *kmer = funcm.cbuffer + READ_MAX_LEN;
	//regions_query_text(sequence_modified, bit_vector, hash_seed, arg.kmer_length, arg.num_hash_func, arg.bit_vector_width, read_length, regions_text);	
    //printf("arg num hash func %d\n", arg.num_hash_func);
    
    //for(int i = 0; i < arg.num_hash_func; i ++)
    //    printf("%d\n", hash_seed[i]);    
    //找到连续的正确的kmer区域
	for(uint32_t it_kmer = 0; it_kmer < num_kmers; it_kmer ++){
        memcpy(kmer, read_content + it_kmer, arg.kmer_length);
		//if(regions_text[2 * READ_MAX_LEN + it_kmer] == 1){
        if(query_text(kmer, bit_vector, hash_seed, arg.kmer_length, arg.num_hash_func, arg.bit_vector_width) == true){ 
			if(is_solid_kmer_prev == false){
				solid_regions1[solid_regions_ind] = it_kmer;
                solid_regions1_org[solid_regions_org_ind] = it_kmer;
				regions_ind_prev = 1;
                is_solid_kmer_prev = true;
			}

		}else{
			if(is_solid_kmer_prev == true){
				solid_regions2[solid_regions_ind ++] = it_kmer - 1;
                solid_regions2_org[solid_regions_org_ind ++] = it_kmer - 1;
				regions_ind_prev = 0;	
                is_solid_kmer_prev = false;
			}
		}
	}

	//last solid regions
	if(is_solid_kmer_prev == true){
		solid_regions2[solid_regions_ind ++] = num_kmers - 1;	
        solid_regions2_org[solid_regions_org_ind ++] = num_kmers - 1;
	}
    
    //0-1
    //根据read的quality调整区域
     
    
    uint32_t *solid_regions1_tmp = funcm.ubuffer + 4 * READ_MAX_LEN;
    uint32_t *solid_regions2_tmp = funcm.ubuffer + 5 * READ_MAX_LEN;
	uint32_t *solid_regions1_org_tmp = funcm.ubuffer + 6 * READ_MAX_LEN;
	uint32_t *solid_regions2_org_tmp = funcm.ubuffer + 7 * READ_MAX_LEN;
    uint32_t flag = 0;
    uint32_t solid_regions_tmp_ind = 0;
    uint32_t solid_regions_org_tmp_ind = 0;
    if(solid_regions_ind > 0){
    if ((solid_regions_ind != 1) || (solid_regions1[0] != 0) || (solid_regions2[0] != (read_length - arg.kmer_length))) {
    if(solid_regions_ind > 1){
        for(uint32_t it = 0; it < solid_regions_ind - 1; it ++){
            if(solid_regions1[it + 1] - solid_regions2[it] < arg.kmer_length){
                flag = 1;         
            }
        }
        if(flag == 1){
            for(uint32_t it_str = 0; it_str < solid_regions_ind; it_str ++){
                uint32_t low_quality_base = 0;
                //uint32_t low_quality_base_prev = solid_regions1[it_str] - 1;
                uint32_t low_quality_base_prev = 0;
                for(uint32_t it_base = solid_regions1[it_str]; it_base < solid_regions2[it_str] + arg.kmer_length; it_base ++){
                    if((uint32_t)read_quality[it_base] - arg.quality_score_offset < arg.quality_score_cutoff){                              
                        low_quality_base +=1;
                        if(low_quality_base == 1){   
                            if(it_base >=  solid_regions1[it_str] + arg.kmer_length){
                                solid_regions1_tmp[solid_regions_tmp_ind] = solid_regions1[it_str];
                                solid_regions2_tmp[solid_regions_tmp_ind ++ ] = it_base - arg.kmer_length;
                                solid_regions1_org_tmp[solid_regions_org_tmp_ind] = solid_regions1_org[it_str];
                                solid_regions2_org_tmp[solid_regions_org_tmp_ind ++] = solid_regions2_org[it_str];
                            }
                        }
                        else if(it_base >  low_quality_base_prev + arg.kmer_length){
                            solid_regions1_tmp[solid_regions_tmp_ind] = low_quality_base_prev + 1;
                            solid_regions2_tmp[solid_regions_tmp_ind ++ ] = it_base - arg.kmer_length;
                            solid_regions1_org_tmp[solid_regions_org_tmp_ind] = solid_regions1_org[it_str];
                            solid_regions2_org_tmp[solid_regions_org_tmp_ind ++] = solid_regions2_org[it_str];
                        }
                        low_quality_base_prev = it_base;
                    }          
                }
                if(low_quality_base > 0){
                    if(solid_regions2[it_str] >= low_quality_base_prev + arg.kmer_length){
                            solid_regions1_tmp[solid_regions_tmp_ind] = low_quality_base_prev + arg.kmer_length;
                            solid_regions2_tmp[solid_regions_tmp_ind ++ ] = solid_regions2[it_str]; 
                            solid_regions1_org_tmp[solid_regions_org_tmp_ind] = solid_regions1_org[it_str];
                            solid_regions2_org_tmp[solid_regions_org_tmp_ind ++] = solid_regions2_org[it_str];     
                    }
                }else{
                    solid_regions1_tmp[solid_regions_tmp_ind] = solid_regions1[it_str];
                    solid_regions2_tmp[solid_regions_tmp_ind ++] = solid_regions2[it_str];
                    solid_regions1_org_tmp[solid_regions_org_tmp_ind] = solid_regions1_org[it_str];
                    solid_regions2_org_tmp[solid_regions_org_tmp_ind ++] = solid_regions2_org[it_str];
                }
            }
            memcpy(solid_regions1, solid_regions1_tmp, solid_regions_tmp_ind * sizeof(uint32_t));
            memcpy(solid_regions2, solid_regions2_tmp, solid_regions_tmp_ind * sizeof(uint32_t));
            solid_regions_ind = solid_regions_tmp_ind;
            memcpy(solid_regions1_org, solid_regions1_org_tmp, solid_regions_org_tmp_ind * sizeof(uint32_t));
            memcpy(solid_regions2_org, solid_regions2_org_tmp, solid_regions_org_tmp_ind * sizeof(uint32_t));
            solid_regions_org_ind = solid_regions_org_tmp_ind;
        
            
        
        }
    } else if(solid_regions_ind == 1){
        uint32_t low_quality_base = 0;
        uint32_t low_quality_base_prev = 0;
        //int low_quality_base_prev = solid_regions1[0] - 1;
        for(uint32_t it_base = solid_regions1[0]; it_base < solid_regions2[0] + arg.kmer_length; it_base ++){
            //printf("%d\n", (uint32_t)read_quality[it_base] - arg.quality_score_offset);
            if((uint32_t)read_quality[it_base] - arg.quality_score_offset < arg.quality_score_cutoff){                              
                low_quality_base +=1;
                if(low_quality_base == 1){   
                    if(it_base  >= solid_regions1[0] + arg.kmer_length + MIN_SOLID_LENGTH - 1){
                        solid_regions1_tmp[solid_regions_tmp_ind] = solid_regions1[0];
                        solid_regions2_tmp[solid_regions_tmp_ind ++ ] = it_base - arg.kmer_length;   
                        solid_regions1_org_tmp[solid_regions_org_tmp_ind] = solid_regions1_org[0];
                        solid_regions2_org_tmp[solid_regions_org_tmp_ind ++] = solid_regions2_org[0];                    
                    }
                }
                else if(it_base  >= low_quality_base_prev + arg.kmer_length + MIN_SOLID_LENGTH){   
                    solid_regions1_tmp[solid_regions_tmp_ind] = low_quality_base_prev + 1;
                    solid_regions2_tmp[solid_regions_tmp_ind ++ ] = it_base - arg.kmer_length;   
                    solid_regions1_org_tmp[solid_regions_org_tmp_ind] = solid_regions1_org[0];
                    solid_regions2_org_tmp[solid_regions_org_tmp_ind ++] = solid_regions2_org[0];                    
                }
                low_quality_base_prev = it_base;
            }          
        }
        if(solid_regions_tmp_ind){ 
            memcpy(solid_regions1, solid_regions1_tmp, solid_regions_tmp_ind * sizeof(uint32_t));
            memcpy(solid_regions2, solid_regions2_tmp, solid_regions_tmp_ind * sizeof(uint32_t));
            solid_regions_ind = solid_regions_tmp_ind;
            memcpy(solid_regions1_org, solid_regions1_org_tmp, solid_regions_org_tmp_ind * sizeof(uint32_t));
            memcpy(solid_regions2_org, solid_regions2_org_tmp, solid_regions_org_tmp_ind * sizeof(uint32_t));
            solid_regions_org_ind = solid_regions_org_tmp_ind;
        } 
    }
    }
    }
    
   //0-2
    //去掉区域长度过短的区域
    if(solid_regions_ind){
        uint32_t ind_tmp = 0;
        for(uint32_t it = 0; it < solid_regions_ind; it ++){
            if(solid_regions2[it] - solid_regions1[it] + 1 >= MIN_SOLID_LENGTH){
                if(ind_tmp < it){
                    solid_regions1[ind_tmp] = solid_regions1[it];
                    solid_regions2[ind_tmp] = solid_regions2[it];
                    solid_regions1_org[ind_tmp] = solid_regions1_org[it];
                    solid_regions2_org[ind_tmp ++] = solid_regions2_org[it];
                }else{
                    ind_tmp ++;
                }
            }
        }
        solid_regions_ind = ind_tmp;
        solid_regions_org_ind = ind_tmp;
    } 

    
    // 0-3
    // remove short non-solid regions
    if(solid_regions_ind > 1){
        uint32_t ind_tmp = 1;
        for(uint32_t it = 1; it < solid_regions_ind; it ++){
            if(solid_regions1[it] - solid_regions2[it -1] < MIN_SOLID_LENGTH + 1){
                if(ind_tmp - 1 < it)
                    solid_regions2[ind_tmp - 1] = solid_regions2[it];
            }else{
                if(ind_tmp < it){
                    solid_regions1[ind_tmp] = solid_regions1[it];
                    solid_regions2[ind_tmp] = solid_regions2[it];
                    solid_regions1_org[ind_tmp] = solid_regions1_org[it];
                    solid_regions2_org[ind_tmp ++] = solid_regions2_org[it];
                }else{
                    ind_tmp += 1;
                }       
            }
        }
        solid_regions_ind = ind_tmp;
        solid_regions_org_ind = ind_tmp;
    }
    
    
    //0-4
    //reduce size of solid regions
    
    if(solid_regions_ind > 1){
        for(uint32_t it = 1; it < solid_regions_ind; it ++){
            if(solid_regions1[it] - solid_regions2[it - 1] - 1 < arg.kmer_length && solid_regions1[it] - solid_regions2[it - 1] - 1 >= arg.kmer_length - FP_SUSPECT_LENGTH){
                if(solid_regions2[it] - solid_regions1[it] + 1 > FP_SUSPECT_LENGTH){
                    solid_regions1[it] += FP_SUSPECT_LENGTH;
                }               
            
                if(solid_regions2[it - 1] - solid_regions1[it - 1] + 1 > FP_SUSPECT_LENGTH){
                    solid_regions2[it -1] -= FP_SUSPECT_LENGTH;
                }
            }
        }
    }
    
    //0-5
    //remove a solid region that makes a non-solid reiong shorter than k
    
    if(solid_regions_ind == 2){
        if(solid_regions1[0] == 0){
            if(solid_regions1[1] - solid_regions2[0] < arg.kmer_length + 1){
                solid_regions_ind = 1;
            }
        }else if(solid_regions2[1] == read_length - arg.kmer_length){
            if(solid_regions1[1] - solid_regions2[0] < arg.kmer_length + 1 ){
                if(solid_regions2[1] - solid_regions1[1] + 1 >= read_length * 0.1){
                    if(solid_regions2[0] - solid_regions1[0] + 1 < read_length * 0.1){
                        solid_regions1[0] = solid_regions1[1];
                        solid_regions2[0] = solid_regions2[1];
                        solid_regions_ind = 1;
                    }
                }
            }
        }
    }
    
    //0-6
    //check the quality scores of right side of each solid k-mer region
    
    if(solid_regions_ind > 0){
        for(uint32_t it_str = 0; it_str < solid_regions_ind - 1; it_str ++){
            if(solid_regions2[it_str] - solid_regions1[it_str] > SOLID_REGION_ADJUST_RANGE){
                for(uint32_t it_adjust = solid_regions2[it_str]; it_adjust > solid_regions2[it_str] - SOLID_REGION_ADJUST_RANGE; it_adjust --){
                    if(((uint32_t)read_quality[it_adjust + arg.kmer_length - 1] - arg.quality_score_offset < arg.quality_score_cutoff) 
                            || ((uint32_t)read_quality[it_adjust] - arg.quality_score_offset < arg.quality_score_cutoff)){
                       solid_regions2[it_str] = it_adjust - 1;
                        break; 
                    } 
                }
            }
        }    

        uint32_t index_solid_region = solid_regions_ind - 1;
        if(solid_regions2[index_solid_region] < read_length - arg.kmer_length){
             if(solid_regions2[index_solid_region] - solid_regions1[index_solid_region] > SOLID_REGION_ADJUST_RANGE){
                for(uint32_t it_adjust = solid_regions2[index_solid_region]; it_adjust > solid_regions2[index_solid_region] - SOLID_REGION_ADJUST_RANGE; it_adjust --){
                    if(read_quality[it_adjust + arg.kmer_length - 1] - arg.quality_score_offset < arg.quality_score_cutoff){
                        solid_regions2[index_solid_region] = it_adjust - 1;
                        break;
                    }
                }
             }
        }
        if(solid_regions1[0] > 0){
            if(solid_regions2[0] - solid_regions1[0] > SOLID_REGION_ADJUST_RANGE){
                for(uint32_t it_adjust = solid_regions1[0]; it_adjust < solid_regions1[0] + SOLID_REGION_ADJUST_RANGE; it_adjust ++){
                    if((uint32_t)read_quality[it_adjust + arg.kmer_length - 1] - arg.quality_score_offset < arg.quality_score_cutoff){
                        solid_regions1[0] = it_adjust + 1;
                        break; 
                    }
                } 
            } 
        }
    }
   
    //0-7
    //check whether a non-solid region < k still exists
    
    bool short_non_solid_region = false;
    if(solid_regions_ind > 1){
        for(uint32_t it_str = 1; it_str < solid_regions_ind - 1; it_str ++){
            if(solid_regions1[it_str] - solid_regions2[it_str - 1] <= arg.kmer_length){
                short_non_solid_region = true;
                break; 
            }
        }
    }
    /*
    for(int it_str = 0; it_str < solid_regions_ind;it_str ++){
        printf("%d %d\n", solid_regions1[it_str], solid_regions2[it_str]);
        printf("%d %d\n", solid_regions1_org[it_str], solid_regions2_org[it_str]);
    }
    */
      if(solid_regions_ind > 0 && short_non_solid_region == false){
        // STEP 1-1: Correct errors between solid regions
        if(solid_regions_ind > 1){
            for(uint32_t it_region = 1; it_region < solid_regions_ind; it_region ++){
                if((solid_regions1[it_region] - 1) - (solid_regions2[it_region - 1] + 1) + 1 >= arg.kmer_length){
                    if(solid_regions2[it_region - 1] + arg.kmer_length < read_length - trim_3_end && solid_regions1[it_region] > trim_5_end){
                       //printf("kk %d\n", ssss + 1);
					  	 
                        cpu_correct_errors_between_solid_regions(
                                read_content,
                                sequence_modified,
                                read_quality,
                                solid_regions1[it_region - 1],
                                solid_regions2[it_region - 1] + 1,
                                solid_regions1[it_region] - 1,
                                solid_regions2[it_region],
                                solid_regions2_org[it_region - 1] - 1,
                                solid_regions1_org[it_region] - 1,
                                sequence_modification,
                                trim_5_end,
                                trim_3_end,
                                solid_regions_ind,
                                read_length,
                                max_trimmed_bases,
                                min_check_length,
                                num_corrected_errors_local,
                                bit_vector,
                                hash_seed,
                                arg,
                                funcm  
                                );
								
                        //printf("kkkk%d\n", ssss + 1);
                        
                    }
                }
            }
        }
		//printf("%d %d\n", solid_regions2[solid_regions_ind - 1], trim_3_end);
	//	printf("ssss == %d\n", ssss);
        if(solid_regions_ind >= 1){
            if(solid_regions1[0] > 0){
                if(solid_regions1[0] > trim_5_end){
                    //printf("jj%d\n", ssss + 1);
						
                    cpu_correct_errors_5_prime_end(
                            read_content,
                            sequence_modified,
                            read_quality,
                            solid_regions1[0] - 1,
                            sequence_modification,
                            trim_5_end,
                            trim_3_end,
                            solid_regions1_org[0] - 1,
                            read_length,
                            max_trimmed_bases,
                            min_check_length,
                            num_corrected_errors_local,
                            bit_vector,
                            hash_seed,
                            arg,
                            funcm
                            );
                    
                    //printf("jjjj%d\n", ssss + 1);
                }
            }
        }
			//printf("%d %d\n", solid_regions2[solid_regions_ind - 1], trim_3_end);
        if(solid_regions_ind >= 1){
            if(solid_regions2[solid_regions_ind - 1] < read_length - arg.kmer_length){
                if(solid_regions2[solid_regions_ind - 1] + arg.kmer_length < read_length - trim_3_end){
                    //printf("mm%d\n", ssss + 1);
					
                    cpu_correct_errors_3_prime_end(
                            read_content,
                            sequence_modified,
                            read_quality,
                            solid_regions2[solid_regions_ind - 1] + 1,
                            sequence_modification,
                            trim_5_end,
                            trim_3_end,
                            solid_regions2_org[solid_regions_ind - 1] + 1,
                            read_length,
                            max_trimmed_bases,
                            min_check_length,
                            num_corrected_errors_local,
                            bit_vector,
                            hash_seed,
                            arg,
                            funcm
                            );
							
                    //printf("mmmmm%d\n", ssss + 1);
                }
            }
        }
	        if(num_corrected_errors_local == 0 && trim_5_end == 0 && trim_3_end == 0){
            trim_5_end = solid_regions1[0];
            for(uint32_t it_solid_short = 0; it_solid_short < solid_regions_ind; it_solid_short ++){
                if(trim_5_end + (read_length - arg.kmer_length - solid_regions2[it_solid_short]) <= max_trimmed_bases){
                    trim_3_end = read_length - arg.kmer_length - solid_regions2[it_solid_short];
                    break;
                }
            }
        }
       //printf("hhh%d\n", ssss + 1); 
    }
   else{
        //2 - 1
        //correct errors in the first kmer
            
        uint32_t *pos_path_vec_first = funcm.ubuffer + (2 * MAX_CANDIDATE_PATHS + 2 + 4) * READ_MAX_LEN;
        char *c_path_vec_first = funcm.cbuffer + (2 * MAX_CANDIDATE_PATHS + 2 + 1) * READ_MAX_LEN;
        uint32_t c_path_size_first = 0;
            //printf("mm %d\n", ssss + 1); 
        cpu_correct_errors_first_kmer(
                sequence_modified,
                read_quality,
                sequence_modification,
                pos_path_vec_first,
                c_path_vec_first,
                c_path_size_first,
                bit_vector,
                hash_seed,
                arg,
                funcm
                ); 
            
            //printf("******\n");
        
                   
        uint32_t c_path_num = 0;
        if(c_path_size_first > 0){
            for(uint32_t it_candidates = 0; it_candidates < c_path_size_first; it_candidates ++){
                if(pos_path_vec_first[it_candidates * 4] == 0){
                    if(c_path_num != it_candidates){
                        memcpy(pos_path_vec_first + c_path_num * 4, pos_path_vec_first + it_candidates * 4, 4 * sizeof(uint32_t));
                        memcpy(c_path_vec_first + c_path_num * 4, c_path_vec_first + it_candidates * 4, 4);
                    }
                    c_path_num += 1;
                }else{
                    bool really_modified = false;
                    for(uint32_t it_mod_bases = 0; it_mod_bases < pos_path_vec_first[it_candidates * 4]; it_mod_bases ++){
                        if(read_content[pos_path_vec_first[it_candidates * 4 + it_mod_bases + 1]] != c_path_vec_first[it_candidates * 4 + it_mod_bases + 1]){
                            really_modified = true;
                        }
                    }
                    if(really_modified == true){
                        if(pos_path_vec_first[it_candidates * 4 + 1] < arg.kmer_length - 1){
                            bool extension_success = false;
                            //printf("%d \n", ssss + 1);
                            cpu_solid_first_kmer(
                                    pos_path_vec_first + it_candidates * 4,
                                    c_path_vec_first + it_candidates * 4,
                                    sequence_modified,
                                    extension_success,
                                    bit_vector,
                                    hash_seed,
                                    arg,
                                    funcm
                                    );
                            //printf("*****\n");
                            if(extension_success == true){
                                if(c_path_num != it_candidates){
                                    memcpy(pos_path_vec_first + c_path_num * 4, pos_path_vec_first + it_candidates * 4, 4 * sizeof(uint32_t));
                                    memcpy(c_path_vec_first + c_path_num * 4, c_path_vec_first + it_candidates * 4, 4);
                                }
                                c_path_num += 1;
                            }
                        }
                        else{
                            if(c_path_num != it_candidates){
                                memcpy(pos_path_vec_first + c_path_num * 4, pos_path_vec_first + it_candidates * 4, 4 * sizeof(uint32_t));
                                memcpy(c_path_vec_first + c_path_num * 4, c_path_vec_first + it_candidates * 4, 4);
                            }
                            c_path_num += 1;
                        }
 

                    }
                } 
            }
        }
        
        c_path_size_first = c_path_num;

    		//printf("start\n");	
		uint32_t *pos_path_vec = funcm.ubuffer + 238 * READ_MAX_LEN;
        char *c_path_vec = funcm.cbuffer + 236 * READ_MAX_LEN;
        uint32_t c_path_size = 0;  
        if(c_path_size_first > 0){
            bool run_exploration = true;
            for(uint32_t it_candidates = 0; it_candidates < c_path_size_first; it_candidates ++){
				uint32_t *pos_path_vec_tmp = pos_path_vec + c_path_size * READ_MAX_LEN;
				char* c_path_vec_tmp = c_path_vec + c_path_size * READ_MAX_LEN;
				uint32_t c_path_size_tmp = 0;
                if(run_exploration == true){
                    cpu_extend_first_kmer_to_right(
                            sequence_modified,
                            read_quality,
                            pos_path_vec_first + it_candidates * 4,
                            c_path_vec_first + it_candidates * 4,
                            pos_path_vec_tmp,
                            c_path_vec_tmp,
                            c_path_size_tmp,
                            read_length,
                            bit_vector,
                            hash_seed,
                            arg,
                            funcm,
                            run_exploration
                            );
                }
				c_path_size += c_path_size_tmp;
            }
            if(run_exploration == false){
                c_path_size = 0;
            }
        }
		
              c_path_num = 0;
        for(uint32_t it_candidate = 0; it_candidate < c_path_size; it_candidate ++){
            bool really_modified = false;
            for(uint32_t it_mod_base = 0; it_mod_base < pos_path_vec[it_candidate * READ_MAX_LEN]; it_mod_base++){
                if(read_content[pos_path_vec[it_candidate * READ_MAX_LEN + it_mod_base + 1]] != c_path_vec[it_candidate * READ_MAX_LEN + it_mod_base + 1]){
                    really_modified = true;
                }
            }
            if(really_modified){
				if(c_path_num != it_candidate){
                	memcpy(pos_path_vec + c_path_num * READ_MAX_LEN, pos_path_vec + it_candidate * READ_MAX_LEN, (pos_path_vec[it_candidate * READ_MAX_LEN]  + 1) * sizeof(uint32_t));
                	memcpy(c_path_vec + c_path_num * READ_MAX_LEN, c_path_vec + it_candidate * READ_MAX_LEN, pos_path_vec[it_candidate * READ_MAX_LEN]  + 1);
				}
                c_path_num ++;
            }
        }
        c_path_size = c_path_num;
	        if(c_path_size > 1){
			
            uint32_t it_path;
            uint32_t it_path_1st;
            uint32_t it_path_2nd;
        
            uint32_t qs_1st = INIT_MIN_QS;
            uint32_t qs_2nd = INIT_MIN_QS;
        
            for(it_path = 0; it_path < c_path_size; it_path ++){
                uint32_t sum_qs = 0;
                for(uint32_t it_mod = 0; it_mod < pos_path_vec[it_path * READ_MAX_LEN]; it_mod ++){
                    if(read_content[pos_path_vec[it_path * READ_MAX_LEN + it_mod + 1]] != c_path_vec[it_path * READ_MAX_LEN + it_mod + 1]){ 
                        sum_qs +=read_quality[pos_path_vec[it_path * READ_MAX_LEN + it_mod + 1]] - arg.quality_score_offset;
                    }
                }
                if(sum_qs <= qs_1st){
                    qs_2nd = qs_1st;
                    qs_1st = sum_qs;
                    it_path_2nd = it_path_1st;
                    it_path_1st = it_path;
                }else if(sum_qs <= qs_2nd){
                    qs_2nd = sum_qs;
                    it_path_2nd = it_path;
                }
            }
            bool too_many_correction = false;
       

            if(read_length - pos_path_vec[it_path_1st * READ_MAX_LEN + 1] >= min_check_length && pos_path_vec[it_path_1st * READ_MAX_LEN] > MAX_MODIFICATION){
                for(uint32_t it_mod_base = 0; it_mod_base < pos_path_vec[it_path_1st * READ_MAX_LEN] - MAX_MODIFICATION; it_mod_base ++){
                    if(pos_path_vec[it_path_1st * READ_MAX_LEN + it_mod_base + 1 + MAX_MODIFICATION] - pos_path_vec[it_path_1st * READ_MAX_LEN + it_mod_base + 1]< min_check_length){
                        if(read_length - pos_path_vec[it_path_1st * READ_MAX_LEN + it_mod_base + 1] <= max_trimmed_bases){
                            trim_3_end = read_length - pos_path_vec[it_path_1st * READ_MAX_LEN + it_mod_base + 1];
                        }else if(pos_path_vec[it_path_1st * READ_MAX_LEN + it_mod_base + 1] + 1 <= max_trimmed_bases){
                            trim_5_end = pos_path_vec[it_path_1st * READ_MAX_LEN + it_mod_base + 1] + 1;
                        }
                        too_many_correction = true;
                        break;
                    }
                }
            }
	            cpu_modify_sequence(read_content, sequence_modified, read_quality, sequence_modification, trim_5_end, trim_3_end, 
                read_length, num_corrected_errors_local, pos_path_vec, c_path_vec, c_path_size, arg, funcm, 
                it_path_1st, qs_1st, qs_2nd, too_many_correction, max_trimmed_bases, 1); 
		
        } else if(c_path_size == 1){
            bool too_many_correction = false;
            if(read_length - pos_path_vec[0 + 1] >= min_check_length && pos_path_vec[0] > MAX_MODIFICATION){
                for(uint32_t it_mod_base = 0; it_mod_base < pos_path_vec[0] - MAX_MODIFICATION; it_mod_base ++){
                    if(pos_path_vec[it_mod_base + 1 + MAX_MODIFICATION] - pos_path_vec[it_mod_base + 1] < min_check_length){
                        if(read_length - pos_path_vec[it_mod_base + 1] <= max_trimmed_bases){
                            trim_3_end = read_length - pos_path_vec[it_mod_base + 1];
                        
                            if(it_mod_base > 0){
                                for(uint32_t it_base = 0; it_base < (it_mod_base - 1); it_base ++){
                                    sequence_modification[pos_path_vec[it_base + 1]] = c_path_vec[it_base + 1];
                                    sequence_modified[pos_path_vec[it_base + 1]] = c_path_vec[it_base + 1];
                                    num_corrected_errors_local ++;
                                }
                            }
                        }else if(pos_path_vec[it_mod_base + 1] + 1 <= max_trimmed_bases){
                            trim_5_end = pos_path_vec[it_mod_base + 1] + 1;
                            if(it_mod_base < pos_path_vec[0] -1 ){
                                for(uint32_t it_base = it_mod_base + 1; it_base < pos_path_vec[0]; it_base ++){
                                    sequence_modification[pos_path_vec[it_base + 1]] = c_path_vec[it_base + 1];
                                    sequence_modified[pos_path_vec[it_base + 1]] = c_path_vec[it_base + 1];
                                    num_corrected_errors_local ++;
                                }
                            }
                        }
                        too_many_correction = true;
                        break;
                    }
                }
            }
            if(too_many_correction == false){
                for(uint32_t it_base = 0; it_base < pos_path_vec[0]; it_base ++){
                    sequence_modification[pos_path_vec[it_base + 1]] = c_path_vec[it_base + 1];
                    sequence_modified[pos_path_vec[it_base + 1]] = c_path_vec[it_base + 1];
                    num_corrected_errors_local ++;
                }
            }
	        } else{
            if(solid_regions_ind > 0){
                trim_5_end = solid_regions1[0];
                for(uint32_t it_solid_short = 0; it_solid_short < solid_regions_ind; it_solid_short ++){
                    if(trim_5_end + read_length - arg.kmer_length - solid_regions2[it_solid_short] <= max_trimmed_bases){
                        trim_3_end = read_length - arg.kmer_length - solid_regions2[it_solid_short];
                        break; 
                    }
                }
            }else if(arg.kmer_length < max_trimmed_bases){
                trim_5_end = arg.kmer_length;
            }
        
        } 
    }
    
 

}


