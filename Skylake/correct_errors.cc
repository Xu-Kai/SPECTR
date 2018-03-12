/*
 * correct_errors.c
 *
 *  Created on: 2015年12月14日
 *      Author: xk
 */


#include "correct_errors.h"
#include "query_text.h"
#include "mic_vec_query_text.h"
#define BLOCK_READ_NUM 2000000
#define READ_MAX_LEN 256
//#define READ_LEN 101
#include "struct.h"

//----------------------------------------------------------------------
// generate_hash_seed
//----------------------------------------------------------------------
void generate_hash_seed(const uint32_t random_seed, uint32_t num_hash_func,
        unsigned int* hash_seed) {

    srand(random_seed);

    for (unsigned short int it_seed = 0; it_seed < num_hash_func; it_seed++) {
        hash_seed[it_seed] = static_cast<unsigned int>(rand());
    }
}

void read_fastq(FILE* pfile, uint32_t& read_num, char* read_name, char* read_plus,
        char* read_content, char* read_quality) {
    int buffer_size = 100 * 1024 * 1024; //100M
    char* buffer = (char*) malloc(buffer_size);
    int read_size = 1;
    uint32_t flag = 0;
    int iter = 0;
    read_num = 0;
    int pointer = 0;
    while (read_num < BLOCK_READ_NUM && read_size != 0) {
        read_size = fread(buffer, 1, buffer_size, pfile);
        //		printf("first char is %c %d\n", buffer[0], read_size);

        for (iter = 0; iter < read_size; iter++) {
            if (flag == 0) {
                read_name[read_num * READ_MAX_LEN + pointer++] = buffer[iter];
                if (buffer[iter] == '\n') {
                    pointer = 0;
                    flag = 1;
                }

            } else if (flag == 1) {
                read_content[read_num * READ_MAX_LEN + pointer++] =
                    buffer[iter];
                if (buffer[iter] == '\n') {
                    pointer = 0;
                    flag = 2;
                }
            } else if (flag == 2) {
                read_plus[read_num * READ_MAX_LEN  + pointer++] = buffer[iter];
                if (buffer[iter] == '\n') {
                    pointer = 0;
                    flag = 3;
                }
            } else if (flag == 3) {
                read_quality[read_num * READ_MAX_LEN + pointer++] =
                    buffer[iter];
                if (buffer[iter] == '\n') {
                    flag = 0;
                    pointer = 0;
                    read_num ++;
                    //					if(read_num == 1){
                    //						printf("%s \n", read_name);
                    //						printf("%s \n", read_content);
                    //						printf("%s \n", read_plus);
                    //						printf("%s \n", read_quality);
                    //					}
                    if (read_num == BLOCK_READ_NUM)
                        break;
                }
            }
        }
        if(read_num == BLOCK_READ_NUM)
            break;
    }

    fseek(pfile, iter - read_size + 1, SEEK_CUR);
    free(buffer);
}



//----------------------------------------------------------------------
// extend_a_kmer
//----------------------------------------------------------------------
__ONMIC__ void extend_a_kmer(const char* kmer, const char* sequence, const uint32_t& index_kmer, const uint32_t& index_last_mod,
        C_candidate_path& current_path, vector_candidate_path& candidate_path_vector, const uint32_t& org_boundary_left,
        const uint32_t& org_boundary_right, const char* quality_score, const unsigned char* bit_vector,
        const uint32_t* hash_seed, bool& run_exploration, init_args& arg, func_mem& funcm) {



    if (run_exploration == true) {
        // generate a new k-mer
        //      std::string kmer_new(kmer.substr(1, kmer_length - 1));
        //      kmer_new.push_back(sequence[index_kmer + kmer_length]);

        //	   char* kmer_new = (char*)malloc(arg.kmer_length + 2);
        //	   memcpy( kmer_new, kmer + 1, arg.kmer_length - 1);
        char* kmer_new = funcm.kmer + funcm.kmer_ind * funcm.kmer_length;
        copy_kmer(funcm, kmer + 1, arg.kmer_length - 1);
        kmer_new[arg.kmer_length - 1] =  sequence[index_kmer + arg.kmer_length];
        path_query_text(kmer_new, bit_vector,hash_seed,  arg.kmer_length, arg.num_hash_func, arg.bit_vector_width, 1);
        // this was the real boundary between solid k-mers and weak k-mers
        // check all possible cases
        if ((index_kmer == (org_boundary_left - 1)) || (index_kmer == (org_boundary_right - 1)) || (((unsigned short int)quality_score[index_kmer + arg.kmer_length] - arg.quality_score_offset) <= arg.extremely_low_quality_score)) {
            // each alternative neocletide
            for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
                // make a change
                kmer_new[arg.kmer_length - 1] = NEOCLEOTIDE[it_alter];

                // kmer_new is solid
//                if (query_text(kmer_new, bit_vector, &hash_seed[0],arg.kmer_length,arg.num_hash_func,arg.bit_vector_width) == true) {
                if(kmer_new[funcm.kmer_length + it_alter]){
                    // generate a new path
                    //               C_candidate_path temporary_path;
                    //               memcpy(&temporary_path,&current_path, sizeof(C_candidate_path));
                    //            	C_candidate_path *temporary_path ;
                    int increment = 0;
                    if (sequence[index_kmer + arg.kmer_length] != NEOCLEOTIDE[it_alter]) {
                        pair_uint_char pair_tmp;
                        pair_tmp.first  = index_kmer + arg.kmer_length;
                        pair_tmp.second = NEOCLEOTIDE[it_alter];

                        push_back(current_path, pair_tmp);
                        increment ++;
                    }

                    // if this k-mer is the last k-mer that can be modified
                    // running extend_a_kmer_right is not needed any more
                    if ((index_kmer + 1) == index_last_mod) {
                        push_back(candidate_path_vector, current_path);

                        // too many candidate paths
                        if (candidate_path_vector.size > MAX_CANDIDATE_PATHS) {
                            run_exploration = false;
                        }
                    }
                    else {
                        // trace  this kmer recursively and update candidate_path_vector
                        extend_a_kmer(
                                kmer_new,
                                sequence,
                                index_kmer + 1,
                                index_last_mod,
                                current_path,
                                candidate_path_vector,
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
                    current_path.size -= increment;
                }
            }
        }
        else {

            // kmer_new is a solid k-mer
            if (query_text(kmer_new, bit_vector, &hash_seed[0],arg.kmer_length,arg.num_hash_func,arg.bit_vector_width) == true) {

                // if this k-mer is the last k-mer that can be modified
                // running extend_a_kmer_right is not needed any more
                if ((index_kmer + 1) == index_last_mod) {
                    push_back(candidate_path_vector,current_path);

                    // too many candidate paths
                    if (candidate_path_vector.size > MAX_CANDIDATE_PATHS) {
                        run_exploration = false;
                    }
                }
                else {
                    extend_a_kmer(
                            kmer_new,
                            sequence,
                            index_kmer + 1,
                            index_last_mod,
                            current_path,
                            candidate_path_vector,
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
            else {
                path_query_text(kmer_new, bit_vector,hash_seed,  arg.kmer_length, arg.num_hash_func, arg.bit_vector_width, 1);
                // each alternative neocletide
                for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
                    // not equal to the original character
                    if (sequence[index_kmer + arg.kmer_length] != NEOCLEOTIDE[it_alter]) {
                        // make a change
                        kmer_new[arg.kmer_length - 1] = NEOCLEOTIDE[it_alter];

                        // kmer_new is solid
//                        if (query_text(kmer_new, bit_vector, &hash_seed[0],arg.kmer_length,arg.num_hash_func,arg.bit_vector_width) == true) {
                        if(kmer_new[funcm.kmer_length + it_alter]){
                            // generate a new path
                            //                     C_candidate_path temporary_path;
                            //                     memcpy(&temporary_path,&current_path, sizeof(C_candidate_path));
                            int increment = 0;
                            pair_uint_char pair_tmp;
                            pair_tmp.first  = index_kmer + arg.kmer_length;
                            pair_tmp.second = NEOCLEOTIDE[it_alter];

                            push_back(current_path, pair_tmp);

                            increment += 1;
                            // if this k-mer is the last k-mer that can be modified
                            // running extend_a_kmer_right is not needed any more
                            if ((index_kmer + 1) == index_last_mod) {
                                push_back(candidate_path_vector, current_path);

                                // too many candidate paths
                                if (candidate_path_vector.size > MAX_CANDIDATE_PATHS) {
                                    run_exploration = false;
                                }
                            }
                            else {
                                // trace  this kmer recursively and update candidate_path_vector
                                extend_a_kmer(
                                        kmer_new,
                                        sequence,
                                        index_kmer + 1,
                                        index_last_mod,
                                        current_path,
                                        candidate_path_vector,
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
                            current_path.size -= increment;
                        }
                    }
                }
            }
        }
        //      free(kmer_new);
        funcm.kmer_ind --;
    }
}


//----------------------------------------------------------------------
// base_intersection
//----------------------------------------------------------------------
__ONMIC__ inline void base_intersection(pair_uint_char* in1_begin, const pair_uint_char* in1_end, pair_uint_char* in2_begin, const pair_uint_char* in2_end, vector_pair_uint_char& out) {
    while ((in1_begin != in1_end) && (in2_begin != in2_end)) {
        if ((*in1_begin).first < (*in2_begin).first) {
            in1_begin++;
        }
        else if ((*in2_begin).first < (*in1_begin).first) {
            in2_begin++;
        }
        else {
            if ((*in1_begin).second == (*in2_begin).second) {
                push_back(out, *in1_begin);
            }

            in1_begin++;
            in2_begin++;
        }
    }
}

//----------------------------------------------------------------------
// base_union
//----------------------------------------------------------------------
__ONMIC__ inline void base_union(pair_uint_char* in1_begin, const pair_uint_char* in1_end, pair_uint_char* in2_begin, const pair_uint_char* in2_end, vector_pair_uint_char& out) {
    while (true) {
        if (in1_begin == in1_end) {
            insert(out.uint_char[out.size], in2_begin, in2_end);
            return;
        }
        else if (in2_begin == in2_end) {
            //         out.insert(out.end(), in1_begin, in1_end);
            insert(out.uint_char[out.size], in2_begin, in2_end);
            return;
        }
        else {
            if (((*in1_begin).first) < ((*in2_begin).first)) {
                push_back(out,*in1_begin);
                in1_begin++;
            }
            else if (((*in2_begin).first) < ((*in1_begin).first)) {
                push_back(out,*in2_begin);
                in2_begin++;
            }
            else {
                push_back(out, *in1_begin);
                in1_begin++;
                in2_begin++;
            }
        }
    }
}



//----------------------------------------------------------------------
// base_difference
//----------------------------------------------------------------------
__ONMIC__ inline void base_difference(pair_uint_char* in1_begin, const pair_uint_char* in1_end, pair_uint_char* in2_begin, const pair_uint_char* in2_end, vector_pair_uint_char& out) {
    while ((in1_begin != in1_end) && (in2_begin != in2_end)) {
        if ((*in1_begin).first < (*in2_begin).first) {
            push_back(out, *in1_begin);
            in1_begin++;
        }
        else if ((*in2_begin).first < (*in1_begin).first) {
            in2_begin++;
        }
        else {
            in1_begin++;
            in2_begin++;
        }
    }

    insert(out.uint_char[out.size], in1_begin, in1_end);
}



//----------------------------------------------------------------------
// correct_errors_between_solid_regions
//----------------------------------------------------------------------
__ONMIC__ void correct_errors_between_solid_regions(const char* org_sequence, char* sequence, const char* quality_score,
        const uint32_t left_first, const uint32_t index_start, const uint32_t index_end, const uint32_t right_second,
        const uint32_t org_boundary_left, const uint32_t org_boundary_right, char* sequence_modification,
        uint32_t& trim_5_end, uint32_t& trim_3_end, const uint32_t& num_solid_islands, const uint32_t& read_length,
        const uint32_t& max_trimmed_bases, const uint32_t& min_check_length,
        uint32_t& num_corrected_errors_local, const uint8_t * bit_vector,
        const uint32_t* hash_seed, init_args& arg, func_mem& funcm) {
    //--------------------------------------------------
    // from i-th region to (i + 1)-th region
    //--------------------------------------------------
    // i-th solid region | non-solid | (i+1)-th solid region
    // --------------------------------------- read
    //              |-----|                    (index_start)-th k-mer
    //                              |-----|    (index_end)-th k-mer: last non-solid k-mer
    //                         |-----|         (index_last_mod)-th k-mer = (index_end - kmer_length + 1)-th k-mer: last k-mer that can be modified
    //                               |----|    This regions should not be modified
    //--------------------------------------------------
    // list of candidate paths

    vector_candidate_path candidate_path_vector_tmp;
    init_vector_candidate_path(candidate_path_vector_tmp);


    // index of the k-mer that can be modified
    // k-mers that are overlapped with a solid regioin cannot be modified
    uint32_t index_last_mod=(index_end - arg.kmer_length + 1);

    char* kmer_initial = funcm.kmer + funcm.kmer_ind * funcm.kmer_length;
    // make an initial k-mer
    //   std::string kmer_initial(sequence.substr(index_start, kmer_length));
    //  memcpy(kmer_initial, sequence + index_start, arg.kmer_length);
    copy_kmer(funcm, sequence + index_start, arg.kmer_length);
    path_query_text(kmer_initial, bit_vector,hash_seed,  arg.kmer_length, arg.num_hash_func, arg.bit_vector_width, 1);

    // each alternative neocletide
    bool run_exploration = (true);

    for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
        // make a change
        kmer_initial[arg.kmer_length - 1] = NEOCLEOTIDE[it_alter];

        // kmer_initial is solid
//        if ( query_text(kmer_initial, bit_vector, &hash_seed[0],arg.kmer_length,arg.num_hash_func,arg.bit_vector_width) == true) {
        if(kmer_initial[funcm.kmer_length + it_alter]){

            // generate a new path
            C_candidate_path candidate_path;
            candidate_path.size = 0;
            candidate_path.sum_qs = 0;
            if (sequence[index_start + arg.kmer_length - 1] != NEOCLEOTIDE[it_alter]) {
                pair_uint_char pair_tmp;
                pair_tmp.first  = index_start + arg.kmer_length - 1;
                pair_tmp.second = NEOCLEOTIDE[it_alter];

                push_back(candidate_path, pair_tmp);
            }

            // if this k-mer is the last k-mer that can be modified
            // running extend_a_kmer_right is not needed any more
            if (index_start == index_last_mod) {
                //		 candidate_path_vector_tmp.push_back(candidate_path);
                push_back(candidate_path_vector_tmp, candidate_path);
            }
            else {
                // trace  this kmer recursively and update candidate_path_vector_tmp
                extend_a_kmer(
                        kmer_initial,
                        sequence,
                        index_start,
                        index_last_mod,
                        candidate_path,
                        candidate_path_vector_tmp,
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
    // complete exploration was not done because there are too many candidata paths
    // remove all the paths in candidate_path_vector_tmp
    if (run_exploration == false) {
        clear(candidate_path_vector_tmp);
    }

    vector_candidate_path candidate_path_vector;
    init_vector_candidate_path(candidate_path_vector);

    // check the solidness of k-mers between index_last_mod and index_end
    bool all_solid_wo_modification = (false);

    // each candidate path
    for (C_candidate_path* it_path = candidate_path_vector_tmp.c_path; it_path != candidate_path_vector_tmp.c_path + candidate_path_vector_tmp.size; it_path++) {
        if ((*it_path).size == 0) {
            all_solid_wo_modification = true;
            break;
        }
        else {
            // checking is needed
            uint32_t index_last_modified_base = ((*it_path).modified_bases[(*it_path).size - 1].first);

            if (index_last_modified_base > index_last_mod) {
                // generate a temporary sequence
                //            std::string sequence_tmp(sequence);
                char* sequence_tmp;// = (char*)malloc(READ_MAX_LEN);
                init_sequence(funcm, &sequence_tmp);
                //			   strcpy(sequence_tmp, sequence);
                memcpy(sequence_tmp,sequence, READ_MAX_LEN);
                for (unsigned int it_base = 0; it_base < (*it_path).size; it_base++) {
                    sequence_tmp[(*it_path).modified_bases[it_base].first] = (*it_path).modified_bases[it_base].second;
                }

                // check k-mers
                uint32_t num_success = (0);
                uint32_t num_fails = (0);
                for (uint32_t it_check = index_last_mod; it_check <= index_last_modified_base; it_check++) {
                    //               std::string kmer_current(sequence_tmp.substr(it_check, kmer_length));
                    //				   char* kmer_current = (char*) malloc(arg.kmer_length + 2);
                    //				   memcpy(kmer_current, sequence_tmp + it_check, arg.kmer_length);
                    char* kmer_current = funcm.kmer + funcm.kmer_ind * funcm.kmer_length;
                    copy_kmer(funcm, sequence_tmp + it_check, arg.kmer_length);

                    if (query_text(kmer_current, bit_vector, &hash_seed[0],arg.kmer_length,arg.num_hash_func,arg.bit_vector_width) == true) {
                        num_success++;
                    }
                    else {
                        num_fails++;

                        if (num_fails > NUM_ALLOWABLE_FAILS) {
                            break;
                        }
                    }
                    funcm.kmer_ind --;
                }

                if (num_success >= (index_last_modified_base - index_last_mod + 1 - NUM_ALLOWABLE_FAILS)) {
                    push_back(candidate_path_vector,*it_path);
                }
                funcm.sequence_ind --;
            }
            // checking is not needed
            else {
                push_back(candidate_path_vector, *it_path);
            }
        }
    }

    // remain only really modified paths
    //   candidate_path_vector_tmp.clear();
    clear(candidate_path_vector_tmp);
    // each path
    for (uint32_t it_candidate = 0; it_candidate < candidate_path_vector.size; it_candidate++) {
        // each modification
        bool really_modified = (false);
        for (unsigned int it_mod_base = 0; it_mod_base < candidate_path_vector.c_path[it_candidate].size; it_mod_base++) {
            if (org_sequence[candidate_path_vector.c_path[it_candidate].modified_bases[it_mod_base].first] != candidate_path_vector.c_path[it_candidate].modified_bases[it_mod_base].second) {
                really_modified = true;
            }
        }

        if (really_modified) {
            push_back(candidate_path_vector_tmp, candidate_path_vector.c_path[it_candidate]);
        }
    }

    //   candidate_path_vector = candidate_path_vector_tmp;
    vector_candidate_path_assignment(&candidate_path_vector, &candidate_path_vector_tmp);

    // all k-mers are solid without any modification
    if (all_solid_wo_modification == true) {
        // do nothing
    }
    // compare quality scores of candidate paths
    // if the number of paths in candidate_path_vector is larger than 1
    else if (candidate_path_vector.size > 1) {
        // each path
        C_candidate_path*  it_path;
        C_candidate_path* it_path_1st;
        C_candidate_path* it_path_2nd;

        uint32_t qs_1st = (INIT_MIN_QS);
        uint32_t qs_2nd = (INIT_MIN_QS);

        for (it_path = candidate_path_vector.c_path; it_path != candidate_path_vector.c_path + candidate_path_vector.size; it_path++) {
            // each modification
            for (pair_uint_char* it_mod = (*it_path).modified_bases; it_mod != (*it_path).modified_bases + (*it_path).size; it_mod++) {
                // add quality scores of modified bases
                (*it_path).sum_qs += ((unsigned short int)quality_score[(*it_mod).first] - arg.quality_score_offset);
            }

            // compare quality scores of each path
            if ((*it_path).sum_qs <= qs_1st) {
                qs_2nd = qs_1st;
                qs_1st = (*it_path).sum_qs;

                it_path_2nd = it_path_1st;
                it_path_1st = it_path;
            }
            else if ((*it_path).sum_qs <= qs_2nd) {
                qs_2nd = (*it_path).sum_qs;

                it_path_2nd = it_path;
            }
        }

        // check whether too many bases are modified in the first path, which implies indel may exist in the original read
        // if too many modifications exist, the read is just trimmed
        bool too_many_corrections = (false);

        // this is done only when the number of solid islands is two
        if (num_solid_islands == 2) {
            // the 1st island is small && the second island is big
            if (((index_start - left_first) > MAX_MODIFICATION) && ((right_second - index_end + 2) <= MAX_MODIFICATION)) {
                // at least over MAX_MODIFICATION times of modifications from the first modified position
                if ((((*it_path_1st).modified_bases[(*it_path_1st).size - 1].first -arg.kmer_length - index_start + 2) >= min_check_length) &&
                        ((*it_path_1st).size > MAX_MODIFICATION)) {
                    // each modified base in a right range
                    //std::uint32_t partial_qs((*it_path_1st).sum_qs);

                    for (uint32_t it_mod_base = ((*it_path_1st).size - 1); it_mod_base >= MAX_MODIFICATION; it_mod_base--) {
                        // at least MAX_MODIFICATION times of modifications within min_check_length bases
                        if (((*it_path_1st).modified_bases[it_mod_base].first - (*it_path_1st).modified_bases[it_mod_base - MAX_MODIFICATION].first) < min_check_length) {
                            //if (it_mod_base < ((*it_path_1st).modified_bases.size() - 1)) {
                            //   partial_qs -= ((unsigned short int)quality_score[(*it_path_1st).modified_bases[it_mod_base + 1].first] - quality_score_offset);
                            //}

                            // average quality score of the modified bases is too high
                            //if (1.0 * partial_qs / (it_mod_base + 1) > quality_score_cutoff) {
                            // trim_5_end or trim_3_end
                            if ((read_length - (*it_path_1st).modified_bases[it_mod_base].first) <= max_trimmed_bases) {
                                trim_3_end = read_length - (*it_path_1st).modified_bases[it_mod_base].first;
                            }
                            else if (((*it_path_1st).modified_bases[it_mod_base].first + 1) <= max_trimmed_bases) {
                                trim_5_end = (*it_path_1st).modified_bases[it_mod_base].first + 1;
                            }

                            too_many_corrections = true;
                            break;
                            //}
                        }
                    }
                }
            }
            // the 1st island is big && the second island is small
            else if (((index_start - left_first) <= MAX_MODIFICATION) && ((right_second - index_end + 2) > MAX_MODIFICATION)) {
                // at least over MAX_MODIFICATION times of modifications from the first modified position
                if (((index_end - (*it_path_1st).modified_bases[0].first + 1) >= min_check_length) &&
                        ((*it_path_1st).size > MAX_MODIFICATION)) {
                    // each modified base in a right range
                    //std::uint32_t partial_qs((*it_path_1st).sum_qs);

                    for (unsigned int it_mod_base = 0; it_mod_base < ((*it_path_1st).size - MAX_MODIFICATION); it_mod_base++) {
                        // at least MAX_MODIFICATION times of modifications within min_check_length bases
                        if (((*it_path_1st).modified_bases[it_mod_base + MAX_MODIFICATION].first - (*it_path_1st).modified_bases[it_mod_base].first) < min_check_length) {
                            //if (it_mod_base > 0) {
                            //   partial_qs -= ((unsigned short int)quality_score[(*it_path_1st).modified_bases[it_mod_base - 1].first] - quality_score_offset);
                            //}

                            // average quality score of the modified bases is too high
                            //if (1.0 * partial_qs / ((*it_path_1st).modified_bases.size() - it_mod_base) > quality_score_cutoff) {
                            // trim_5_end or trim_3_end
                            if ((read_length - (*it_path_1st).modified_bases[it_mod_base].first) <= max_trimmed_bases) {
                                trim_3_end = read_length - (*it_path_1st).modified_bases[it_mod_base].first;
                            }
                            else if (((*it_path_1st).modified_bases[it_mod_base].first + 1) <= max_trimmed_bases) {
                                trim_5_end = (*it_path_1st).modified_bases[it_mod_base].first + 1;
                            }

                            too_many_corrections = true;
                            break;
                            //}
                        }
                    }
                }
            }
        }

        // use the 1st path if not too many corrections are made AND if the 1st path has a sufficiently low score
        // if an indel exists using quality scores is not a good way to choose the best path
        if ((too_many_corrections == false) && (qs_1st + MIN_QS_DIFF <= qs_2nd)) {
            // each modification
            for (pair_uint_char* it_base = (*it_path_1st).modified_bases; it_base != (*it_path_1st).modified_bases + (*it_path_1st).size; it_base++) {
                // update sequence_modification
                sequence_modification[(*it_base).first] = (*it_base).second;
                sequence[(*it_base).first] = (*it_base).second;
                num_corrected_errors_local++;
            }
        }
        // hard to choose one path
        else {
            // A AND B
            vector_pair_uint_char base_vector_intersection;
            base_vector_intersection.size = 0;
            base_vector_intersection.sum_qs = 0;
            // A OR B
            vector_pair_uint_char base_vector_union;
            base_vector_union.size = 0;
            base_vector_union.sum_qs = 0;
            // temporary vectors
            vector_pair_uint_char base_vector_intersection_prev;


            memcpy(&base_vector_intersection_prev,&(candidate_path_vector.c_path[0]),sizeof(vector_pair_uint_char));
            vector_pair_uint_char base_vector_union_prev;
            memcpy(&base_vector_union_prev,&(candidate_path_vector.c_path[0]),sizeof(vector_pair_uint_char));
            //		   (candidate_path_vector[0].modified_bases);

            // each candidate path
            for (uint32_t it_p = 1; it_p < candidate_path_vector.size; it_p++) {
                //			   base_vector_intersection.clear();
                clear(base_vector_intersection);

                //			   base_vector_union.clear();
                clear(base_vector_union);

                base_intersection(base_vector_intersection_prev.uint_char, base_vector_intersection_prev.uint_char + base_vector_intersection_prev.size,
                        candidate_path_vector.c_path[it_p].modified_bases, candidate_path_vector.c_path[it_p].modified_bases + candidate_path_vector.c_path[it_p].size, base_vector_intersection);
                //			   base_union       (base_vector_union_prev.begin(),        base_vector_union_prev.end(),        candidate_path_vector[it_p].modified_bases.begin(), candidate_path_vector[it_p].modified_bases.end(), base_vector_union);
                base_union (base_vector_intersection_prev.uint_char, base_vector_intersection_prev.uint_char + base_vector_intersection_prev.size,
                        candidate_path_vector.c_path[it_p].modified_bases, candidate_path_vector.c_path[it_p].modified_bases + candidate_path_vector.c_path[it_p].size, base_vector_intersection);
                //			   base_vector_intersection_prev = base_vector_intersection;
                //			   base_vector_union_prev        = base_vector_union;
                memcpy(&base_vector_intersection, &base_vector_intersection, sizeof(vector_pair_uint_char));
                memcpy(&base_vector_union_prev, &base_vector_union, sizeof(vector_pair_uint_char));
            }

            // A - B
            vector_pair_uint_char base_vector_difference;
            base_vector_difference.size = 0;
            base_vector_difference.sum_qs = 0;
            base_difference(base_vector_union.uint_char, base_vector_union.uint_char + base_vector_union.size , base_vector_intersection.uint_char, base_vector_intersection.uint_char + base_vector_intersection.size, base_vector_difference);

            // find trimmed region
            // correcting the 5'-end and 3'-end was not done
            // therefore the total number of trimmed bases is 0 yet
            if (base_vector_difference.size > 0) {
                uint32_t vector_index_leftmost = (0);
                uint32_t vector_index_rightmost = (base_vector_difference.size - 1);

                bool keep_going = (true);

                while (keep_going) {
                    // # of trimmed bases at the 5' end: base_vector_difference[vector_index_leftmost].first + 1
                    // # of trimmed bases at the 3' end: read_length - base_vector_difference[vector_index_rightmost].first
                    // the 5'-end is smaller
                    if ((base_vector_difference.uint_char[vector_index_leftmost].first + 1) < (read_length - base_vector_difference.uint_char[vector_index_rightmost].first)) {
                        // check the total number of trimmed bases
                        if ((base_vector_difference.uint_char[vector_index_leftmost].first + 1 + trim_3_end) <= max_trimmed_bases) {
                            if ((base_vector_difference.uint_char[vector_index_leftmost].first + 1) > trim_5_end) {
                                trim_5_end = base_vector_difference.uint_char[vector_index_leftmost].first + 1;
                            }

                            // two points are met
                            if (vector_index_leftmost == vector_index_rightmost) {
                                keep_going = false;
                            }
                            else {
                                vector_index_leftmost++;
                            }
                        }
                        // no need for more check
                        else {
                            keep_going = false;
                        }
                    }
                    // the 3'-end is smaller
                    else {
                        // check the total number of trimmed bases
                        if ((read_length - base_vector_difference.uint_char[vector_index_rightmost].first) <= max_trimmed_bases) {
                            if ((read_length - base_vector_difference.uint_char[vector_index_rightmost].first) > trim_3_end) {
                                trim_3_end = read_length - base_vector_difference.uint_char[vector_index_rightmost].first;
                            }

                            // two points are met
                            if (vector_index_leftmost == vector_index_rightmost) {
                                keep_going = false;
                            }
                            else {
                                vector_index_rightmost--;
                            }
                        }
                        // no need for more check
                        else {
                            keep_going = false;
                        }
                    }
                }
            }

            // find consensus modifications
            for (uint32_t it_inter = 0; it_inter < base_vector_intersection.size; it_inter++) {
                // check whether the base is not in the trimmed regions
                if ((base_vector_intersection.uint_char[it_inter].first <  (read_length - trim_3_end)) &&
                        (base_vector_intersection.uint_char[it_inter].first >= trim_5_end)) {
                    // filter out the bases that are equal to the original ones
                    if (sequence[base_vector_intersection.uint_char[it_inter].first] != base_vector_intersection.uint_char[it_inter].second) {
                        sequence_modification[base_vector_intersection.uint_char[it_inter].first] = base_vector_intersection.uint_char[it_inter].second;
                        sequence[base_vector_intersection.uint_char[it_inter].first] = base_vector_intersection.uint_char[it_inter].second;

                        num_corrected_errors_local++;
                    }
                }
            }
        }
    }

    // only one path
    else if (candidate_path_vector.size == 1) {
        // check whether too many bases are modified in the first path, which implies indel may exist in the original read
        // if too many modifications exist, the read is just trimmed
        bool too_many_corrections(false);

        // this is done only when the number of solid islands is two
        if (num_solid_islands == 2) {
            // the 1st island is small && the second island is big
            if (((index_start - left_first) > MAX_MODIFICATION) && ((right_second - index_end + 2) <= MAX_MODIFICATION)) {
                // at least over MAX_MODIFICATION times of modifications from the first modified position
                if (((candidate_path_vector.c_path[0].modified_bases[candidate_path_vector.c_path[0].size - 1].first -arg.kmer_length - index_start + 2) >= min_check_length) &&
                        (candidate_path_vector.c_path[0].size > MAX_MODIFICATION)) {
                    // calculate sum_qs
                    //for (unsigned int it_mod_base = 0; it_mod_base < candidate_path_vector[0].modified_bases.size(); it_mod_base++) {
                    //   candidate_path_vector[0].sum_qs += ((unsigned short int)quality_score[candidate_path_vector[0].modified_bases[it_mod_base].first] - quality_score_offset);
                    //}

                    // each modified base in a right range
                    //std::uint32_t partial_qs(candidate_path_vector[0].sum_qs);

                    for (uint32_t it_mod_base = (candidate_path_vector.c_path[0].size - 1); it_mod_base >= MAX_MODIFICATION; it_mod_base--) {
                        // at least MAX_MODIFICATION times of modifications within min_check_length bases
                        if ((candidate_path_vector.c_path[0].modified_bases[it_mod_base].first - candidate_path_vector.c_path[0].modified_bases[it_mod_base - MAX_MODIFICATION].first) < min_check_length) {
                            //if (it_mod_base < (candidate_path_vector[0].modified_bases.size() - 1)) {
                            //   partial_qs -= ((unsigned short int)quality_score[candidate_path_vector[0].modified_bases[it_mod_base + 1].first] - quality_score_offset);
                            //}

                            // average quality score of the modified bases is too high
                            //if (1.0 * partial_qs / (it_mod_base + 1) > quality_score_cutoff) {
                            // trim_5_end or trim_3_end
                            if ((read_length - candidate_path_vector.c_path[0].modified_bases[it_mod_base].first) <= max_trimmed_bases) {
                                trim_3_end = read_length - candidate_path_vector.c_path[0].modified_bases[it_mod_base].first;

                                // update sequence_modification for the non-trimmed corrections
                                if (it_mod_base > 0) {
                                    for (unsigned int it_base = 0; it_base < (it_mod_base - 1); it_base++) {
                                        sequence_modification[(candidate_path_vector.c_path[0]).modified_bases[it_base].first] = (candidate_path_vector.c_path[0]).modified_bases[it_base].second;
                                        sequence[(candidate_path_vector.c_path[0]).modified_bases[it_base].first] = (candidate_path_vector.c_path[0]).modified_bases[it_base].second;
                                        num_corrected_errors_local++;
                                    }
                                }
                            }
                            else if ((candidate_path_vector.c_path[0].modified_bases[it_mod_base].first + 1) <= max_trimmed_bases) {
                                trim_5_end = candidate_path_vector.c_path[0].modified_bases[it_mod_base].first + 1;

                                // update sequence_modification for the non-trimmed corrections
                                if (it_mod_base < (candidate_path_vector.c_path[0].size - 1)) {
                                    for (unsigned int it_base = (it_mod_base + 1); it_base < candidate_path_vector.c_path[0].size; it_base++) {
                                        sequence_modification[(candidate_path_vector.c_path[0]).modified_bases[it_base].first] = (candidate_path_vector.c_path[0]).modified_bases[it_base].second;
                                        sequence[(candidate_path_vector.c_path[0]).modified_bases[it_base].first] = (candidate_path_vector.c_path[0]).modified_bases[it_base].second;
                                        num_corrected_errors_local++;
                                    }
                                }
                            }

                            too_many_corrections = true;
                            break;
                            //}
                        }
                    }
                }
            }
            // the 1st island is big && the second island is small
            else if (((index_start - left_first) <= MAX_MODIFICATION) && ((right_second - index_end + 2) > MAX_MODIFICATION)) {
                // at least over MAX_MODIFICATION times of modifications from the first modified position
                if (((index_end - candidate_path_vector.c_path[0].modified_bases[0].first + 1) >= min_check_length) &&
                        (candidate_path_vector.c_path[0].size > MAX_MODIFICATION)) {
                    // each modified base in a right range
                    //std::uint32_t partial_qs(candidate_path_vector[0].sum_qs);

                    for (unsigned int it_mod_base = 0; it_mod_base < (candidate_path_vector.c_path[0].size - MAX_MODIFICATION); it_mod_base++) {
                        // at least MAX_MODIFICATION times of modifications within min_check_length bases
                        if ((candidate_path_vector.c_path[0].modified_bases[it_mod_base + MAX_MODIFICATION].first - candidate_path_vector.c_path[0].modified_bases[it_mod_base].first) < min_check_length) {
                            //if (it_mod_base > 0) {
                            //   partial_qs -= ((unsigned short int)quality_score[candidate_path_vector[0].modified_bases[it_mod_base - 1].first] - quality_score_offset);
                            //}

                            // average quality score of the modified bases is too high
                            //if (1.0 * partial_qs / (candidate_path_vector[0].modified_bases.size() - it_mod_base) > quality_score_cutoff) {
                            // trim_5_end or trim_3_end
                            if ((read_length - candidate_path_vector.c_path[0].modified_bases[it_mod_base].first) <= max_trimmed_bases) {
                                trim_3_end = read_length - candidate_path_vector.c_path[0].modified_bases[it_mod_base].first;
                            }
                            else if ((candidate_path_vector.c_path[0].modified_bases[it_mod_base].first + 1) <= max_trimmed_bases) {
                                trim_5_end = candidate_path_vector.c_path[0].modified_bases[it_mod_base].first + 1;
                            }

                            too_many_corrections = true;
                            break;
                            //}
                        }
                    }
                }
            }
        }

        // not too many modified bases
        if (too_many_corrections == false) {
            // each modification
            for (unsigned int it_base = 0; it_base < candidate_path_vector.c_path[0].size; it_base++) {
                // update sequence_modification
                sequence_modification[candidate_path_vector.c_path[0].modified_bases[it_base].first] = candidate_path_vector.c_path[0].modified_bases[it_base].second;
                sequence[candidate_path_vector.c_path[0].modified_bases[it_base].first] = candidate_path_vector.c_path[0].modified_bases[it_base].second;

                num_corrected_errors_local++;
            }
        }


    }

    // no path
    // if indels exist between two solid k-mer islands, checking the solidness of k-mers between index_last_mod and index_end will always fails
    // this kind of errors cannot be corrected without modifying existing solid k-mer islands
    // one of the both the sides should be trimmed
    else if (candidate_path_vector.size == 0) {
        // trim all the bases to the right of the 1st solid k-mer island
        if ((read_length - (index_start + arg.kmer_length - 1)) <= max_trimmed_bases) {
            trim_3_end = read_length - (index_start + arg.kmer_length - 1);
        }
        // trim all the bases to the left of the 2nd solid k-mer island
        else if ((index_end + 1) <= max_trimmed_bases) {
            trim_5_end = index_end + 1;
        }
        // trim all the bases to the right of the 1st solid k-mer island
        else if ((read_length - (org_boundary_left + arg.kmer_length - 1)) <= max_trimmed_bases) {
            trim_3_end = read_length - (org_boundary_left + arg.kmer_length - 1);
        }
        // trim all the bases to the left of the 2nd solid k-mer island
        else if ((org_boundary_right + 1) <= max_trimmed_bases) {
            trim_5_end = org_boundary_right + 1;
        }
    }

    //   free(kmer_initial);
    funcm.kmer_ind --;
    free(candidate_path_vector_tmp.c_path);
    free(candidate_path_vector.c_path);

}

//----------------------------------------------------------------------
// sort_indexes
//----------------------------------------------------------------------
//inline uint32_t*  sort_indexes(const char* in_string, const std::uint32_t& str_length) {
//	std::vector<unsigned int> index_vector(str_length);
//	for (std::uint32_t it = 0; it != index_vector.size(); it++) {
//		index_vector[it] = it;
//	}
//
//	std::sort(index_vector.begin(), index_vector.end(), [&in_string](std::uint32_t i1, std::uint32_t i2) {return (unsigned short int)in_string[i1] < (unsigned short int)in_string[i2];});
//
//	return index_vector;
//}
__ONMIC__ inline uint32_t* sort_indexes(const char* in_string, const uint32_t& str_length) {
    uint32_t*  index_vector = (uint32_t*)malloc(str_length * sizeof(uint32_t));
    for (uint32_t it = 0; it != str_length; it++) {
        index_vector[it] = it;
    }

    //	std::sort(index_vector.begin(), index_vector.end(), [&in_string](std::uint32_t i1, std::uint32_t i2) {return (unsigned short int)in_string[i1] < (unsigned short int)in_string[i2];});

    for(uint32_t it = 0; it != str_length; it++){
        for(uint32_t it2 = it; it2 > 0; it2--){
            if(in_string[index_vector[it2]] < in_string[index_vector[it2-1]]){
                uint32_t t = index_vector[it2];
                index_vector[it2] = index_vector[it2-1];
                index_vector[it2-1] = t;

            }
        }
    }
    return index_vector;
}

//----------------------------------------------------------------------
// extend_a_kmer_5_prime_end
//----------------------------------------------------------------------
__ONMIC__ void extend_a_kmer_5_prime_end(const char* kmer, const char* sequence,
        const uint32_t& index_kmer, C_candidate_path& current_path, vector_candidate_path& candidate_path_vector,
        const uint32_t& org_boundary, const char* quality_score, const unsigned char* bit_vector,
        const uint32_t* hash_seed, bool& run_exploration, init_args& arg, func_mem& funcm) {
    if (run_exploration == true) {
        // generate a new k-mer
        //      std::string kmer_new(kmer.substr(0, kmer_length - 1));
        //      kmer_new = sequence[index_kmer - 1] + kmer_new;
        //		char* kmer_new = (char*)malloc(arg.kmer_length);
        //		memcpy(kmer_new, kmer + 1, arg.kmer_length-1);
        char* kmer_new = funcm.kmer + funcm.kmer_ind * funcm.kmer_length;
        copy_kmer(funcm, kmer, arg.kmer_length-1, 1);
        kmer_new[0] = sequence[index_kmer - 1];


        // this was the real boundary between solid k-mers and weak k-mers
        // check all possible cases
        if ((index_kmer == (org_boundary + 1)) || (((unsigned short int)quality_score[index_kmer - 1] - arg.quality_score_offset) <= arg.extremely_low_quality_score)) {
            path_query_text(kmer_new, bit_vector,hash_seed,  arg.kmer_length, arg.num_hash_func, arg.bit_vector_width, 0);
        	// each alternative neocletide
            for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
                // make a change
                kmer_new[0] = NEOCLEOTIDE[it_alter];

                // kmer_new is solid
//                if (query_text(kmer_new, bit_vector, &hash_seed[0],arg.kmer_length, arg.num_hash_func, arg.bit_vector_width)  == true) {
                if(kmer_new[funcm.kmer_length + it_alter]){
                    // generate a new path
                    //					C_candidate_path temporary_path;
                    //					memcpy(&temporary_path, &current_path, sizeof(C_candidate_path));

                    int increment = 0;
                    // not equal to the original character
                    if (sequence[index_kmer - 1] != NEOCLEOTIDE[it_alter]) {
                        pair_uint_char pair_tmp;
                        pair_tmp.first  = index_kmer - 1;
                        pair_tmp.second = NEOCLEOTIDE[it_alter];

                        push_back(current_path,pair_tmp);
                        increment ++;
                    }

                    // if this k-mer is the first k-mer in a read
                    // running extend_a_kmer_5_prime_end is not needed any more
                    if ((index_kmer - 1) == 0) {
                        push_back(candidate_path_vector, current_path);

                        // too many candidate paths
                        if (candidate_path_vector.size > MAX_CANDIDATE_PATHS) {
                            run_exploration = false;
                        }
                    }
                    else if ((index_kmer - 1) > 0) {
                        // trace  this kmer recursively and update candidate_path_vector
                        extend_a_kmer_5_prime_end(
                                kmer_new,
                                sequence,
                                index_kmer - 1,
                                current_path,
                                candidate_path_vector,
                                org_boundary,
                                quality_score,
                                bit_vector,
                                hash_seed,
                                run_exploration,
                                arg,
                                funcm
                                );
                    }
                    current_path.size -= increment;
                }
            }
        }
        else {
            // kmer_new is a solid k-mer
            if (query_text(kmer_new, bit_vector, &hash_seed[0],arg.kmer_length, arg.num_hash_func, arg.bit_vector_width) == true) {
                // if this k-mer is the first k-mer in a read
                // running extend_a_kmer_5_prime_end is not needed any more
                if ((index_kmer - 1) == 0) {
                    push_back(candidate_path_vector,current_path);

                    // too many candidate paths
                    if (candidate_path_vector.size > MAX_CANDIDATE_PATHS) {
                        run_exploration = false;
                    }
                }
                else if ((index_kmer - 1) > 0) {
                    extend_a_kmer_5_prime_end(
                            kmer_new,
                            sequence,
                            index_kmer - 1,
                            current_path,
                            candidate_path_vector,
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
            else {
                path_query_text(kmer_new, bit_vector,hash_seed,  arg.kmer_length, arg.num_hash_func, arg.bit_vector_width, 0);
                // each alternative neocletide
                for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
                    // not equal to the original character
                    if (sequence[index_kmer - 1] != NEOCLEOTIDE[it_alter]) {
                        // make a change
                        kmer_new[0] = NEOCLEOTIDE[it_alter];

                        // kmer_new is solid
//                        if (query_text(kmer_new, bit_vector, &hash_seed[0],arg.kmer_length, arg.num_hash_func, arg.bit_vector_width) == true) {
                        if(kmer_new[funcm.kmer_length + it_alter]){
                            // generate a new path
                            //							C_candidate_path temporary_path(current_path);

                            pair_uint_char pair_tmp;
                            pair_tmp.first  = index_kmer - 1;
                            pair_tmp.second = NEOCLEOTIDE[it_alter];

                            push_back(current_path, pair_tmp);

                            // if this k-mer is the first k-mer in a read
                            // running extend_a_kmer_5_prime_end is not needed any more
                            if ((index_kmer - 1) == 0) {
                                push_back(candidate_path_vector, current_path);

                                // too many candidate paths
                                if (candidate_path_vector.size > MAX_CANDIDATE_PATHS) {
                                    run_exploration = false;
                                }
                            }
                            else if ((index_kmer - 1) > 0) {
                                // trace  this kmer recursively and update candidate_path_vector
                                extend_a_kmer_5_prime_end(
                                        kmer_new,
                                        sequence,
                                        index_kmer - 1,
                                        current_path,
                                        candidate_path_vector,
                                        org_boundary,
                                        quality_score,
                                        bit_vector,
                                        hash_seed,
                                        run_exploration,
                                        arg,
                                        funcm
                                        );
                            }
                            current_path.size --;
                        }
                    }
                }
            }
        }
        //		free(kmer_new);
        funcm.kmer_ind --;
    }
}
//----------------------------------------------------------------------
// extend_out_left
//----------------------------------------------------------------------
__ONMIC__ inline void extend_out_left(const char* kmer, const uint32_t& num_extend, const uint32_t& extend_amount,
        bool& extension_success, const unsigned char* bit_vector, const uint32_t* hash_seed, init_args& arg, func_mem& funcm) {
    // generate a new k-mer
    //   std::string kmer_new(kmer.substr(0, kmer_length - 1));
    //   kmer_new = '0' + kmer_new;
    //	char* kmer_new = (char*) malloc(arg.kmer_length + 2);
    //	memcpy(kmer_new + 1, kmer , arg.kmer_length -1);
    char* kmer_new = funcm.kmer + funcm.kmer_ind * funcm.kmer_length;
    copy_kmer(funcm, kmer , arg.kmer_length -1, 1);
    kmer_new[0] = 'A';
    path_query_text(kmer_new, bit_vector,hash_seed,  arg.kmer_length, arg.num_hash_func, arg.bit_vector_width, 0);
    kmer_new[0] = '0';
    // each alternative neocletide
    for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
        // generate kmer_new
        kmer_new[0] = NEOCLEOTIDE[it_alter];

        // kmer_new is solid
//        if (query_text(kmer_new, bit_vector, &hash_seed[0],arg.kmer_length,arg.num_hash_func,arg.bit_vector_width) == true) {
        if(kmer_new[funcm.kmer_length + it_alter]){
            // if current num_extend = extend_amount
            // running extend_out_left is not needed any more
            if ((num_extend + 1) == extend_amount) {
                extension_success = true;
                break;
            }
            else {
                // trace  this kmer recursively
                extend_out_left(
                        kmer_new,
                        num_extend + 1,
                        extend_amount,
                        extension_success,
                        bit_vector,
                        hash_seed,
                        arg,
                        funcm
                        );
            }
        }
    }
    //	free(kmer_new);
    funcm.kmer_ind --;
}


//----------------------------------------------------------------------
// correct_errors_5_prime_end
//----------------------------------------------------------------------
__ONMIC__ inline void correct_errors_5_prime_end(const char* org_sequence, char* sequence, const char* quality_score,
        const uint32_t& index_start, char* sequence_modification, uint32_t& trim_5_end, uint32_t& trim_3_end,
        const uint32_t& org_boundary, const uint32_t& read_length, const uint32_t& max_trimmed_bases,
        const uint32_t& min_check_length, uint32_t& num_corrected_errors_local,
        const uint8_t* bit_vector, const uint32_t* hash_seed, init_args& arg, func_mem& funcm) {
    // |  non-solid  | 1st solid region
    // |--------------------------------------| read
    //         |-----|                          (index_start)-th k-mer
    //--------------------------------------------------
    // list of candidate paths
    vector_candidate_path candidate_path_vector_tmp;
    init_vector_candidate_path(candidate_path_vector_tmp);


    //   // make an initial k-mer
    //   std::string kmer_initial(sequence.substr(index_start, kmer_length));

    //	char* kmer_initial = (char*) malloc(arg.kmer_length + 2);
    //	memcpy(kmer_initial, sequence + index_start, arg.kmer_length);
    char* kmer_initial = funcm.kmer + funcm.kmer_ind * funcm.kmer_length;
    copy_kmer(funcm, sequence + index_start, arg.kmer_length);
    // each alternative neocletide
    bool run_exploration(true);
    path_query_text(kmer_initial, bit_vector,hash_seed,  arg.kmer_length, arg.num_hash_func, arg.bit_vector_width, 0);
    for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
        // make a change
        kmer_initial[0] = NEOCLEOTIDE[it_alter];

        // kmer_initial is solid
//        if ( query_text(kmer_initial, bit_vector, &hash_seed[0],arg.kmer_length,arg.num_hash_func,arg.bit_vector_width)  == true) {
        if(kmer_initial[funcm.kmer_length + it_alter]){
            // generate a new path
            C_candidate_path candidate_path;
            candidate_path.size = 0;
            candidate_path.sum_qs = 0;

            if (sequence[index_start] != NEOCLEOTIDE[it_alter]) {
                pair_uint_char pair_tmp;
                pair_tmp.first  = index_start;
                pair_tmp.second = NEOCLEOTIDE[it_alter];

                push_back(candidate_path,pair_tmp);
                //				candidate_path.modified_bases.
            }

            // if this k-mer is the first k-mer in a read
            // running extend_a_kmer_5_prime_end is not needed any more
            if (index_start == 0) {
                push_back(candidate_path_vector_tmp,candidate_path);
            }
            else if (index_start > 0) {

                // trace  this kmer recursively and update candidate_path_vector_tmp
                extend_a_kmer_5_prime_end(
                        kmer_initial,
                        sequence,
                        index_start,
                        candidate_path,
                        candidate_path_vector_tmp,
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

    // complete exploration was not done because there are too many candidata paths
    // remove all the paths in candidate_path_vector_tmp
    if (run_exploration == false) {
        //		candidate_path_vector_tmp.clear();
        clear(candidate_path_vector_tmp);
    }

    vector_candidate_path candidate_path_vector_tmp_tmp;
    init_vector_candidate_path(candidate_path_vector_tmp_tmp);

    // check the solidness of the leftmost k-mers of each modified base
    bool all_solid_wo_modification(false);

    // each candidate path
    for (C_candidate_path* it_path = candidate_path_vector_tmp.c_path; it_path != candidate_path_vector_tmp.c_path + candidate_path_vector_tmp.size; it_path++) {
        if ((*it_path).size == 0) {
            all_solid_wo_modification = true;
            break;
        }
        else {
            //         std::string sequence_tmp(sequence);
            char* sequence_tmp;//= (char*)malloc(READ_MAX_LEN);
            init_sequence(funcm, &sequence_tmp);
            //			strcpy(sequence_tmp, sequence);
            memcpy(sequence_tmp,sequence, READ_MAX_LEN);

            // index_smallest_modified
            uint32_t index_smallest_modified((*it_path).modified_bases[(*it_path).size - 1].first);

            // number of bases that should be extended
            uint32_t extend_amount;

            // calculate extend_amount
            // no extension is needed
            // kmer_length = 11, max_extension = 5
            // |0|0|0|0|0|0|0|0|0|0|1|1|1|-
            // |0|1|2|3|4|5|6|7|8|9|0|1|2|-
            // |<------------------->|      k = 11
            // |--------------------------- read
            //                     |<------ index_smallest_modified >= 10
            if (index_smallest_modified >= arg.kmer_length - 1) {
                push_back(candidate_path_vector_tmp_tmp, *it_path);
            }
            // extension is needed
            else {
                // applied the modified bases to sequence_tmp
                for (unsigned int it_base = 0; it_base < (*it_path).size; it_base++) {
                    // modify sequence_tmp
                    sequence_tmp[(*it_path).modified_bases[it_base].first] = (*it_path).modified_bases[it_base].second;
                }

                // determine the number of extensions
                // extension amount = kmer_length - index_smallest_modified - 1
                // kmer_length = 11, max_extension = 5
                // |0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|1|1|1|-
                // |5|4|3|2|1|0|1|2|3|4|5|6|7|8|9|0|1|2|-
                //           |<------------------->|      k = 11
                //           |--------------------------- read
                //                     |<------->|        (index_smallest_modified < 10) AND (index_smallest_modified >= 5)
                //     |<------------------->|            index_smallest_modified = 7 -> extend_amount = 3
                if (index_smallest_modified >= arg.kmer_length - arg.max_extension - 1) {
                    extend_amount = arg.kmer_length - index_smallest_modified - 1;
                }
                // kmer_length = 11, max_extension = 5
                // |0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|1|1|1|-
                // |5|4|3|2|1|0|1|2|3|4|5|6|7|8|9|0|1|2|-
                //           |<------------------->|      k = 11
                //           |--------------------------- read
                //           |<------->|                  index_smallest_modified < 5
                else {
                    extend_amount = arg.max_extension;
                }

                bool extension_success = (false);

                //            // generate an initial k-mer
                //            std::string kmer_initial(sequence_tmp.substr(0, kmer_length - 1));
                //            kmer_initial = '0' + kmer_initial;

                //				char* kmer_initial = (char*) malloc(arg.kmer_length + 2);
                //				memcpy(kmer_initial + 1, sequence_tmp , arg.kmer_length -1);
                char* kmer_initial = funcm.kmer + funcm.kmer_ind * funcm.kmer_length;
                copy_kmer(funcm, sequence_tmp , arg.kmer_length -1, 1);
                kmer_initial[0] = 'A';
                path_query_text(kmer_initial, bit_vector,hash_seed,  arg.kmer_length, arg.num_hash_func, arg.bit_vector_width, 0);
                kmer_initial[0] = '0';
                // each alternative neocletide
                for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
                    // make a change
                    kmer_initial[0] = NEOCLEOTIDE[it_alter];

                    // kmer_initial is solid
//                    if (query_text(kmer_initial, bit_vector, &hash_seed[0],arg.kmer_length,arg.num_hash_func,arg.bit_vector_width) == true) {
                        // if extend_amount == 1
                    if(kmer_initial[funcm.kmer_length + it_alter]){
                        // running extend_out_left is not needed any more
                        if (extend_amount == 1) {
                            extension_success = true;
                            break;
                        }
                        else {
                            // trace  this kmer recursively and update candidate_path_vector_tmp
                            extend_out_left(
                                    kmer_initial,
                                    1,
                                    extend_amount,
                                    extension_success,
                                    bit_vector,
                                    hash_seed,
                                    arg,
                                    funcm
                                    );
                        }
                    }
                }

                if (extension_success == true) {
                    push_back(candidate_path_vector_tmp_tmp, *it_path);
                }
                funcm.kmer_ind --;
            }
            funcm.sequence_ind --;
        }
    }

    // remain only really modified paths
    clear(candidate_path_vector_tmp);
    // each path
    for (uint32_t it_candidate = 0; it_candidate < candidate_path_vector_tmp_tmp.size; it_candidate++) {
        // each modification
        bool really_modified(false);
        for (unsigned int it_mod_base = 0; it_mod_base < candidate_path_vector_tmp_tmp.c_path[it_candidate].size; it_mod_base++) {
            if (org_sequence[candidate_path_vector_tmp_tmp.c_path[it_candidate].modified_bases[it_mod_base].first] != candidate_path_vector_tmp_tmp.c_path[it_candidate].modified_bases[it_mod_base].second) {
                really_modified = true;
            }
        }

        if (really_modified) {
            push_back(candidate_path_vector_tmp, candidate_path_vector_tmp_tmp.c_path[it_candidate]);
        }
    }

    //	candidate_path_vector_tmp_tmp = candidate_path_vector_tmp;
    vector_candidate_path_assignment(&candidate_path_vector_tmp_tmp, &candidate_path_vector_tmp);

    // all k-mers are solid without any modification
    // do nothing
    if (all_solid_wo_modification == true) {
    }
    // compare quality scores of candidate paths
    // if the number of paths in candidate_path_vector_tmp_tmp is larger than 1
    else if (candidate_path_vector_tmp_tmp.size > 1) {
        // each path
        C_candidate_path* it_path;
        C_candidate_path* it_path_1st;
        C_candidate_path* it_path_2nd;

        uint32_t qs_1st = (INIT_MIN_QS);
        uint32_t qs_2nd = (INIT_MIN_QS);

        // each candidate path
        for (it_path = candidate_path_vector_tmp_tmp.c_path ; it_path != candidate_path_vector_tmp_tmp.c_path + candidate_path_vector_tmp_tmp.size; it_path++) {
            // each modification
            for (uint32_t it_mod = 0; it_mod < (*it_path).size; it_mod++) {
                // add quality scores of modified bases
                if (sequence[(*it_path).modified_bases[it_mod].first] != (*it_path).modified_bases[it_mod].second) {
                    (*it_path).sum_qs += ((unsigned short int)quality_score[(*it_path).modified_bases[it_mod].first] - arg.quality_score_offset);
                }
            }

            // compare quality scores of each path
            if ((*it_path).sum_qs <= qs_1st) {
                qs_2nd = qs_1st;
                qs_1st = (*it_path).sum_qs;

                it_path_2nd = it_path_1st;
                it_path_1st = it_path;
            }
            else if ((*it_path).sum_qs <= qs_2nd) {
                qs_2nd = (*it_path).sum_qs;

                it_path_2nd = it_path;
            }
        }

        // check whether too many bases are modified in the first path, which implies indel may exist in the original read
        // if too many modifications exist, the read is just trimmed
        //bool keep_going(true);
        bool too_many_corrections(false);

        // at least over MAX_MODIFICATION times of modifications from the first modified position
        if ((((*it_path_1st).modified_bases[(*it_path_1st).size - 1].first + 1) >= min_check_length) &&
                ((*it_path_1st).size > MAX_MODIFICATION)) {
            // each modified base in a right range
            //std::uint32_t partial_qs((*it_path_1st).sum_qs);

            for (uint32_t it_mod_base = ((*it_path_1st).size - 1); it_mod_base >= MAX_MODIFICATION; it_mod_base--) {
                // at least MAX_MODIFICATION times of modifications within min_check_length bases
                if (((*it_path_1st).modified_bases[it_mod_base].first - (*it_path_1st).modified_bases[it_mod_base - MAX_MODIFICATION].first) < min_check_length) {
                    //if (it_mod_base < ((*it_path_1st).modified_bases.size() - 1)) {
                    //   partial_qs -= ((unsigned short int)quality_score[(*it_path_1st).modified_bases[it_mod_base + 1].first] - quality_score_offset);
                    //}

                    // average quality score of the modified bases is too high
                    //if (1.0 * partial_qs / (it_mod_base + 1) > quality_score_cutoff) {
                    // trim_5_end or trim_3_end
                    if ((read_length - (*it_path_1st).modified_bases[it_mod_base].first) <= max_trimmed_bases) {
                        trim_3_end = read_length - (*it_path_1st).modified_bases[it_mod_base].first;
                    }
                    else if (((*it_path_1st).modified_bases[it_mod_base].first + 1) <= max_trimmed_bases) {
                        trim_5_end = (*it_path_1st).modified_bases[it_mod_base].first + 1;
                    }

                    too_many_corrections = true;
                    break;
                    //}
                }
            }
        }

        // use the 1st path if not too many corrections are made AND if the 1st path has a sufficiently low score
        // if an indel exists using quality scores is not a good way to choose the best path
        if ((too_many_corrections == false) && (qs_1st + MIN_QS_DIFF <= qs_2nd)) {
            // update sequence_modification
            for (unsigned int it_base = 0; it_base < (*it_path_1st).size; it_base++) {
                sequence_modification[(*it_path_1st).modified_bases[it_base].first] = (*it_path_1st).modified_bases[it_base].second;

                sequence[(*it_path_1st).modified_bases[it_base].first] = (*it_path_1st).modified_bases[it_base].second;
                num_corrected_errors_local++;
            }
        }
        // hard to choose one path
        else {
            // A AND B
            vector_pair_uint_char base_vector_intersection;
            base_vector_intersection.size = 0;
            base_vector_intersection.sum_qs = 0;
            // A OR B
            vector_pair_uint_char base_vector_union;
            base_vector_union.size = 0;
            base_vector_union.sum_qs = 0;
            // temporary vectors
            vector_pair_uint_char base_vector_intersection_prev;
            memcpy(&base_vector_intersection_prev, &(candidate_path_vector_tmp_tmp.c_path[0]), sizeof(vector_pair_uint_char));
            vector_pair_uint_char base_vector_union_prev; //(candidate_path_vector_tmp_tmp.c_path[0]);
            memcpy(&base_vector_union_prev, &(candidate_path_vector_tmp_tmp.c_path[0]), sizeof(vector_pair_uint_char));

            // each candidate path
            for (uint32_t it_p = 1; it_p < candidate_path_vector_tmp_tmp.size; it_p++) {
                clear(base_vector_intersection);
                clear(base_vector_union);

                base_intersection(base_vector_intersection_prev.uint_char , base_vector_intersection_prev.uint_char + base_vector_intersection_prev.size, candidate_path_vector_tmp_tmp.c_path[it_p].modified_bases,
                        candidate_path_vector_tmp_tmp.c_path[it_p].modified_bases + candidate_path_vector_tmp_tmp.c_path[it_p].size, base_vector_intersection);
                base_union       (base_vector_union_prev.uint_char,        base_vector_union_prev.uint_char + base_vector_union_prev.size,        candidate_path_vector_tmp_tmp.c_path[it_p].modified_bases,
                        candidate_path_vector_tmp_tmp.c_path[it_p].modified_bases + candidate_path_vector_tmp_tmp.c_path[it_p].size, base_vector_union);

                //				base_vector_intersection_prev = base_vector_intersection;
                //				base_vector_union_prev        = base_vector_union;
                memcpy(&base_vector_intersection_prev, &base_vector_intersection, sizeof(vector_pair_uint_char));
                memcpy(&base_vector_union_prev, &base_vector_union, sizeof(vector_pair_uint_char));
            }

            // A - B
            vector_pair_uint_char base_vector_difference;
            base_vector_difference.size = 0;
            base_vector_difference.sum_qs = 0;
            base_difference(base_vector_union.uint_char, base_vector_union.uint_char + base_vector_union.size , base_vector_intersection.uint_char,base_vector_intersection.uint_char + base_vector_intersection.size, base_vector_difference);

            // find trimmed region
            // correcting the 5'-end and 3'-end was not done
            // therefore the total number of trimmed bases is 0 yet
            if (base_vector_difference.size > 0) {
                uint32_t vector_index_leftmost = (0);
                uint32_t vector_index_rightmost = (base_vector_difference.size - 1);

                bool keep_going(true);

                while (keep_going) {
                    // # of trimmed bases at the 5' end: base_vector_difference[vector_index_leftmost].first + 1
                    // # of trimmed bases at the 3' end: read_length - base_vector_difference[vector_index_rightmost].first
                    // the 5'-end is smaller
                    if ((base_vector_difference.uint_char[vector_index_leftmost].first + 1) < (read_length - base_vector_difference.uint_char[vector_index_rightmost].first)) {
                        // check the total number of trimmed bases
                        if ((base_vector_difference.uint_char[vector_index_leftmost].first + 1 + trim_3_end) <= max_trimmed_bases) {
                            if ((base_vector_difference.uint_char[vector_index_leftmost].first + 1) > trim_5_end) {
                                trim_5_end = base_vector_difference.uint_char[vector_index_leftmost].first + 1;
                            }

                            // two points are met
                            if (vector_index_leftmost == vector_index_rightmost) {
                                keep_going = false;
                            }
                            else {
                                vector_index_leftmost++;
                            }
                        }
                        // no need for more check
                        else {
                            keep_going = false;
                        }
                    }
                    // the 3'-end is smaller
                    else {
                        // check the total number of trimmed bases
                        if ((read_length - base_vector_difference.uint_char[vector_index_rightmost].first) <= max_trimmed_bases) {
                            if ((read_length - base_vector_difference.uint_char[vector_index_rightmost].first) > trim_3_end) {
                                trim_3_end = read_length - base_vector_difference.uint_char[vector_index_rightmost].first;
                            }

                            // two points are met
                            if (vector_index_leftmost == vector_index_rightmost) {
                                keep_going = false;
                            }
                            else {
                                vector_index_rightmost--;
                            }
                        }
                        // no need for more check
                        else {
                            keep_going = false;
                        }
                    }
                }
            }

            // find consensus modifications
            for (uint32_t it_inter = 0; it_inter < base_vector_intersection.size; it_inter++) {
                // check whether the base is not in the trimmed regions
                if ((base_vector_intersection.uint_char[it_inter].first <  (read_length - trim_3_end)) &&
                        (base_vector_intersection.uint_char[it_inter].first >= trim_5_end)) {
                    // filter out the bases that are equal to the original ones
                    if (sequence[base_vector_intersection.uint_char[it_inter].first] != base_vector_intersection.uint_char[it_inter].second) {
                        sequence_modification[base_vector_intersection.uint_char[it_inter].first] = base_vector_intersection.uint_char[it_inter].second;
                        sequence[base_vector_intersection.uint_char[it_inter].first] = base_vector_intersection.uint_char[it_inter].second;

                        num_corrected_errors_local++;
                    }
                }
            }
        }
    }
    // only one path
    else if (candidate_path_vector_tmp_tmp.size == 1) {
        // check whether too many bases are modified in the first path, which implies indel may exist in the original read
        // if too many modifications exist, the read is just trimmed
        //bool keep_going(true);
        bool too_many_corrections(false);

        // at least over MAX_MODIFICATION times of modifications from the first modified position
        if (((candidate_path_vector_tmp_tmp.c_path[0].modified_bases[0].first + 1) >= min_check_length) &&
                (candidate_path_vector_tmp_tmp.c_path[0].size > MAX_MODIFICATION)) {
            // calculate sum_qs
            //for (unsigned int it_mod_base = 0; it_mod_base < candidate_path_vector_tmp_tmp[0].modified_bases.size(); it_mod_base++) {
            //   candidate_path_vector_tmp_tmp[0].sum_qs += ((unsigned short int)quality_score[candidate_path_vector_tmp_tmp[0].modified_bases[it_mod_base].first] - quality_score_offset);
            //}

            // each modified base in a right range
            //std::uint32_t partial_qs(candidate_path_vector_tmp_tmp[0].sum_qs);

            for (uint32_t it_mod_base = (candidate_path_vector_tmp_tmp.c_path[0].size - 1); it_mod_base >= MAX_MODIFICATION; it_mod_base--) {
                // at least MAX_MODIFICATION times of modifications within min_check_length bases
                if ((candidate_path_vector_tmp_tmp.c_path[0].modified_bases[it_mod_base].first - candidate_path_vector_tmp_tmp.c_path[0].modified_bases[it_mod_base - MAX_MODIFICATION].first) < min_check_length) {
                    //if (it_mod_base < (candidate_path_vector_tmp_tmp[0].modified_bases.size() - 1)) {
                    //   partial_qs -= ((unsigned short int)quality_score[candidate_path_vector_tmp_tmp[0].modified_bases[it_mod_base + 1].first] - quality_score_offset);
                    //}

                    // average quality score of the modified bases is too high
                    //if (1.0 * partial_qs / (it_mod_base + 1) > quality_score_cutoff) {
                    // trim_5_end or trim_3_end
                    if ((read_length - candidate_path_vector_tmp_tmp.c_path[0].modified_bases[it_mod_base].first) <= max_trimmed_bases) {
                        trim_3_end = read_length - candidate_path_vector_tmp_tmp.c_path[0].modified_bases[it_mod_base].first;

                        // update sequence_modification for the non-trimmed corrections
                        if (it_mod_base > 0) {
                            for (unsigned int it_base = 0; it_base < (it_mod_base - 1); it_base++) {
                                sequence_modification[(candidate_path_vector_tmp_tmp.c_path[0]).modified_bases[it_base].first] = (candidate_path_vector_tmp_tmp.c_path[0]).modified_bases[it_base].second;
                                sequence[(candidate_path_vector_tmp_tmp.c_path[0]).modified_bases[it_base].first] = (candidate_path_vector_tmp_tmp.c_path[0]).modified_bases[it_base].second;
                                num_corrected_errors_local++;
                            }
                        }
                    }
                    else if ((candidate_path_vector_tmp_tmp.c_path[0].modified_bases[it_mod_base].first + 1) <= max_trimmed_bases) {
                        trim_5_end = candidate_path_vector_tmp_tmp.c_path[0].modified_bases[it_mod_base].first + 1;

                        // update sequence_modification for the non-trimmed corrections
                        if (it_mod_base < (candidate_path_vector_tmp_tmp.c_path[0].size - 1)) {
                            for (unsigned int it_base = (it_mod_base + 1); it_base < candidate_path_vector_tmp_tmp.c_path[0].size; it_base++) {
                                sequence_modification[(candidate_path_vector_tmp_tmp.c_path[0]).modified_bases[it_base].first] = (candidate_path_vector_tmp_tmp.c_path[0]).modified_bases[it_base].second;
                                sequence[(candidate_path_vector_tmp_tmp.c_path[0]).modified_bases[it_base].first] = (candidate_path_vector_tmp_tmp.c_path[0]).modified_bases[it_base].second;
                                num_corrected_errors_local++;
                            }
                        }
                    }

                    too_many_corrections = true;
                    break;
                    //}
                }
            }
        }

        // not too many modified bases
        if (too_many_corrections == false) {
            // update sequence_modification
            for (unsigned int it_base = 0; it_base < candidate_path_vector_tmp_tmp.c_path[0].size; it_base++) {
                sequence_modification[candidate_path_vector_tmp_tmp.c_path[0].modified_bases[it_base].first] = candidate_path_vector_tmp_tmp.c_path[0].modified_bases[it_base].second;

                sequence[candidate_path_vector_tmp_tmp.c_path[0].modified_bases[it_base].first] = candidate_path_vector_tmp_tmp.c_path[0].modified_bases[it_base].second;

                num_corrected_errors_local++;
            }
        }
    }
    //	free(kmer_initial);
    funcm.kmer_ind --;
    free(candidate_path_vector_tmp.c_path);
    free(candidate_path_vector_tmp_tmp.c_path);
}

//----------------------------------------------------------------------
// extend_a_kmer_3_prime_end
//----------------------------------------------------------------------
__ONMIC__  void extend_a_kmer_3_prime_end(const char* kmer, const char* sequence, const uint32_t& index_kmer,
        C_candidate_path& current_path, vector_candidate_path& candidate_path_vector, const uint32_t& org_boundary,
        const char* quality_score, const unsigned char* bit_vector, const uint32_t* hash_seed, bool& run_exploration, uint32_t read_length,init_args& arg, func_mem& funcm) {


    if (run_exploration == true) {
        //		int tt;
        //		printf("extend_a_kmer_3_prime_end\n");
        //		// generate a new k-mer
        //		//      std::string kmer_new(kmer.substr(1, kmer_length - 1));
        //		//      kmer_new = kmer_new + sequence[index_kmer + kmer_length];
        //		char* kmer_new = (char*) malloc(arg.kmer_length + 2);
        //		memcpy(kmer_new, kmer + 1, arg.kmer_length);
        //		kmer_new[arg.kmer_length - 1] = sequence[index_kmer + arg.kmer_length];
        char* kmer_new = funcm.kmer + funcm.kmer_ind * funcm.kmer_length;
        copy_kmer(funcm, kmer + 1, arg.kmer_length - 1);
        kmer_new[arg.kmer_length - 1] = sequence[index_kmer + arg.kmer_length];
        //		printf("kmer %s \n", kmer_new);
        //		printf("sequence %s \n", sequence);
        //		printf("index_kmer %d \n", index_kmer);
        //		printf("quality_score %s \n", quality_score);
        //		for(int i = 0; i < current_path.size; i++){
        //			printf("current_path %d %c \n", current_path.modified_bases[i].first, current_path.modified_bases[i].second);
        //		}
        //		printf("org_boundary %d \n", org_boundary);
        //		printf("quality_score %d \n", ((unsigned short int)quality_score[index_kmer + arg.kmer_length] - arg.quality_score_offset));
        //		const char* sequence, const uint32_t& index_kmer,
        //				C_candidate_path& current_path, vector_candidate_path& candidate_path_vector, const uint32_t& org_boundary,
        //printf(" read %d %d \n",  arg.quality_score_offset, arg.extremely_low_quality_score);

        //		// this was the real boundary between solid k-mers and weak k-mers
        //		// check all possible cases
        if ((index_kmer == (org_boundary - 1)) || (((unsigned short int)quality_score[index_kmer + arg.kmer_length] - arg.quality_score_offset) <= arg.extremely_low_quality_score)) {
            // each alternative neocletide
        	path_query_text(kmer_new, bit_vector,hash_seed,  arg.kmer_length, arg.num_hash_func, arg.bit_vector_width, 1);
            for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
                // make a change
                kmer_new[arg.kmer_length - 1] = NEOCLEOTIDE[it_alter];

                // kmer_new is solid
//                if (query_text(kmer_new, bit_vector, &hash_seed[0],arg.kmer_length,arg.num_hash_func,arg.bit_vector_width) == true) {
                if(kmer_new[funcm.kmer_length + it_alter]){
                    // generate a new path
                    //					C_candidate_path temporary_path;
                    //					memcpy(&temporary_path, &current_path, sizeof(C_candidate_path));

                    int increment = 0;
                    // not equal to the original character
                    if (sequence[index_kmer + arg.kmer_length] != NEOCLEOTIDE[it_alter]) {
                        pair_uint_char pair_tmp;
                        pair_tmp.first  = index_kmer + arg.kmer_length;
                        pair_tmp.second = NEOCLEOTIDE[it_alter];

                        push_back(current_path, pair_tmp);

                        increment ++;
                    }

                    // if this k-mer is the last k-mer in a read
                    // running extend_a_kmer_3_prime_end is not needed any more
                    if ((index_kmer + 1) == (read_length - arg.kmer_length)) {
                        push_back(candidate_path_vector, current_path);
                        // too many candidate paths
                        if (candidate_path_vector.size > MAX_CANDIDATE_PATHS) {

                            run_exploration = false;
                        }
                    }
                    else if ((index_kmer + 1) < (read_length - arg.kmer_length)) {
                        //						scanf("%d",&tt);
                        //						printf("1\n");
                        // trace  this kmer recursively and update candidate_path_vector
                        extend_a_kmer_3_prime_end(
                                kmer_new,
                                sequence,
                                index_kmer + 1,
                                current_path,
                                candidate_path_vector,
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
                    current_path.size -= increment;
                }
            }
        }
        else {

            // kmer_new is a solid k-mer
            if (query_text(kmer_new, bit_vector, &hash_seed[0],arg.kmer_length,arg.num_hash_func,arg.bit_vector_width) == true) {
                // if this k-mer is the last k-mer in a read
                // running extend_a_kmer_3_prime_end is not needed any more
                if ((index_kmer + 1) == (read_length - arg.kmer_length)) {
                    push_back(candidate_path_vector, current_path);
                    // too many candidate paths
                    if (candidate_path_vector.size > MAX_CANDIDATE_PATHS) {

                        run_exploration = false;
                    }
                }
                else if ((index_kmer + 1) < (read_length - arg.kmer_length)) {
                    //					scanf("%d",&tt);
                    //					printf("2\n");
                    extend_a_kmer_3_prime_end(
                            kmer_new,
                            sequence,
                            index_kmer + 1,
                            current_path,
                            candidate_path_vector,
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
            else {
            	path_query_text(kmer_new, bit_vector,hash_seed,  arg.kmer_length, arg.num_hash_func, arg.bit_vector_width, 1);
                // each alternative neocletide
                for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
                    // not equal to the original character
                    if (sequence[index_kmer + arg.kmer_length] != NEOCLEOTIDE[it_alter]) {
                        // make a change
                        kmer_new[arg.kmer_length - 1] = NEOCLEOTIDE[it_alter];

                        // kmer_new is solid
//                        if (query_text(kmer_new, bit_vector, &hash_seed[0],arg.kmer_length,arg.num_hash_func,arg.bit_vector_width) == true) {
                        if(kmer_new[funcm.kmer_length + it_alter]){
                        // generate a new path
                            //							C_candidate_path temporary_path(current_path);
                            //							C_candidate_path temporary_path;
                            //							memcpy(&temporary_path, &current_path, sizeof(C_candidate_path));

                            pair_uint_char pair_tmp;
                            pair_tmp.first  = index_kmer + arg.kmer_length;
                            pair_tmp.second = NEOCLEOTIDE[it_alter];

                            push_back(current_path, pair_tmp);

                            // if this k-mer is the last k-mer in a read
                            // running extend_a_kmer_3_prime_end is not needed any more
                            if ((index_kmer + 1) == (read_length - arg.kmer_length)) {
                                push_back(candidate_path_vector, current_path);
                                // too many candidate paths
                                if (candidate_path_vector.size > MAX_CANDIDATE_PATHS) {

                                    run_exploration = false;
                                }
                            }
                            else if ((index_kmer + 1) < (read_length - arg.kmer_length)) {
                                //								scanf("%d",&tt);
                                // trace  this kmer recursively and update candidate_path_vector
                                //								printf("3\n");
                                extend_a_kmer_3_prime_end(
                                        kmer_new,
                                        sequence,
                                        index_kmer + 1,
                                        current_path,
                                        candidate_path_vector,
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
                            current_path.size --;
                        }
                    }
                }
            }
        }
        //		free(kmer_new);
        funcm.kmer_ind --;
    }

}

//----------------------------------------------------------------------
// check_first_kmer_sorted
//----------------------------------------------------------------------
__ONMIC__ inline void check_first_kmer_sorted(const char*  kmer, C_candidate_path& candidate_path_in,
        const uint32_t* low_qs_indexes_sorted, vector_candidate_path& candidate_path_vector,
        const uint32_t& index, const unsigned char* bit_vector, const uint32_t*  hash_seed, init_args& arg, func_mem& funcm) {
    //   std::string new_kmer(kmer);
    //	char* new_kmer = (char*) malloc(arg.kmer_length + 2);
    //	memcpy(new_kmer, kmer, arg.kmer_length);
    char* new_kmer = funcm.kmer + funcm.kmer_ind * funcm.kmer_length;
    copy_kmer(funcm, kmer, arg.kmer_length);
    for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
        // make a new k-mer
        new_kmer[low_qs_indexes_sorted[index]] = NEOCLEOTIDE[it_alter];

        //		C_candidate_path candidate_path_next;
        //		memcpy(&candidate_path_next, &candidate_path_in, sizeof(C_candidate_path));

        int increment = 0;
        if (kmer[low_qs_indexes_sorted[index]] != NEOCLEOTIDE[it_alter]) {
            pair_uint_char pair_tmp;
            pair_tmp.first = low_qs_indexes_sorted[index];
            pair_tmp.second = NEOCLEOTIDE[it_alter];
            push_back(candidate_path_in, pair_tmp);
            increment ++;
        }

        // the max number of low quality base is reached
        // add the path to the vector
        // the original k-mer is also included
        if (index == MAX_LOW_QS_BASES - 1) {
            if (query_text(new_kmer, bit_vector, &hash_seed[0],arg.kmer_length,arg.num_hash_func,arg.bit_vector_width) == true) {
                push_back(candidate_path_vector, candidate_path_in);
            }
        }
        // the rightmost low quality base is not reached
        // do recursively
        else {
            check_first_kmer_sorted(
                    new_kmer,
                    candidate_path_in,
                    low_qs_indexes_sorted,
                    candidate_path_vector,
                    index + 1,
                    bit_vector,
                    hash_seed,
                    arg,
                    funcm
                    );
        }
        candidate_path_in.size -= increment;
    }
    //	free(new_kmer);
    funcm.kmer_ind --;
}


//----------------------------------------------------------------------
// check_first_kmer
//----------------------------------------------------------------------
__ONMIC__ void check_first_kmer(const char* kmer, C_candidate_path& candidate_path_in,
        const uint32_t* low_qs_indexes, vector_candidate_path& candidate_path_vector,
        const uint32_t& index, const unsigned char* bit_vector, const uint32_t*  hash_seed, init_args& arg, func_mem& funcm) {
    //   std::string new_kmer(kmer);
    //	char* new_kmer = (char*) malloc(arg.kmer_length + 2);
    //	memcpy(new_kmer, kmer, arg.kmer_length);

    char* new_kmer = funcm.kmer + funcm.kmer_ind * funcm.kmer_length;
    copy_kmer(funcm, kmer, arg.kmer_length);

    for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
        // make a new k-mer
        new_kmer[low_qs_indexes[index]] = NEOCLEOTIDE[it_alter];

        //		C_candidate_path candidate_path_next;
        //		memcpy(&candidate_path_next, &candidate_path_in, sizeof(C_candidate_path));

        int increment = 0;
        if (kmer[low_qs_indexes[index]] != NEOCLEOTIDE[it_alter]) {
            pair_uint_char pair_tmp;
            pair_tmp.first = low_qs_indexes[index];
            pair_tmp.second = NEOCLEOTIDE[it_alter];
            push_back(candidate_path_in, pair_tmp);
            increment ++;
        }

        // the rightmost low quality base is reached
        // add the path to the vector
        // the original k-mer is also included
        if (index == arg.kmer_length - 1) {
            if (query_text(new_kmer, bit_vector, &hash_seed[0],arg.kmer_length,arg.num_hash_func,arg.bit_vector_width) == true) {
                push_back(candidate_path_vector,candidate_path_in);
            }
        }
        // the rightmost low quality base is not reached
        // do recursively
        else {
            check_first_kmer(
                    new_kmer,
                    candidate_path_in,
                    low_qs_indexes,
                    candidate_path_vector,
                    index + 1,
                    bit_vector,
                    hash_seed,
                    arg,
                    funcm
                    );
        }
        candidate_path_in.size -= increment;
    }
    //	free(new_kmer);
    funcm.kmer_ind --;
}



//----------------------------------------------------------------------
// correct_errors_first_kmer
//----------------------------------------------------------------------
__ONMIC__ inline void correct_errors_first_kmer(const char*  sequence, const char* quality_score, char* sequence_modification,
        vector_candidate_path& candidate_path_vector, const unsigned char* bit_vector,
        const uint32_t*  hash_seed, init_args& arg, func_mem& funcm) {
    //   std::string first_kmer(sequence.substr(0, kmer_length));

    //	char* first_kmer = (char*) malloc(arg.kmer_length + 2);
    //	memcpy(first_kmer, sequence, arg.kmer_length);
    char* first_kmer = funcm.kmer + funcm.kmer_ind * funcm.kmer_length;
    copy_kmer(funcm, sequence, arg.kmer_length);
    //	std::vector<unsigned int> low_qs_indexes;
    uint32_t* low_qs_indexes = (uint32_t*)malloc(arg.kmer_length * sizeof(uint32_t));
    uint32_t low_qs_indexes_size = 0;
    for (unsigned int it_bases = 0; it_bases < arg.kmer_length; it_bases++) {
        if (((unsigned short int)quality_score[it_bases] - arg.quality_score_offset) < arg.quality_score_cutoff) {
            //			low_qs_indexes.push_back(it_bases);
            low_qs_indexes[low_qs_indexes_size++] = it_bases;
        }
    }

    // low quality bases exist
    if (low_qs_indexes_size == 0) {
        // the first k-mer is solid
        if (query_text(first_kmer, bit_vector, &hash_seed[0],arg.kmer_length,arg.num_hash_func,arg.bit_vector_width) == true) {
            C_candidate_path candidate_path;
            candidate_path.size = 0;
            candidate_path.sum_qs = 0;

            push_back(candidate_path_vector, candidate_path);
        }
        // the first k-mer is not solid
        else {
            // change the bases with the lowest MAX_LOW_QS_BASES quality scores
            uint32_t* low_qs_indexes_sorted = (sort_indexes(quality_score, arg.kmer_length));

            C_candidate_path candidate_path_sorted;
            candidate_path_sorted.size = 0;
            candidate_path_sorted.sum_qs = 0;
            check_first_kmer_sorted(
                    first_kmer,
                    candidate_path_sorted,
                    low_qs_indexes_sorted,
                    candidate_path_vector,
                    0,
                    bit_vector,
                    hash_seed,
                    arg,
                    funcm
                    );
            // change each base in the first k-mer and check whether it is solid
            // it is still needed because there could be more than MAX_LOW_QS_BASES bases that have
            // lower quality scores than real problematic ones
            for (unsigned int it_bases = 0; it_bases < arg.kmer_length; it_bases++) {
                //            std::string kmer_tmp(first_kmer);
                //				char* kmer_tmp = (char*) malloc(arg.kmer_length + 2);
                //				memcpy(kmer_tmp, first_kmer , arg.kmer_length);
                char* kmer_tmp = funcm.kmer + funcm.kmer_ind * funcm.kmer_length;
                copy_kmer(funcm, first_kmer , arg.kmer_length);
                // each alternative neocletide
                for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
                    // not equal to the original character
                    if (first_kmer[it_bases] != NEOCLEOTIDE[it_alter]) {
                        // generate a new k-mer
                        kmer_tmp[it_bases] = NEOCLEOTIDE[it_alter];

                        // add kmer_tmp to candidate_path_tmp if it is solid
                        if (query_text(kmer_tmp, bit_vector, &hash_seed[0],arg.kmer_length,arg.num_hash_func,arg.bit_vector_width) == true) {
                            // generate a new candidate path
                            C_candidate_path candidate_path;
                            candidate_path.size = 0;
                            candidate_path.sum_qs = 0;

                            pair_uint_char pair_tmp;
                            pair_tmp.first = it_bases;
                            pair_tmp.second = NEOCLEOTIDE[it_alter];
                            push_back(candidate_path, pair_tmp);

                            push_back(candidate_path_vector, candidate_path);
                        }
                    }
                }
                //				free(kmer_tmp);
                funcm.kmer_ind --;
            }
        }
    }
    // low quality bases exist
    else if (arg.kmer_length > 0) {
        // the number of low-quality bases is smaller than the threshold
        if (arg.kmer_length <= MAX_LOW_QS_BASES) {
            C_candidate_path candidate_path;
            candidate_path.size = 0;
            candidate_path.sum_qs = 0;
            check_first_kmer(
                    first_kmer,
                    candidate_path,
                    low_qs_indexes,
                    candidate_path_vector,
                    0,
                    bit_vector,
                    hash_seed,
                    arg,
                    funcm
                    );
            // no candidate path is found
            if (candidate_path_vector.size == 0) {
                // the first k-mer is solid
                if (query_text(first_kmer, bit_vector, &hash_seed[0],arg.kmer_length,arg.num_hash_func,arg.bit_vector_width) == true) {
                    clear_path(candidate_path);
                    push_back(candidate_path_vector, candidate_path);
                }
                else {
                    // change each base in the first k-mer and check whether it is solid
                    for (unsigned int it_bases = 0; it_bases < arg.kmer_length; it_bases++) {
                        //                  std::string kmer_tmp(first_kmer);
                        //						char* kmer_tmp = (char*) malloc(arg.kmer_length + 2);
                        //						memcpy(kmer_tmp, first_kmer, arg.kmer_length);

                        char* kmer_tmp = funcm.kmer + funcm.kmer_ind * funcm.kmer_length;
                        copy_kmer(funcm, first_kmer, arg.kmer_length);

                        // each alternative neocletide
                        for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
                            // not equal to the original character
                            if (first_kmer[it_bases] != NEOCLEOTIDE[it_alter]) {
                                // generate a new k-mer
                                kmer_tmp[it_bases] = NEOCLEOTIDE[it_alter];

                                // add kmer_tmp to candidate_path_tmp if it is solid
                                if (query_text(kmer_tmp, bit_vector, &hash_seed[0],arg.kmer_length,arg.num_hash_func,arg.bit_vector_width) == true) {
                                    // generate a new candidate path
                                    clear_path(candidate_path);

                                    pair_uint_char pair_tmp;
                                    pair_tmp.first = it_bases;
                                    pair_tmp.second = NEOCLEOTIDE[it_alter];
                                    push_back(candidate_path, pair_tmp);

                                    push_back(candidate_path_vector, candidate_path);
                                }
                            }
                        }
                        //						free(kmer_tmp);
                        funcm.kmer_ind --;
                    }
                }
            }
        }
        // too many low-quality bases
        else {
            // change the bases with the lowest MAX_LOW_QS_BASES quality scores
            uint32_t* low_qs_indexes_sorted = (sort_indexes(quality_score, arg.kmer_length));
            uint32_t low_qs_indexes_sorted_size = 0;
            C_candidate_path candidate_path_sorted;
            candidate_path_sorted.size = 0;
            candidate_path_sorted.sum_qs = 0;
            check_first_kmer_sorted(
                    first_kmer,
                    candidate_path_sorted,
                    low_qs_indexes_sorted,
                    candidate_path_vector,
                    0,
                    bit_vector,
                    hash_seed,
                    arg,
                    funcm
                    );
            // change each base in the first k-mer and check whether it is solid
            C_candidate_path candidate_path;

            for (unsigned int it_bases = 0; it_bases < arg.kmer_length; it_bases++) {
                //            std::string kmer_tmp(first_kmer);
                //				char* kmer_tmp = (char*) malloc(arg.kmer_length + 2);
                //				memcpy(kmer_tmp, first_kmer , arg.kmer_length);
                char* kmer_tmp = funcm.kmer + funcm.kmer_ind * funcm.kmer_length;
                copy_kmer(funcm, first_kmer, arg.kmer_length);

                // each alternative neocletide
                for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
                    // not equal to the original character
                    if (first_kmer[it_bases] != NEOCLEOTIDE[it_alter]) {
                        // generate a new k-mer
                        kmer_tmp[it_bases] = NEOCLEOTIDE[it_alter];

                        // add kmer_tmp to candidate_path_tmp if it is solid
                        if (query_text(kmer_tmp, bit_vector, &hash_seed[0],arg.kmer_length,arg.num_hash_func,arg.bit_vector_width) == true) {
                            // generate a new candidate path
                            clear_path(candidate_path);

                            pair_uint_char pair_tmp;
                            pair_tmp.first = it_bases;
                            pair_tmp.second = NEOCLEOTIDE[it_alter];
                            push_back(candidate_path,pair_tmp);

                            push_back(candidate_path_vector, candidate_path);
                        }
                    }
                }
                //				free(kmer_tmp);
                funcm.kmer_ind --;
            }

            // the first k-mer is solid
            if (query_text(first_kmer, bit_vector, &hash_seed[0],arg.kmer_length,arg.num_hash_func,arg.bit_vector_width) == true) {
                clear_path(candidate_path);

                push_back(candidate_path_vector, candidate_path);
            }
        }
    }
    //	free(first_kmer);
    funcm.kmer_ind --;
}


//----------------------------------------------------------------------
// extend_out_right
//----------------------------------------------------------------------
__ONMIC__ void extend_out_right(const char* kmer, const uint32_t& num_extend, const uint32_t& extend_amount,
        bool& extension_success, const unsigned char* bit_vector, const uint32_t* hash_seed, init_args& arg, func_mem& funcm) {
    // generate a new k-mer
    //   std::string kmer_new(kmer.substr(1, kmer_length - 1));
    //   kmer_new = kmer_new + '0';
    //	char* kmer_new = (char* ) malloc(arg.kmer_length + 2);
    //	memcpy(kmer_new, kmer + 1, arg.kmer_length - 1);
    //	kmer_new[arg.kmer_length - 1] = '0';
    char* kmer_new = funcm.kmer + funcm.kmer_ind * funcm.kmer_length;
    //	memcpy(kmer_new, kmer + 1, arg.kmer_length - 1);
    copy_kmer(funcm, kmer + 1, arg.kmer_length -1);
    kmer_new[arg.kmer_length - 1] = 'A';
    path_query_text(kmer_new, bit_vector,hash_seed,  arg.kmer_length, arg.num_hash_func, arg.bit_vector_width, 1);
    kmer_new[arg.kmer_length - 1] = '0';
    // each alternative neocletide
    for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
        // generate kmer_new
        kmer_new[arg.kmer_length - 1] = NEOCLEOTIDE[it_alter];

        // kmer_new is solid
//        if (query_text(kmer_new, bit_vector, &hash_seed[0],arg.kmer_length,arg.num_hash_func,arg.bit_vector_width) == true) {
        if(kmer_new[funcm.kmer_length + it_alter]){
            // if current num_extend = extend_amount
            // running extend_out_right is not needed any more
            if ((num_extend + 1) == extend_amount) {
                extension_success = true;
                break;
            }
            else {
                // trace  this kmer recursively
                extend_out_right(
                        kmer_new,
                        num_extend + 1,
                        extend_amount,
                        extension_success,
                        bit_vector,
                        hash_seed,
                        arg,
                        funcm
                        );
            }
        }
    }
    //	free(kmer_new);
    funcm.kmer_ind --;
}

//----------------------------------------------------------------------
// extend_first_kmer_to_right
//----------------------------------------------------------------------
__ONMIC__ void extend_first_kmer_to_right(const char*  sequence, const char*  quality_score, C_candidate_path& candidate_path_in,
        vector_candidate_path& candidate_path_vector_all, const uint32_t& read_length, const unsigned char* bit_vector,
        const uint32_t*  hash_seed, bool& run_exploration, init_args& arg, func_mem& funcm) {
    // generate the first k-mer
    //   std::string first_kmer(sequence.substr(0, kmer_length + 1));
    //	 char* first_kmer = (char*) malloc(arg.kmer_length + 3);
    //
    //	 memcpy(first_kmer, sequence, arg.kmer_length + 1);
    char* first_kmer = funcm.kmer + funcm.kmer_ind * funcm.kmer_length;
    copy_kmer(funcm, sequence, arg.kmer_length + 1);

    for (unsigned int it_base = 0; it_base < candidate_path_in.size; it_base++) {
        first_kmer[candidate_path_in.modified_bases[it_base].first] = candidate_path_in.modified_bases[it_base].second;
    }

    // generate the second k-mer
    //   std::string second_kmer(first_kmer.substr(1, kmer_length));
    //	 char* second_kmer = (char*) malloc(arg.kmer_length + 3);
    //	 memcpy(second_kmer, first_kmer + 1, arg.kmer_length);
    char* second_kmer = funcm.kmer + funcm.kmer_ind * funcm.kmer_length;
    copy_kmer(funcm, first_kmer + 1, arg.kmer_length);

    // list of candidate paths
    vector_candidate_path candidate_path_vector_tmp;
    init_vector_candidate_path(candidate_path_vector_tmp);
    // second_kmer is solid
    if (query_text(second_kmer, bit_vector, &hash_seed[0],arg.kmer_length,arg.num_hash_func,arg.bit_vector_width) == true) {

        // if this k-mer is the last k-mer in a read
        // running extend_a_kmer_3_prime_end is not needed any more
        if ((read_length - arg.kmer_length) == 1) {
            push_back(candidate_path_vector_tmp, candidate_path_in);
        }
        else if ((read_length - arg.kmer_length) > 1) {
            // trace  this kmer recursively and update candidate_path_vector_tmp
            // org_boundary is not needed in this case: read_length is being used as a dummy number
            extend_a_kmer_3_prime_end(
                    second_kmer,
                    sequence,
                    1,
                    candidate_path_in,
                    candidate_path_vector_tmp,
                    read_length,
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
    // second_kmer is not solid
    else {
    	path_query_text(second_kmer, bit_vector,hash_seed,  arg.kmer_length, arg.num_hash_func, arg.bit_vector_width, 1);
        // each alternative neocletide
        for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
            // not equal to the original character
            if (sequence[arg.kmer_length] != NEOCLEOTIDE[it_alter]) {
                // make a change
                second_kmer[arg.kmer_length - 1] = NEOCLEOTIDE[it_alter];

                // new second_kmer is solid
//                if (query_text(second_kmer, bit_vector, &hash_seed[0],arg.kmer_length,arg.num_hash_func,arg.bit_vector_width) == true) {
                if(second_kmer[funcm.kmer_length + it_alter]){
                    // generate a new path
                    C_candidate_path new_candidate_path;
                    memcpy(&new_candidate_path, &candidate_path_in, sizeof(C_candidate_path));

                    pair_uint_char pair_tmp;
                    pair_tmp.first  = arg.kmer_length;
                    pair_tmp.second = NEOCLEOTIDE[it_alter];

                    push_back(new_candidate_path, pair_tmp);

                    // if this k-mer is the last k-mer in a read
                    // running extend_a_kmer_3_prime_end is not needed any more
                    if ((read_length - arg.kmer_length) == 1) {
                        push_back(candidate_path_vector_tmp, new_candidate_path);
                    }
                    else if ((read_length - arg.kmer_length) > 1) {
                        // trace  this kmer recursively and update candidate_path_vector_tmp
                        // org_boundary is not needed in this case: read_length is being used as a dummy number
                        extend_a_kmer_3_prime_end(
                                second_kmer,
                                sequence,
                                1,
                                new_candidate_path,
                                candidate_path_vector_tmp,
                                read_length,
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
        }
    }

    // complete exploration was not done because there are too many candidata paths
    // remove all the paths in candidate_path_vector_tmp
    if (run_exploration == false) {
        //      candidate_path_vector_tmp.clear();
        clear(candidate_path_vector_tmp);
    }
    char* sequence_tmp;//  = (char*) malloc(READ_MAX_LEN);
    init_sequence(funcm, &sequence_tmp);
    // each candidate path
    for (C_candidate_path* it_path = candidate_path_vector_tmp.c_path; it_path != candidate_path_vector_tmp.c_path + candidate_path_vector_tmp.size; it_path++) {
        //      std::string sequence_tmp(sequence);

        memcpy(sequence_tmp, sequence, read_length);

        // index_largest_modified
        uint32_t index_largest_modified = ((*it_path).modified_bases[(*it_path).size - 1].first);

        // number of bases that should be extended
        uint32_t extend_amount;

        // calculate extend_amount
        // no extension is needed
        // sequence.length() = 20, kmer_length = 11, max_extension = 5
        // |0|0|0|1|1|1|1|1|1|1|1|1|1|
        // |7|8|9|0|1|2|3|4|5|6|7|8|9|
        //     |<------------------->| k = 11
        // --------------------------| read
        // ----->|                     index_largest_modified <= 9
        if (index_largest_modified <= read_length - arg.kmer_length) {
            push_back(candidate_path_vector_all, *it_path);
        }
        // extension is needed
        else {
            // applied the modified bases to sequence_tmp
            for (unsigned int it_base = 0; it_base < (*it_path).size; it_base++) {
                // modify sequence_tmp
                sequence_tmp[(*it_path).modified_bases[it_base].first] = (*it_path).modified_bases[it_base].second;
            }

            // determine the number of extensions
            // sequence.length() = 20, kmer_length = 11, max_extension = 5
            // |0|0|0|1|1|1|1|1|1|1|1|1|1|0|0|0|0|0|
            // |7|8|9|0|1|2|3|4|5|6|7|8|9|1|2|3|4|5|
            //     |<------------------->|           k = 11
            // --------------------------|           read
            //       |<------->|                     (index_largest_modified > 10) AND (index_largest_modified <= 14)
            //           |<------------------->|     index_largest_modified = 12 -> extend_amout = 3
            if (index_largest_modified <= read_length + arg.max_extension - arg.kmer_length) {
                extend_amount = arg.kmer_length - (read_length - index_largest_modified);
            }
            // sequence.length() = 20, kmer_length = 11, max_extension = 5
            // |0|0|0|1|1|1|1|1|1|1|1|1|1|0|0|0|0|0|
            // |7|8|9|0|1|2|3|4|5|6|7|8|9|1|2|3|4|5|
            //     |<------------------->|           k = 11
            // --------------------------|           read
            //                 |<------->|           index_largest_modified > 15
            else {
                extend_amount = arg.max_extension;
            }

            bool extension_success(false);

            // generate an initial k-mer
            // sequence.length() = 20, kmer_length = 11
            // |0|0|0|1|1|1|1|1|1|1|1|1|1|
            // |7|8|9|0|1|2|3|4|5|6|7|8|9|
            //       |<----------------->| kmer_length - 1 = 10
            // --------------------------| read
            //       |-|                   20 - 11 + 1 = 10
            //         std::string kmer_initial(sequence_tmp.substr(sequence.length() - kmer_length + 1, kmer_length - 1));
            //         kmer_initial = kmer_initial + '0';

            //         char* kmer_initial = (char*) malloc(arg.kmer_length + 2);
            //         memcpy(kmer_initial, sequence_tmp + read_length- arg.kmer_length + 1, arg.kmer_length - 1);
            //         kmer_initial[arg.kmer_length - 1] = '0';
            char* kmer_initial = funcm.kmer + funcm.kmer_ind * funcm.kmer_length;
            copy_kmer(funcm, sequence_tmp + read_length- arg.kmer_length + 1, arg.kmer_length - 1);
            kmer_initial[arg.kmer_length - 1] = 'A';
        	path_query_text(kmer_initial, bit_vector,hash_seed,  arg.kmer_length, arg.num_hash_func, arg.bit_vector_width, 1);
            kmer_initial[arg.kmer_length - 1] = '0';
            // each alternative neocletide
            for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
                // make a change
                kmer_initial[arg.kmer_length - 1] = NEOCLEOTIDE[it_alter];
                // kmer_initial is solid
//                if (query_text(kmer_initial, bit_vector, &hash_seed[0],arg.kmer_length,arg.num_hash_func,arg.bit_vector_width) == true) {
                if(kmer_initial[funcm.kmer_length + it_alter]){
                    // if extend_amount == 1
                    // running extend_out_right is not needed any more
                    if (extend_amount == 1) {
                        extension_success = true;
                        break;
                    }
                    else {
                        // trace  this kmer recursively and update candidate_path_vector_tmp
                        extend_out_right(
                                kmer_initial,
                                1,
                                extend_amount,
                                extension_success,
                                bit_vector,
                                hash_seed,
                                arg,
                                funcm
                                );
                    }
                }
            }

            if (extension_success == true) {
                push_back(candidate_path_vector_all, *it_path);
            }
            funcm.kmer_ind --;
        }
    }
    //   free(sequence_tmp);
    //   free(first_kmer);
    //   free(second_kmer);
    funcm.kmer_ind --;
    funcm.kmer_ind --;
    funcm.sequence_ind --;
    free(candidate_path_vector_tmp.c_path);
}

//----------------------------------------------------------------------
// solid_first_kmer
//----------------------------------------------------------------------
__ONMIC__ inline void solid_first_kmer(const C_candidate_path& candidate_path, const char* sequence, bool& extension_success,
        const unsigned char* bit_vector, const uint32_t* hash_seed, init_args& arg, func_mem& funcm) {
    // index_smallest_modified
    uint32_t index_smallest_modified = (candidate_path.modified_bases[0].first);

    // number of bases that should be extended
    uint32_t extend_amount;

    // applied the modified bases to first_kmer
    //   std::string first_kmer(sequence.substr(0, kmer_length));

    //	char* first_kmer = (char*) malloc(arg.kmer_length + 2);
    //	memcpy(first_kmer, sequence , arg.kmer_length);
    char* first_kmer = funcm.kmer + funcm.kmer_ind * funcm.kmer_length;
    copy_kmer(funcm, sequence , arg.kmer_length);

    for (unsigned int it_base = 0; it_base < candidate_path.size; it_base++) {
        first_kmer[candidate_path.modified_bases[it_base].first] = candidate_path.modified_bases[it_base].second;
    }

    // determine the number of extensions
    // extension amount = kmer_length - index_smallest_modified - 1
    // kmer_length = 11, max_extension = 5
    // |0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|1|1|1|-
    // |5|4|3|2|1|0|1|2|3|4|5|6|7|8|9|0|1|2|-
    //           |<------------------->|      k = 11
    //           |--------------------------- read
    //                     |<------->|        (index_smallest_modified < 10) AND (index_smallest_modified >= 5)
    //     |<------------------->|            index_smallest_modified = 7 -> extend_amount = 3
    if (index_smallest_modified >= arg.kmer_length - arg.max_extension - 1) {
        extend_amount = arg.kmer_length - index_smallest_modified - 1;
    }
    // kmer_length = 11, max_extension = 5
    // |0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|1|1|1|-
    // |5|4|3|2|1|0|1|2|3|4|5|6|7|8|9|0|1|2|-
    //           |<------------------->|      k = 11
    //           |--------------------------- read
    //           |<------->|                  index_smallest_modified < 5
    else {
        extend_amount = arg.max_extension;
    }

    // generate an initial k-mer
    //   std::string kmer_initial(first_kmer.substr(0, kmer_length - 1));
    //   kmer_initial = '0' + kmer_initial;

    //	char* kmer_initial = (char*) malloc(arg.kmer_length + 2);
    //	memcpy(kmer_initial + 1, first_kmer, arg.kmer_length - 1);
    char* kmer_initial = funcm.kmer + funcm.kmer_ind * funcm.kmer_length;;
    copy_kmer(funcm, first_kmer, arg.kmer_length - 1, 1);
    kmer_initial[0] = '0';
    // each alternative neocletide
    for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
        // make a change
        kmer_initial[0] = NEOCLEOTIDE[it_alter];

        // kmer_initial is solid
        if (query_text(kmer_initial, bit_vector, &hash_seed[0],arg.kmer_length,arg.num_hash_func,arg.bit_vector_width) == true) {
            // if extend_amount == 1
            // running extend_out_left is not needed any more
            if (extend_amount == 1) {
                extension_success = true;
                break;
            }
            else {
                // trace  this kmer recursively and update candidate_path_vector_tmp
                extend_out_left(
                        kmer_initial,
                        1,
                        extend_amount,
                        extension_success,
                        bit_vector,
                        hash_seed,
                        arg,
                        funcm
                        );
            }
        }
    }
    //	free(first_kmer);
    //	free(kmer_initial);
    funcm.kmer_ind --;
    funcm.kmer_ind --;
}


//----------------------------------------------------------------------
// correct_errors_3_prime_end
//----------------------------------------------------------------------
__ONMIC__ void correct_errors_3_prime_end(const char*  org_sequence, char*  sequence, const char* quality_score,
        const uint32_t& index_start, char* sequence_modification, uint32_t& trim_5_end, uint32_t& trim_3_end,
        const uint32_t& org_boundary, const uint32_t& read_length, const uint32_t& max_trimmed_bases,
        const uint32_t& min_check_length, uint32_t& num_corrected_errors_local, const unsigned char* bit_vector,
        const uint32_t*  hash_seed, init_args& arg, func_mem& funcm) {
    //  last solid region | non-solid region |
    // --------------------------------------| read
    //               |-----|                   (index_start)-th k-mer
    //--------------------------------------------------
    // list of candidate paths
    vector_candidate_path candidate_path_vector_tmp;
    init_vector_candidate_path(candidate_path_vector_tmp);

    // make an initial k-mer
    //   std::string kmer_initial(sequence.substr(index_start, kmer_length));
    //	char* kmer_initial = (char*) malloc(arg.kmer_length + 2);
    //	memcpy(kmer_initial, sequence + index_start, arg.kmer_length);

    char* kmer_initial = funcm.kmer + funcm.kmer_ind * funcm.kmer_length;
    copy_kmer(funcm, sequence + index_start, arg.kmer_length);

    //	   std::cout << "index start" << index_start << std::endl;
    //	   std::cout << "kmer initial" << kmer_initial << std::endl;

    //	printf("%d \n", index_start);
    //	printf("read length %d \n", read_length);
    //	printf("%s \n", kmer_initial);
    // each alternative neocletide
    bool run_exploration = (true);
    path_query_text(kmer_initial, bit_vector,hash_seed,  arg.kmer_length, arg.num_hash_func, arg.bit_vector_width, 1);
    for (unsigned short int it_alter = G; it_alter <= G; it_alter++) {
        // make a change
        kmer_initial[arg.kmer_length - 1] = NEOCLEOTIDE[it_alter];

        // kmer_initial is solid
//        if (query_text(kmer_initial, bit_vector, &hash_seed[0],arg.kmer_length,arg.num_hash_func,arg.bit_vector_width) == true) {
        if(kmer_initial[funcm.kmer_length + it_alter]){
            // generate a new path
            C_candidate_path candidate_path;
            candidate_path.size = 0;
            candidate_path.sum_qs = 0;
            //			printf(" candidate_path %d \n", candidate_path.size);
            if (sequence[index_start + arg.kmer_length - 1] != NEOCLEOTIDE[it_alter]) {
                pair_uint_char pair_tmp;
                pair_tmp.first  = index_start + arg.kmer_length - 1;
                pair_tmp.second = NEOCLEOTIDE[it_alter];

                push_back(candidate_path,pair_tmp);
            }

            // if this k-mer is the last k-mer in a read
            // running extend_a_kmer_3_prime_end is not needed any more
            if (index_start == (read_length - arg.kmer_length)) {
                push_back(candidate_path_vector_tmp,candidate_path);
            }
            else if (index_start < (read_length - arg.kmer_length)) {
                // trace  this kmer recursively and update candidate_path_vector_tmp
                extend_a_kmer_3_prime_end(
                        kmer_initial,
                        sequence,
                        index_start,
                        candidate_path,
                        candidate_path_vector_tmp,
                        org_boundary,
                        quality_score,
                        bit_vector,
                        hash_seed,
                        run_exploration,
                        read_length,
                        arg,
                        funcm
                        );
                //				printf("this is done \n");
            }
        }
    }
    //	for(int i = 0; i < candidate_path_vector_tmp.size; i ++){
    //		for(int j = 0; j < candidate_path_vector_tmp.c_path[i].size; j++){
    //			printf("%d %c\n", candidate_path_vector_tmp.c_path[i].modified_bases[j].first,
    //					candidate_path_vector_tmp.c_path[i].modified_bases[j].second);
    //		}
    //	}

    // complete exploration was not done because there are too many candidata paths
    // remove all the paths in candidate_path_vector_tmp
    if (run_exploration == false) {
        clear(candidate_path_vector_tmp);
    }

    vector_candidate_path candidate_path_vector_tmp_tmp;
    init_vector_candidate_path(candidate_path_vector_tmp_tmp);

    // check the solidness of the rightmost k-mers of each modified base
    bool all_solid_wo_modification = (false);
    char* sequence_tmp ; //= (char*) malloc(READ_MAX_LEN);
    init_sequence(funcm, &sequence_tmp);
    // each candidate path
    for (C_candidate_path* it_path = candidate_path_vector_tmp.c_path; it_path != candidate_path_vector_tmp.c_path +candidate_path_vector_tmp.size ; it_path++) {

        if ((*it_path).size== 0) {
            all_solid_wo_modification = true;
            break;
        }
        else {
            //         std::string sequence_tmp(sequence);

            memcpy(sequence_tmp, sequence, read_length);
            // index_largest_modified
            uint32_t index_largest_modified = ((*it_path).modified_bases[(*it_path).size - 1].first);

            // number of bases that should be extended
            uint32_t extend_amount;

            // calculate extend_amount
            // no extension is needed
            // sequence.length() = 20, kmer_length = 11, max_extension = 5
            // |0|0|0|1|1|1|1|1|1|1|1|1|1|
            // |7|8|9|0|1|2|3|4|5|6|7|8|9|
            //     |<------------------->| k = 11
            // --------------------------| read
            // ----->|                     index_largest_modified <= 9
            if (index_largest_modified <= read_length - arg.kmer_length) {
                push_back(candidate_path_vector_tmp_tmp, *it_path);
            }
            // extension is needed
            else {
                // applied the modified bases to sequence_tmp
                for (unsigned int it_base = 0; it_base < (*it_path).size; it_base++) {
                    // modify sequence_tmp
                    sequence_tmp[(*it_path).modified_bases[it_base].first] = (*it_path).modified_bases[it_base].second;
                }

                // determine the number of extensions
                // sequence.length() = 20, kmer_length = 11, max_extension = 5
                // |0|0|0|1|1|1|1|1|1|1|1|1|1|0|0|0|0|0|
                // |7|8|9|0|1|2|3|4|5|6|7|8|9|1|2|3|4|5|
                //     |<------------------->|           k = 11
                // --------------------------|           read
                //       |<------->|                     (index_largest_modified > 10) AND (index_largest_modified <= 14)
                //           |<------------------->|     index_largest_modified = 12 -> extend_amout = 3
                if (index_largest_modified <= read_length + arg.max_extension - arg.kmer_length) {
                    extend_amount = arg.kmer_length - (read_length - index_largest_modified);
                }
                // sequence.length() = 20, kmer_length = 11, max_extension = 5
                // |0|0|0|1|1|1|1|1|1|1|1|1|1|0|0|0|0|0|
                // |7|8|9|0|1|2|3|4|5|6|7|8|9|1|2|3|4|5|
                //     |<------------------->|           k = 11
                // --------------------------|           read
                //                 |<------->|           index_largest_modified > 15
                else {
                    extend_amount =arg.max_extension;
                }

                bool extension_success = (false);

                // generate an initial k-mer
                // sequence.length() = 20, kmer_length = 11
                // |0|0|0|1|1|1|1|1|1|1|1|1|1|
                // |7|8|9|0|1|2|3|4|5|6|7|8|9|
                //       |<----------------->| kmer_length - 1 = 10
                // --------------------------| read
                //       |-|                   20 - 11 + 1 = 10
                //            std::string kmer_initial(sequence_tmp.substr(sequence.length() - kmer_length + 1, kmer_length - 1));
                //            kmer_initial = kmer_initial + '0';

                //				char* kmer_initial = (char*) malloc(arg.kmer_length + 2);
                //				memcpy(kmer_initial, sequence_tmp + read_length - arg.kmer_length + 1, arg.kmer_length - 1);
                //				kmer_initial[arg.kmer_length - 1] = '0';
                char* kmer_initial = funcm.kmer + funcm.kmer_ind * funcm.kmer_length;
                copy_kmer(funcm, sequence_tmp + read_length - arg.kmer_length + 1, arg.kmer_length - 1);
                kmer_initial[arg.kmer_length - 1] = 'A';
                path_query_text(kmer_initial, bit_vector,hash_seed,  arg.kmer_length, arg.num_hash_func, arg.bit_vector_width, 1);
                kmer_initial[arg.kmer_length - 1] = '0';
                // each alternative neocletide
                for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
                    // make a change
                    kmer_initial[arg.kmer_length - 1] = NEOCLEOTIDE[it_alter];

                    // kmer_initial is solid
//                    if (query_text(kmer_initial, bit_vector, &hash_seed[0],arg.kmer_length,arg.num_hash_func,arg.bit_vector_width) == true) {
                    if(kmer_initial[funcm.kmer_length + it_alter]){
                        // if extend_amount == 1
                        // running extend_out_right is not needed any more
                        if (extend_amount == 1) {
                            extension_success = true;
                            break;
                        }
                        else {
                            // trace  this kmer recursively and update candidate_path_vector_tmp
                            extend_out_right(
                                    kmer_initial,
                                    1,
                                    extend_amount,
                                    extension_success,
                                    bit_vector,
                                    hash_seed,
                                    arg,
                                    funcm
                                    );
                        }
                    }
                }

                if (extension_success == true) {
                    push_back(candidate_path_vector_tmp_tmp, *it_path);
                }
                //				free(kmer_initial);
                funcm.kmer_ind --;
            }
        }
    }
    //	free(sequence_tmp);
    funcm.sequence_ind --;
    // remain only really modified paths
    clear(candidate_path_vector_tmp);

    // each path
    for (uint32_t it_candidate = 0; it_candidate < candidate_path_vector_tmp_tmp.size; it_candidate++) {
        // each modification
        bool really_modified = (false);
        for (uint32_t it_mod_base = 0; it_mod_base < candidate_path_vector_tmp_tmp.c_path[it_candidate].size; it_mod_base++) {
            //			printf("")
            //			printf("%c \n", org_sequence[candidate_path_vector_tmp_tmp.c_path[it_candidate].modified_bases[it_mod_base].first]);
            //			printf("%c \n", candidate_path_vector_tmp_tmp.c_path[it_candidate].modified_bases[it_mod_base].second);
            if (org_sequence[candidate_path_vector_tmp_tmp.c_path[it_candidate].modified_bases[it_mod_base].first] != candidate_path_vector_tmp_tmp.c_path[it_candidate].modified_bases[it_mod_base].second) {
                really_modified = true;
            }
        }

        if (really_modified) {
            push_back(candidate_path_vector_tmp, candidate_path_vector_tmp_tmp.c_path[it_candidate]);
        }
    }


    //	for(uint32_t i = 0; i < candidate_path_vector_tmp.size; i ++){
    //		for(uint32_t j = 0; j < candidate_path_vector_tmp.c_path[i].size; j ++){
    //			printf("%d %c\n",candidate_path_vector_tmp.c_path[i].modified_bases[j].first, candidate_path_vector_tmp.c_path[i].modified_bases[j].second);
    //		}
    //		printf("\n");
    //	}
    //printf("%d \n", 22222);
    //	candidate_path_vector_tmp_tmp = candidate_path_vector_tmp;
    vector_candidate_path_assignment(&candidate_path_vector_tmp_tmp, &candidate_path_vector_tmp);
    //		for(uint32_t i = 0; i < candidate_path_vector_tmp_tmp.size; i ++){
    //			for(uint32_t j = 0; j < candidate_path_vector_tmp_tmp.c_path[i].size; j ++){
    //				printf("%d %c\n",candidate_path_vector_tmp_tmp.c_path[i].modified_bases[j].first, candidate_path_vector_tmp_tmp.c_path[i].modified_bases[j].second);
    //			}
    //			printf("\n");
    //		}
    // all k-mers are solid without any modification
    // do nothing

    if (all_solid_wo_modification == true) {
    }
    // compare quality scores of candidate paths
    // if the number of paths in candidate_path_vector_tmp_tmp is larger than 1
    else if (candidate_path_vector_tmp_tmp.size > 1) {
        // each path
        C_candidate_path* it_path;
        C_candidate_path* it_path_1st;
        C_candidate_path* it_path_2nd;

        uint32_t qs_1st = (INIT_MIN_QS);
        uint32_t qs_2nd = (INIT_MIN_QS);

        // each candidate path
        for (it_path = candidate_path_vector_tmp_tmp.c_path; it_path != candidate_path_vector_tmp_tmp.c_path + candidate_path_vector_tmp_tmp.size; it_path++) {
            // each modification
            for (uint32_t it_mod = 0; it_mod < (*it_path).size; it_mod++) {
                // add quality scores of modified bases
                (*it_path).sum_qs += ((unsigned short int)quality_score[(*it_path).modified_bases[it_mod].first] - arg.quality_score_offset);
            }

            // compare quality scores of each path
            if ((*it_path).sum_qs <= qs_1st) {
                qs_2nd = qs_1st;
                qs_1st = (*it_path).sum_qs;

                it_path_2nd = it_path_1st;
                it_path_1st = it_path;
            }
            else if ((*it_path).sum_qs <= qs_1st) {
                qs_2nd = (*it_path).sum_qs;

                it_path_2nd = it_path;
            }
        }

        // check whether too many bases are modified in the first path, which implies indel may exist in the original read
        // if too many modifications exist, the read is just trimmed
        //bool keep_going(true);
        bool too_many_corrections = (false);
        //		printf("candidate_path_vector_tmp_tmp size %d \n", candidate_path_vector_tmp_tmp.size);
        //		for(uint32_t i = 0; i < candidate_path_vector_tmp_tmp.size; i ++){
        //			for(uint32_t j = 0; j < candidate_path_vector_tmp_tmp.c_path[i].size; j ++){
        //				printf("%d %c\n",candidate_path_vector_tmp_tmp.c_path[i].modified_bases[j].first, candidate_path_vector_tmp_tmp.c_path[i].modified_bases[j].second);
        //			}
        //			printf("\n");
        //		}
        // at least over MAX_MODIFICATION times of modifications from the first modified position
        if (((read_length - (*it_path_1st).modified_bases[0].first) >= min_check_length) &&
                ((*it_path_1st).size > MAX_MODIFICATION)) {
            // each modified base in a right range
            //std::uint32_t partial_qs((*it_path_1st).sum_qs);

            for (unsigned int it_mod_base = 0; it_mod_base < ((*it_path_1st).size - MAX_MODIFICATION); it_mod_base++) {
                // at least MAX_MODIFICATION times of modifications within min_check_length bases
                if (((*it_path_1st).modified_bases[it_mod_base + MAX_MODIFICATION].first - (*it_path_1st).modified_bases[it_mod_base].first) < min_check_length) {
                    //if (it_mod_base > 0) {
                    //   partial_qs -= ((unsigned short int)quality_score[(*it_path_1st).modified_bases[it_mod_base - 1].first] - quality_score_offset);
                    //}

                    // average quality score of the modified bases is too high
                    //if (1.0 * partial_qs / ((*it_path_1st).modified_bases.size() - it_mod_base) > quality_score_cutoff) {
                    // trim_5_end or trim_3_end
                    if ((read_length - (*it_path_1st).modified_bases[it_mod_base].first) <= max_trimmed_bases) {
                        trim_3_end = read_length - (*it_path_1st).modified_bases[it_mod_base].first;
                    }
                    else if (((*it_path_1st).modified_bases[it_mod_base].first + 1) <= max_trimmed_bases) {
                        trim_5_end = (*it_path_1st).modified_bases[it_mod_base].first + 1;
                    }

                    too_many_corrections = true;
                    break;
                    //}
                }
            }
        }
        // use the 1st path if not too many corrections are made AND if the 1st path has a sufficiently low score
        // if an indel exists using quality scores is not a good way to choose the best path
        if ((too_many_corrections == false) && (qs_1st + MIN_QS_DIFF <= qs_2nd)) {

            // update sequence_modification
            for (unsigned int it_base = 0; it_base < (*it_path_1st).size; it_base++) {
                sequence_modification[(*it_path_1st).modified_bases[it_base].first] = (*it_path_1st).modified_bases[it_base].second;
                sequence[(*it_path_1st).modified_bases[it_base].first] = (*it_path_1st).modified_bases[it_base].second;

                num_corrected_errors_local++;
            }

        }
        // hard to choose one path
        else {

            // A AND B
            vector_pair_uint_char base_vector_intersection;
            base_vector_intersection.size = 0;
            base_vector_intersection.sum_qs = 0;
            // A OR B
            vector_pair_uint_char base_vector_union;
            base_vector_union.size = 0;
            base_vector_union.sum_qs = 0;

            // temporary vectors
            vector_pair_uint_char base_vector_intersection_prev;
            memcpy(&base_vector_intersection_prev, &(candidate_path_vector_tmp_tmp.c_path[0]), sizeof(vector_pair_uint_char));
            vector_pair_uint_char base_vector_union_prev;
            memcpy(&base_vector_union_prev, &(candidate_path_vector_tmp_tmp.c_path[0]), sizeof(vector_pair_uint_char));

            // each candidate path
            for (uint32_t it_p = 1; it_p < candidate_path_vector_tmp_tmp.size; it_p++) {
                clear(base_vector_intersection);
                clear(base_vector_union);

                base_intersection(base_vector_intersection_prev.uint_char, base_vector_intersection_prev.uint_char + base_vector_intersection_prev.size,
                        candidate_path_vector_tmp_tmp.c_path[it_p].modified_bases , candidate_path_vector_tmp_tmp.c_path[it_p].modified_bases + candidate_path_vector_tmp_tmp.c_path[it_p].size, base_vector_intersection);
                base_union       (base_vector_union_prev.uint_char,        base_vector_union_prev.uint_char + base_vector_union_prev.size,
                        candidate_path_vector_tmp_tmp.c_path[it_p].modified_bases, candidate_path_vector_tmp_tmp.c_path[it_p].modified_bases +candidate_path_vector_tmp_tmp.c_path[it_p].size , base_vector_union);

                //				base_vector_intersection_prev = base_vector_intersection;
                //				base_vector_union_prev        = base_vector_union;
                memcpy(&base_vector_intersection_prev, &base_vector_intersection, sizeof(vector_pair_uint_char));
                memcpy(&base_vector_union_prev, &base_vector_union, sizeof(vector_pair_uint_char));
            }

            // A - B
            vector_pair_uint_char base_vector_difference;
            base_vector_difference.size = 0;
            base_vector_difference.sum_qs = 0;
            base_difference(base_vector_union.uint_char, base_vector_union.uint_char + base_vector_union.size,
                    base_vector_intersection.uint_char,base_vector_intersection.uint_char + base_vector_intersection.size, base_vector_difference);

            // find trimmed region
            // correcting the 5'-end and 3'-end was not done
            // therefore the total number of trimmed bases is 0 yet
            if (base_vector_difference.size > 0) {
                uint32_t vector_index_leftmost = (0);
                uint32_t vector_index_rightmost = (base_vector_difference.size - 1);

                bool keep_going = (true);

                while (keep_going) {
                    // # of trimmed bases at the 5' end: base_vector_difference[vector_index_leftmost].first + 1
                    // # of trimmed bases at the 3' end: read_length - base_vector_difference[vector_index_rightmost].first
                    // the 5'-end is smaller
                    if ((base_vector_difference.uint_char[vector_index_leftmost].first + 1) < (read_length - base_vector_difference.uint_char[vector_index_rightmost].first)) {
                        // check the total number of trimmed bases
                        if ((base_vector_difference.uint_char[vector_index_leftmost].first + 1 + trim_3_end) <= max_trimmed_bases) {
                            if ((base_vector_difference.uint_char[vector_index_leftmost].first + 1) > trim_5_end) {
                                trim_5_end = base_vector_difference.uint_char[vector_index_leftmost].first + 1;
                            }

                            // two points are met
                            if (vector_index_leftmost == vector_index_rightmost) {
                                keep_going = false;
                            }
                            else {
                                vector_index_leftmost++;
                            }
                        }
                        // no need for more check
                        else {
                            keep_going = false;
                        }
                    }
                    // the 3'-end is smaller
                    else {
                        // check the total number of trimmed bases
                        if ((read_length - base_vector_difference.uint_char[vector_index_rightmost].first) <= max_trimmed_bases) {
                            if ((read_length - base_vector_difference.uint_char[vector_index_rightmost].first) > trim_3_end) {
                                trim_3_end = read_length - base_vector_difference.uint_char[vector_index_rightmost].first;
                            }

                            // two points are met
                            if (vector_index_leftmost == vector_index_rightmost) {
                                keep_going = false;
                            }
                            else {
                                vector_index_rightmost--;
                            }
                        }
                        // no need for more check
                        else {
                            keep_going = false;
                        }
                    }
                }
            }

            // find consensus modifications
            for (uint32_t it_inter = 0; it_inter < base_vector_intersection.size; it_inter++) {
                // check whether the base is not in the trimmed regions
                if ((base_vector_intersection.uint_char[it_inter].first <  (read_length - trim_3_end)) &&
                        (base_vector_intersection.uint_char[it_inter].first >= trim_5_end)) {
                    // filter out the bases that are equal to the original ones
                    if (sequence[base_vector_intersection.uint_char[it_inter].first] != base_vector_intersection.uint_char[it_inter].second) {
                        //						printf("%c \n", sequence_modification[base_vector_intersection.uint_char[it_inter].first]);
                        //						printf("%c \n", base_vector_intersection.uint_char[it_inter].second);
                        sequence_modification[base_vector_intersection.uint_char[it_inter].first] = base_vector_intersection.uint_char[it_inter].second;
                        sequence[base_vector_intersection.uint_char[it_inter].first] = base_vector_intersection.uint_char[it_inter].second;

                        num_corrected_errors_local++;
                    }
                }
            }

        }
    }
    // only one path
    else if (candidate_path_vector_tmp_tmp.size == 1) {

        // check whether too many bases are modified in the first path, which implies indel may exist in the original read
        // if too many modifications exist, the read is just trimmed
        //bool keep_going(true);
        bool too_many_corrections = (false);

        // at least over MAX_MODIFICATION times of modifications from the first modified position
        if (((read_length - (candidate_path_vector_tmp_tmp.c_path[0]).modified_bases[0].first) >= min_check_length) &&
                ((candidate_path_vector_tmp_tmp.c_path[0]).size > MAX_MODIFICATION)) {
            // calculate sum_qs
            //for (unsigned int it_mod_base = 0; it_mod_base < candidate_path_vector_tmp_tmp[0].modified_bases.size(); it_mod_base++) {
            //   candidate_path_vector_tmp_tmp[0].sum_qs += ((unsigned short int)quality_score[candidate_path_vector_tmp_tmp[0].modified_bases[it_mod_base].first] - quality_score_offset);
            //}

            // each modified base in a right range
            //std::uint32_t partial_qs((candidate_path_vector_tmp_tmp[0]).sum_qs);

            for (unsigned int it_mod_base = 0; it_mod_base < ((candidate_path_vector_tmp_tmp.c_path[0]).size - MAX_MODIFICATION); it_mod_base++) {
                // at least MAX_MODIFICATION times of modifications within min_check_length bases
                if (((candidate_path_vector_tmp_tmp.c_path[0]).modified_bases[it_mod_base + MAX_MODIFICATION].first - (candidate_path_vector_tmp_tmp.c_path[0]).modified_bases[it_mod_base].first) < min_check_length) {
                    //if (it_mod_base > 0) {
                    //   partial_qs -= ((unsigned short int)quality_score[(candidate_path_vector_tmp_tmp[0]).modified_bases[it_mod_base - 1].first] - quality_score_offset);
                    //}

                    // average quality score of the modified bases is too high
                    //if (1.0 * partial_qs / ((candidate_path_vector_tmp_tmp[0]).modified_bases.size() - it_mod_base) > quality_score_cutoff) {
                    // trim_5_end or trim_3_end
                    if ((read_length - (candidate_path_vector_tmp_tmp.c_path[0]).modified_bases[it_mod_base].first) <= max_trimmed_bases) {
                        trim_3_end = read_length - (candidate_path_vector_tmp_tmp.c_path[0]).modified_bases[it_mod_base].first;

                        // update sequence_modification for the non-trimmed corrections
                        if (it_mod_base > 0) {
                            for (unsigned int it_base = 0; it_base < (it_mod_base - 1); it_base++) {
                                sequence_modification[(candidate_path_vector_tmp_tmp.c_path[0]).modified_bases[it_base].first] = (candidate_path_vector_tmp_tmp.c_path[0]).modified_bases[it_base].second;
                                sequence[(candidate_path_vector_tmp_tmp.c_path[0]).modified_bases[it_base].first] = (candidate_path_vector_tmp_tmp.c_path[0]).modified_bases[it_base].second;

                                num_corrected_errors_local++;
                            }
                        }
                    }
                    else if (((candidate_path_vector_tmp_tmp.c_path[0]).modified_bases[it_mod_base].first + 1) <= max_trimmed_bases) {
                        trim_5_end = (candidate_path_vector_tmp_tmp.c_path[0]).modified_bases[it_mod_base].first + 1;

                        // update sequence_modification for the non-trimmed corrections
                        if (it_mod_base < (candidate_path_vector_tmp_tmp.c_path[0].size - 1)) {
                            for (unsigned int it_base = (it_mod_base + 1); it_base < candidate_path_vector_tmp_tmp.c_path[0].size; it_base++) {
                                sequence_modification[(candidate_path_vector_tmp_tmp.c_path[0]).modified_bases[it_base].first] = (candidate_path_vector_tmp_tmp.c_path[0]).modified_bases[it_base].second;
                                sequence[(candidate_path_vector_tmp_tmp.c_path[0]).modified_bases[it_base].first] = (candidate_path_vector_tmp_tmp.c_path[0]).modified_bases[it_base].second;

                                num_corrected_errors_local++;
                            }
                        }
                    }

                    too_many_corrections = true;
                    break;
                    //}
                }
            }
        }

        // not too many modified bases
        if (too_many_corrections == false) {
            // update sequence_modification
            for (unsigned int it_base = 0; it_base < candidate_path_vector_tmp_tmp.c_path[0].size; it_base++) {
                sequence_modification[(candidate_path_vector_tmp_tmp.c_path[0]).modified_bases[it_base].first] = (candidate_path_vector_tmp_tmp.c_path[0]).modified_bases[it_base].second;
                sequence[(candidate_path_vector_tmp_tmp.c_path[0]).modified_bases[it_base].first] = (candidate_path_vector_tmp_tmp.c_path[0]).modified_bases[it_base].second;

                num_corrected_errors_local++;
            }
        }
    }

    //	free(kmer_initial);
    funcm.kmer_ind --;
    free(candidate_path_vector_tmp.c_path);
    free(candidate_path_vector_tmp_tmp.c_path);

}


//__ONMIC__
__ONMIC__ void mic_correct_errors(char* read_content, char* read_quality, char* sequence_modification,
        uint32_t* hash_seed, uint8_t* bit_vector, init_args& arg,
        uint32_t& num_corrected_errors_local,uint32_t trim_5_end, uint32_t trim_3_end,
        uint32_t max_trimmed_bases,uint32_t min_check_length, uint32_t read_length, int ssss, func_mem& funcm) {
    vector_pair_uint_uint solid_regions;
    vector_pair_uint_uint solid_regions_org;
    solid_regions.size = 0;
    solid_regions_org.size = 0;
    //--------------------------------------------------
    // STEP 0-0: find solid k-mers in this read
    //--------------------------------------------------
    // variables
    uint32_t num_kmers =read_length- arg.kmer_length + 1;

    //	std::vector< std::pair<std::uint32_t, std::uint32_t> > solid_regions;
    //	std::vector< std::pair<std::uint32_t, std::uint32_t> > solid_regions_org;

    bool is_solid_kmer_prev = (false);
    // find solid regions
    pair_uint_uint new_solid_region;
    regions_query_text(read_content, bit_vector, hash_seed, arg.kmer_length, arg.num_hash_func,
            arg.bit_vector_width, funcm.kmer, num_kmers, funcm.kmer_index);
//    if(ssss == 0){
//    for(uint32_t i = 0; i < num_kmers; i ++){
//         printf("%d ", funcm.kmer_index[i]);
//    }
//    printf("\n");
//    }

    for (uint32_t it_kmer = 0; it_kmer < num_kmers; it_kmer++) {
        //	      std::string current_kmer(sequence.substr(it_kmer, kmer_length));
        char* current_kmer = read_content + it_kmer;
        // k-mer is solid
//       if (query_text(current_kmer, bit_vector, hash_seed, arg.kmer_length, arg.num_hash_func, arg.bit_vector_width)  == true) {
        if(funcm.kmer_index[it_kmer] >= 600){
        // start point of a solid region
            if (is_solid_kmer_prev == false) {
                new_solid_region.first = it_kmer;
                //	        	 std::cout << new_solid_region.first << "  1  " << new_solid_region.second <<std::endl;
                is_solid_kmer_prev = true;
            }
        }
        else {
            // end point of a solid region
            if (is_solid_kmer_prev == true) {
                new_solid_region.second = it_kmer - 1;
                //	        	 std::cout << new_solid_region.first << "  2  " << new_solid_region.second <<std::endl;

                if (new_solid_region.second < new_solid_region.first) {
                    //					std::cout << std::endl << "ERROR: The second index is smaller than the first" << std::endl << std::endl;
                    printf("%s\n", "ERROR: The second index is smaller than the first");
                }
                push_back(solid_regions,new_solid_region);
                push_back(solid_regions_org,new_solid_region);

                is_solid_kmer_prev = false;
            }
        }
    }

    // last solid region
    if (is_solid_kmer_prev == true) {
        new_solid_region.second = num_kmers - 1;
        //		solid_regions.push_back(new_solid_region);
        //		solid_regions_org.push_back(new_solid_region);
        push_back(solid_regions,new_solid_region);
        push_back(solid_regions_org,new_solid_region);
    }


    vector_pair_uint_uint solid_regions_tmp;
    vector_pair_uint_uint solid_regions_org_tmp;
    //--------------------------------------------------
    // STEP 0-1: adjust solid regioins using quality scores
    //--------------------------------------------------
    // solid_regions_org: indexes that are not modified
    // when the indices are reached, all kinds of modifications (A/C/G/T) should be made and checked
    // at least one solid region
    if (solid_regions.size > 0) {
        // exceptional case: only one solid island that covers the entire read
        if ((solid_regions.size != 1) || (solid_regions.uint_uint[0].first != 0) || (solid_regions.uint_uint[0].second != (read_length - arg.kmer_length))) {
            // at least two solid k-mer islands
            if (solid_regions.size > 1) {
                // check the distance between every two solid k-mer islands
                bool flag_short_distance(false);

                for (uint32_t it_sr = 0; it_sr < solid_regions.size - 1; it_sr++) {
                    if ((solid_regions.uint_uint[it_sr + 1].first - solid_regions.uint_uint[it_sr].second) < arg.kmer_length) {
                        flag_short_distance = true;
                    }
                }

                if (flag_short_distance == true) {
                    solid_regions_tmp.size = 0;
                    solid_regions_org_tmp.size = 0;
                    // each solid island
                    for (uint32_t it_sr = 0; it_sr < solid_regions.size; it_sr++) {
                        // each base in the solid island (0-base)
                        uint32_t num_low_quality_base = (0);
                        // set an initial value to avoid a compilation warning
                        uint32_t index_prev_low_quality_base = (0);
                        for (unsigned int it_base = solid_regions.uint_uint[it_sr].first; it_base < solid_regions.uint_uint[it_sr].second + arg.kmer_length; it_base++) {
                            // a current base has a low quality score
                            if (((unsigned short int)read_quality[it_base] - arg.quality_score_offset) < arg.quality_score_cutoff) {
                                num_low_quality_base++;

                                // first low quality base
                                if (num_low_quality_base == 1) {
                                    // the low quality base is not in the first k-mer of the solid island
                                    if (it_base >= (solid_regions.uint_uint[it_sr].first + arg.kmer_length)) {
                                        // add the left most high quality region to a temporary vector
                                        //										std::pair<std::uint32_t, std::uint32_t> new_solid_region;
                                        pair_uint_uint new_solid_region;

                                        new_solid_region.first  = solid_regions.uint_uint[it_sr].first;
                                        new_solid_region.second = it_base - arg.kmer_length;
                                        //										solid_regions_tmp.push_back(new_solid_region);
                                        //
                                        //										solid_regions_org_tmp.push_back(solid_regions_org[it_sr]);
                                        push_back(solid_regions_tmp,new_solid_region);
                                        push_back(solid_regions_org_tmp,solid_regions_org.uint_uint[it_sr]);
                                    }
                                }
                                // not first low quality base
                                else {
                                    if ((it_base - index_prev_low_quality_base) > arg.kmer_length) {
                                        pair_uint_uint new_solid_region;

                                        new_solid_region.first = index_prev_low_quality_base + 1;
                                        new_solid_region.second = it_base - arg.kmer_length;
                                        //										solid_regions_tmp.push_back(new_solid_region);
                                        //
                                        //										solid_regions_org_tmp.push_back(solid_regions_org[it_sr]);

                                        push_back(solid_regions_tmp,new_solid_region);
                                        push_back(solid_regions_org_tmp,solid_regions_org.uint_uint[it_sr]);
                                    }
                                }

                                index_prev_low_quality_base = it_base;
                            }
                        }

                        // process the bases to the right of the rightmost low quality base
                        if (num_low_quality_base > 0) {
                            if (solid_regions.uint_uint[it_sr].second >= (index_prev_low_quality_base + arg.kmer_length)) {
                                //								std::pair<std::uint32_t, std::uint32_t> new_solid_region;
                                pair_uint_uint new_solid_region;

                                new_solid_region.first  = index_prev_low_quality_base + 1;
                                new_solid_region.second = solid_regions.uint_uint[it_sr].second;
                                //								solid_regions_tmp.push_back(new_solid_region);
                                //
                                //								solid_regions_org_tmp.push_back(solid_regions_org[it_sr]);
                                push_back(solid_regions_tmp,new_solid_region);
                                push_back(solid_regions_org_tmp,solid_regions_org.uint_uint[it_sr]);
                            }
                        }
                        // no low quality base
                        // add the current solid island
                        else {
                            //							solid_regions_tmp.push_back(solid_regions[it_sr]);
                            //
                            //							solid_regions_org_tmp.push_back(solid_regions_org[it_sr]);
                            push_back(solid_regions_tmp,solid_regions.uint_uint[it_sr]);
                            push_back(solid_regions_org_tmp,solid_regions_org.uint_uint[it_sr]);
                        }
                    }

                    //					solid_regions     = solid_regions_tmp;
                    //					solid_regions_org = solid_regions_org_tmp;
                    memcpy(&solid_regions,&solid_regions_tmp,sizeof(vector_pair_uint_uint));
                    memcpy(&solid_regions_org,&solid_regions_org_tmp,sizeof(vector_pair_uint_uint));
                }
            }
            // only one solid k-mer island
            else if (solid_regions.size == 1) {
                //				std::vector< std::pair<std::uint32_t, std::uint32_t> > solid_regions_tmp;
                //				std::vector< std::pair<std::uint32_t, std::uint32_t> > solid_regions_org_tmp;



                solid_regions_tmp.size = 0;
                solid_regions_org_tmp.size = 0;

                uint32_t num_low_quality_base = (0);
                uint32_t prev_low_quality_index = (0);

                // each base in the solid island (0-base)
                for (unsigned int it_base = solid_regions.uint_uint[0].first; it_base < solid_regions.uint_uint[0].second + arg.kmer_length; it_base++) {
                    // a current base has a low quality score
                    if (((unsigned short int)read_quality[it_base] - arg.quality_score_offset) < arg.quality_score_cutoff) {
                        num_low_quality_base++;

                        // first low quality base
                        if (num_low_quality_base == 1) {
                            if ((it_base - solid_regions.uint_uint[0].first) >= (arg.kmer_length + MIN_SOLID_LENGTH - 1)) {
                                //								std::pair<std::uint32_t, std::uint32_t> new_solid_region;
                                pair_uint_uint new_solid_region;
                                new_solid_region.first = solid_regions.uint_uint[0].first;
                                new_solid_region.second = it_base - arg.kmer_length;
                                //								solid_regions_tmp.push_back(new_solid_region);
                                //
                                //								solid_regions_org_tmp.push_back(solid_regions_org[0]);
                                push_back(solid_regions_tmp, new_solid_region);

                                push_back(solid_regions_org_tmp,solid_regions_org.uint_uint[0]);
                            }

                            prev_low_quality_index = it_base;
                        }
                        // not first low quality base
                        else {
                            if ((it_base - prev_low_quality_index) >= (arg.kmer_length + MIN_SOLID_LENGTH)) {
                                //								std::pair<std::uint32_t, std::uint32_t> new_solid_region;
                                pair_uint_uint new_solid_region;

                                new_solid_region.first = prev_low_quality_index + 1;
                                new_solid_region.second = it_base - arg.kmer_length;

                                //								solid_regions_tmp.push_back(new_solid_region);
                                //
                                //								solid_regions_org_tmp.push_back(solid_regions_org[0]);

                                push_back(solid_regions_tmp, new_solid_region);

                                push_back(solid_regions_org_tmp,solid_regions_org.uint_uint[0]);

                            }

                            prev_low_quality_index = it_base;
                        }
                    }
                }

                // the above is done only when this procedure does not remove the only solid island
                if (solid_regions_tmp.size > 0) {
                    //					solid_regions     = solid_regions_tmp;
                    //					solid_regions_org = solid_regions_org_tmp;
                    memcpy(&solid_regions,&solid_regions_tmp,sizeof(vector_pair_uint_uint));
                    memcpy(&solid_regions_org,&solid_regions_org_tmp,sizeof(vector_pair_uint_uint));

                }
            }
        }
    }

    //--------------------------------------------------
    // STEP 0-2: remove short solid regions
    //--------------------------------------------------
    if (solid_regions.size > 0) {
        //		std::vector< std::pair<std::uint32_t, std::uint32_t> > solid_regions_tmp;
        //		std::vector< std::pair<std::uint32_t, std::uint32_t> > solid_regions_org_tmp;

        solid_regions_tmp.size = 0;
        solid_regions_org_tmp.size = 0;

        for (uint32_t it_region = 0; it_region < solid_regions.size; it_region++) {
            if ((solid_regions.uint_uint[it_region].second - solid_regions.uint_uint[it_region].first + 1) >= MIN_SOLID_LENGTH) {
                //				solid_regions_tmp.push_back(solid_regions[it_region]);
                //				solid_regions_org_tmp.push_back(solid_regions_org[it_region]);
                push_back(solid_regions_tmp, solid_regions.uint_uint[it_region]);
                push_back(solid_regions_org_tmp,solid_regions_org.uint_uint[it_region]);

            }
        }

        //		solid_regions     = solid_regions_tmp;
        //		solid_regions_org = solid_regions_org_tmp;
        memcpy(&solid_regions,&solid_regions_tmp,sizeof(vector_pair_uint_uint));
        memcpy(&solid_regions_org,&solid_regions_org_tmp,sizeof(vector_pair_uint_uint));
    }
    //--------------------------------------------------
    // STEP 0-3: remove short non-solid regions
    //--------------------------------------------------
    if (solid_regions.size > 0) {
        //		std::vector< std::pair<std::uint32_t, std::uint32_t> > solid_regions_tmp;
        //		std::vector< std::pair<std::uint32_t, std::uint32_t> > solid_regions_org_tmp;

        solid_regions_tmp.size = 0;
        solid_regions_org_tmp.size = 0;
        //		solid_regions_tmp.push_back(solid_regions[0]);
        //		solid_regions_org_tmp.push_back(solid_regions_org[0]);

        push_back(solid_regions_tmp, solid_regions.uint_uint[0]);

        push_back(solid_regions_org_tmp,solid_regions_org.uint_uint[0]);


        if (solid_regions.size > 1) {
            for (uint32_t it_region = 1; it_region < solid_regions.size; it_region++) {
                if ((solid_regions.uint_uint[it_region].first - solid_regions.uint_uint[it_region - 1].second - 1) < MIN_NON_SOLID_LENGTH) {
                    solid_regions_tmp.uint_uint[solid_regions_tmp.size - 1].second = solid_regions.uint_uint[it_region].second;
                }
                else {
                    //					solid_regions_tmp.push_back(solid_regions[it_region]);
                    //					solid_regions_org_tmp.push_back(solid_regions_org[it_region]);
                    push_back(solid_regions_tmp, solid_regions.uint_uint[it_region]);
                    push_back(solid_regions_org_tmp,solid_regions_org.uint_uint[it_region]);

                }
            }
        }
        //		solid_regions     = solid_regions_tmp;
        //		solid_regions_org = solid_regions_org_tmp;

        memcpy(&solid_regions,&solid_regions_tmp,sizeof(vector_pair_uint_uint));
        memcpy(&solid_regions_org,&solid_regions_org_tmp,sizeof(vector_pair_uint_uint));

    }
    //--------------------------------------------------
    // STEP 0-4: reduce the size of solid regions
    //--------------------------------------------------
    if (solid_regions.size > 1) {
        for (uint32_t it_region = 1; it_region < solid_regions.size; it_region++) {
            // (length of a non-solid region < kmer_length) && (length of a non-solid region >= kmer_length - FP_SUSPECT_LENGTH(default: 1))
            if (((solid_regions.uint_uint[it_region].first - solid_regions.uint_uint[it_region - 1].second - 1) < arg.kmer_length) &&
                    ((solid_regions.uint_uint[it_region].first - solid_regions.uint_uint[it_region - 1].second - 1) >= arg.kmer_length - FP_SUSPECT_LENGTH)) {
                // length of the right solid region > FP_SUSPECT_LENGTH(default: 1)
                if ((solid_regions.uint_uint[it_region].second - solid_regions.uint_uint[it_region].first + 1) > FP_SUSPECT_LENGTH) {
                    solid_regions.uint_uint[it_region].first += FP_SUSPECT_LENGTH;
                }

                // length of the left solid region > FP_SUSPECT_LENGTH(default: 1)
                if ((solid_regions.uint_uint[it_region - 1].second - solid_regions.uint_uint[it_region - 1].first + 1) > FP_SUSPECT_LENGTH) {
                    solid_regions.uint_uint[it_region - 1].second -= FP_SUSPECT_LENGTH;
                }
            }
        }
    }
    //--------------------------------------------------
    // STEP 0-5: remove a solid region that makes a non-solid reiong shorter than k
    //--------------------------------------------------
    if (solid_regions.size == 2) {
        // the first solid region starts from the first k-mer
        if (solid_regions.uint_uint[0].first == 0) {
            // the distance between two regions is shorter than k
            if ((solid_regions.uint_uint[1].first - solid_regions.uint_uint[0].second) < (arg.kmer_length + 1)) {
                // remove the second solid region
                //				solid_regions.erase(solid_regions.begin() + 1);
                uint32_t pos = 1;
                erase(solid_regions, pos);

            }
        }
        // the second solid region ends in the last k-mer
        else if (solid_regions.uint_uint[1].second == (read_length - arg.kmer_length)) {
            // the distance between two regions is shorter than k
            if ((solid_regions.uint_uint[1].first - solid_regions.uint_uint[0].second) < (arg.kmer_length + 1)) {
                // the length of the second solid region is >= 10% of the sequence length
                if ((solid_regions.uint_uint[1].second - solid_regions.uint_uint[1].first + 1) >= (read_length * 0.1)) {
                    // the length of the first solid region is < 10% of the sequence length
                    if ((solid_regions.uint_uint[0].second - solid_regions.uint_uint[0].first + 1) < (read_length * 0.1)) {
                        // remove the second solid region
                        //						solid_regions.erase(solid_regions.begin());
                        uint32_t pos = 0;
                        erase(solid_regions, pos);
                    }
                }
            }
        }
    }
    //--------------------------------------------------
    // STEP 0-6: check the quality scores of right side of each solid k-mer region
    //--------------------------------------------------
    // at least one solid region
    if (solid_regions.size > 0) {
        // 1 - (n - 1) solid region
        for (uint32_t it_sr = 0; it_sr < (solid_regions.size - 1); it_sr++) {
            // sufficient solid regions length
            if ((solid_regions.uint_uint[it_sr].second - solid_regions.uint_uint[it_sr].first) > SOLID_REGION_ADJUST_RANGE) {
                for (uint32_t it_adjust = solid_regions.uint_uint[it_sr].second; it_adjust > (solid_regions.uint_uint[it_sr].second - SOLID_REGION_ADJUST_RANGE); it_adjust--) {
                    // low quality score
                    if ((((unsigned short int)read_quality[it_adjust + arg.kmer_length - 1] - arg.quality_score_offset) < arg.quality_score_cutoff) ||
                            (((unsigned short int)read_quality[it_adjust] - arg.quality_score_offset) < arg.quality_score_cutoff)
                       ) {
                        solid_regions.uint_uint[it_sr].second = it_adjust - 1;
                        break;
                    }
                }
            }
        }

        // last solid region
        uint32_t index_solid_region = (solid_regions.size - 1);

        // non-solid k-mers exist at the 3-prime end
        if (solid_regions.uint_uint[index_solid_region].second < (read_length - arg.kmer_length)) {
            // sufficient solid regions length
            if ((solid_regions.uint_uint[index_solid_region].second - solid_regions.uint_uint[index_solid_region].first) > SOLID_REGION_ADJUST_RANGE) {
                for (uint32_t it_adjust = solid_regions.uint_uint[index_solid_region].second; it_adjust > (solid_regions.uint_uint[index_solid_region].second - SOLID_REGION_ADJUST_RANGE); it_adjust--) {
                    // low quality score
                    if (((unsigned short int)read_quality[it_adjust + arg.kmer_length - 1] - arg.quality_score_offset) < arg.quality_score_cutoff) {
                        solid_regions.uint_uint[index_solid_region].second = it_adjust - 1;
                        break;
                    }
                }
            }
        }

        // non-solid k-mers exist at the 5-prime end
        if (solid_regions.uint_uint[0].first > 0) {
            // sufficient solid regions length
            if ((solid_regions.uint_uint[0].second - solid_regions.uint_uint[0].first) > SOLID_REGION_ADJUST_RANGE) {
                for (uint32_t it_adjust = solid_regions.uint_uint[0].first; it_adjust < (solid_regions.uint_uint[0].first + SOLID_REGION_ADJUST_RANGE); it_adjust++) {
                    // low quality score
                    if (((unsigned short int)read_quality[it_adjust + arg.kmer_length - 1] - arg.quality_score_offset) < arg.quality_score_cutoff) {
                        solid_regions.uint_uint[0].first = it_adjust + 1;
                        break;
                    }
                }
            }
        }
    }
    //--------------------------------------------------
    // STEP 0-7: check whether a non-solid region < k still exists
    //--------------------------------------------------
    bool short_non_solid_region(false);
    if (solid_regions.size > 1) {
        for (uint32_t it_sr = 1; it_sr < (solid_regions.size - 1); it_sr++) {
            if ((solid_regions.uint_uint[it_sr].first - solid_regions.uint_uint[it_sr - 1].second) <= arg.kmer_length) {
                short_non_solid_region = true;
                break;
            }
        }
    }

    //			   if(short_non_solid_region == false && solid_regions.size > 1){
    //			            std::cout << solid_regions.size <<  "HHHHHH" << ssss << std::endl;
    //			            for(int i = 0; i < solid_regions.size; i++ ){
    //			            	std::cout << solid_regions.uint_uint[i].first << "    " << solid_regions.uint_uint[i].second <<std::endl;
    //			            }
    //			   }
    //			   printf("solid region size %d \n", solid_regions.size);
    //			   for(int i = 0; i < solid_regions.size;  i++){
    //				   printf("solid regions %d %d \n", solid_regions.uint_uint[i].first, solid_regions.uint_uint[i].second);
    //				   printf("solid regions org %d %d", solid_regions_org.uint_uint[i].first, solid_regions_org.uint_uint[i].second);
    ////				   std::cout << solid_regions.uint_uint[i].first <<solid_regions.uint_uint[i].second   << std::endl;
    ////				   std::cout << solid_regions_org.uint_uint[i].first <<solid_regions_org.uint_uint[i].second   << std::endl;
    //		   }


    //--------------------------------------------------
    // correct errors
    //--------------------------------------------------
    //	char*  sequence_modified = (char*)malloc(2*READ_MAX_LEN*sizeof(char));
    //	memcpy(sequence_modified, read_content, READ_MAX_LEN);
    char*  sequence_modified;
    init_sequence(funcm, &sequence_modified);
    memcpy(sequence_modified, read_content, READ_MAX_LEN);

    // solid regions exist and none of them is too short
    if ((solid_regions.size > 0) && (short_non_solid_region == false)) {
        //--------------------------------------------------
        // STEP 1-1: Correct errors between solid regions
        //--------------------------------------------------
        if (solid_regions.size > 1) {
            // for each solid region
            for (uint32_t it_region = 1; it_region < solid_regions.size; it_region++) {
                if ((((solid_regions.uint_uint[it_region].first - 1) - (solid_regions.uint_uint[it_region - 1].second + 1)) + 1) >= arg.kmer_length) {
                    // the bases that may be modified: from (solid_regions[it_region - 1].second + kmer_length) to (solid_regions[it_region].first - 1)
                    // they should not be in trimmed regions
                    if (((solid_regions.uint_uint[it_region - 1].second + arg.kmer_length) < (read_length - trim_3_end)) && ((solid_regions.uint_uint[it_region].first) > trim_5_end)) {
                        //		            	   std::cout << solid_regions[it_region - 1].first << "         " << solid_regions[it_region - 1].second + 1 << std::endl;
                        //		            	   std::cout <<  solid_regions[it_region].first - 1 << "         "<<solid_regions[it_region].second << std::endl;
                        //		            	   std::cout <<  solid_regions_org[it_region - 1].second - 1 << "         "<<solid_regions_org[it_region].first - 1 << std::endl;
                        correct_errors_between_solid_regions(
                                read_content,
                                sequence_modified,
                                read_quality,
                                solid_regions.uint_uint[it_region - 1].first,
                                solid_regions.uint_uint[it_region - 1].second + 1,
                                solid_regions.uint_uint[it_region].first - 1,
                                solid_regions.uint_uint[it_region].second,
                                solid_regions_org.uint_uint[it_region - 1].second - 1,
                                solid_regions_org.uint_uint[it_region].first - 1,
                                sequence_modification,
                                trim_5_end,
                                trim_3_end,
                                solid_regions.size,
                                read_length,
                                max_trimmed_bases,
                                min_check_length,
                                num_corrected_errors_local,
                                bit_vector,
                                hash_seed,
                                arg,
                                funcm
                                    );


                    }
                }
            }
        }

        //--------------------------------------------------
        // STEP 1-2: Correct errors in the 5' end
        //--------------------------------------------------
        // number of solid regions is >= 1
        if (solid_regions.size >= 1) {
            // the first solid region does not start from the 0-th k-mer in a read
            if (solid_regions.uint_uint[0].first > 0) {
                // the bases that may be modified: from 0 to (solid_regions[0].first - 1)
                // they should not be in trimmed regions
                if (solid_regions.uint_uint[0].first > trim_5_end) {
                    //					printf("read content %s \n", read_content);
                    correct_errors_5_prime_end(
                            read_content,
                            sequence_modified,
                            read_quality,
                            solid_regions.uint_uint[0].first - 1,
                            sequence_modification,
                            trim_5_end,
                            trim_3_end,
                            solid_regions_org.uint_uint[0].first - 1,
                            read_length,
                            max_trimmed_bases,
                            min_check_length,
                            num_corrected_errors_local,
                            bit_vector,
                            hash_seed,
                            arg,
                            funcm
                            );
                }
            }
        }


        //--------------------------------------------------
        // STEP 1-3: Correct errors in the 3' end
        //--------------------------------------------------
        // number of solid regions is >= 1
        if (solid_regions.size >= 1) {
            // the last solid region does not end in the last k-mer in a read
            if (solid_regions.uint_uint[solid_regions.size - 1].second < (read_length - arg.kmer_length)) {
                // the bases that may be modified: from (solid_regions[solid_regions.size() - 1].second + kmer_length) to (read_length - 1)
                // they should not be in trimmed regions
                if ((solid_regions.uint_uint[solid_regions.size - 1].second + arg.kmer_length) < (read_length - trim_3_end)) {
                    correct_errors_3_prime_end(
                            read_content,
                            sequence_modified,
                            read_quality,
                            solid_regions.uint_uint[solid_regions.size - 1].second + 1,
                            sequence_modification,
                            trim_5_end,
                            trim_3_end,
                            solid_regions_org.uint_uint[solid_regions.size - 1].second + 1,
                            read_length,
                            max_trimmed_bases,
                            min_check_length,
                            num_corrected_errors_local,
                            bit_vector,
                            hash_seed,
                            arg,
                            funcm
                            );
                }
            }
        }

        //
        // check whether any modification was made
        //
        if ((num_corrected_errors_local == 0) && (trim_5_end == 0) && (trim_3_end == 0)) {
            trim_5_end = solid_regions.uint_uint[0].first;

            // trim as many bases as possible
            for (uint32_t it_solid_short = 0; it_solid_short < solid_regions.size; it_solid_short++) {
                if ((trim_5_end + (read_length - arg.kmer_length - solid_regions.uint_uint[it_solid_short].second)) <= max_trimmed_bases) {
                    trim_3_end = read_length - arg.kmer_length - solid_regions.uint_uint[it_solid_short].second;
                    break;
                }
            }
        }

    }

    //--------------------------------------------------
    // no solid region or short weak regions
    //--------------------------------------------------
    else {

        //--------------------------------------------------
        // STEP 2-1: Correct errors in the first k-mer
        //--------------------------------------------------
        // find potentially wrong bases
        vector_candidate_path candidate_path_vector_tmp;
        init_vector_candidate_path(candidate_path_vector_tmp);
        correct_errors_first_kmer(
                sequence_modified,
                read_quality,
                sequence_modification,
                candidate_path_vector_tmp,
                bit_vector,
                hash_seed,
                arg,
                funcm
                );

        // candidiate_path_vector_tmp: differently modified versions of the first k-mer

        // filter some candidates by extending the first k-mer to the left
        vector_candidate_path candidate_path_vector_tmp_tmp;
        init_vector_candidate_path(candidate_path_vector_tmp_tmp);

        if (candidate_path_vector_tmp.size > 0) {
            // each path
            for (uint32_t it_candidates = 0; it_candidates < candidate_path_vector_tmp.size; it_candidates++) {
                // no modified base
                // if the original first k-mer is added to the candidate vector
                if (candidate_path_vector_tmp.c_path[it_candidates].size == 0) {
                    push_back(candidate_path_vector_tmp_tmp, candidate_path_vector_tmp.c_path[it_candidates]);
                }
                // there exist modified bases
                else {
                    bool really_modified(false);
                    for (unsigned int it_mod_base = 0; it_mod_base < candidate_path_vector_tmp.c_path[it_candidates].size; it_mod_base++) {
                        if (read_content[candidate_path_vector_tmp.c_path[it_candidates].modified_bases[it_mod_base].first] != candidate_path_vector_tmp.c_path[it_candidates].modified_bases[it_mod_base].second) {
                            really_modified = true;
                        }
                    }

                    if (really_modified == true) {
                        // check the index of the first modified base
                        // extension is needed
                        if (candidate_path_vector_tmp.c_path[it_candidates].modified_bases[0].first < (arg.kmer_length - 1)) {
                            bool extension_success(false);
                            solid_first_kmer(
                                    candidate_path_vector_tmp.c_path[it_candidates],
                                    sequence_modified,
                                    extension_success,
                                    bit_vector,
                                    hash_seed,
                                    arg,
                                    funcm
                                    );

                            if (extension_success == true) {
                                push_back(candidate_path_vector_tmp_tmp, candidate_path_vector_tmp.c_path[it_candidates]);
                            }
                        }
                        // extension is not needed
                        else {
                            push_back(candidate_path_vector_tmp_tmp, candidate_path_vector_tmp.c_path[it_candidates]);
                        }
                    }
                }
            }
        }

        // candidiate_path_vector_tmp_tmp: solid k-mers in candidate_path_vector_tmp

        // candidates in candidiate_path_vector_tmp_tmp are moved to candidate_path_vector_tmp
        //		candidate_path_vector_tmp = candidate_path_vector_tmp_tmp;
        //		memcpy(&candidate_path_vector_tmp, &candidate_path_vector_tmp_tmp, sizeof(vector_candidate_path));

        vector_candidate_path_assignment(&candidate_path_vector_tmp, &candidate_path_vector_tmp_tmp);
        clear(candidate_path_vector_tmp_tmp);

        //--------------------------------------------------
        // STEP 2-2: extend candidate paths to the right
        //--------------------------------------------------
        if (candidate_path_vector_tmp.size > 0) {
            // each path
            bool run_exploration(true);

            for (uint32_t it_candidates = 0; it_candidates < candidate_path_vector_tmp.size; it_candidates++) {
                if (run_exploration == true) {
                    extend_first_kmer_to_right(
                            sequence_modified,
                            read_quality,
                            candidate_path_vector_tmp.c_path[it_candidates],
                            candidate_path_vector_tmp_tmp,
                            read_length,
                            bit_vector,
                            hash_seed,
                            run_exploration,
                            arg,
                            funcm
                            );
                }
            }

            // complete exploration was not done because there are too many candidata paths
            // remove all the paths in candidate_path_vector_tmp
            if (run_exploration == false) {
                clear(candidate_path_vector_tmp_tmp);
            }
        }


        // candidiate_path_vector_tmp_tmp: successfully right extended candidates

        // remain only really modified paths
        clear(candidate_path_vector_tmp);
        // each path
        for (uint32_t it_candidate = 0; it_candidate < candidate_path_vector_tmp_tmp.size; it_candidate++) {
            // each modification
            bool really_modified(false);
            for (unsigned int it_mod_base = 0; it_mod_base < candidate_path_vector_tmp_tmp.c_path[it_candidate].size; it_mod_base++) {
                if (read_content[candidate_path_vector_tmp_tmp.c_path[it_candidate].modified_bases[it_mod_base].first] != candidate_path_vector_tmp_tmp.c_path[it_candidate].modified_bases[it_mod_base].second) {
                    really_modified = true;
                }
            }

            if (really_modified) {
                push_back(candidate_path_vector_tmp, candidate_path_vector_tmp_tmp.c_path[it_candidate]);
            }
        }

        //		candidate_path_vector_tmp_tmp = candidate_path_vector_tmp;
        vector_candidate_path_assignment(&candidate_path_vector_tmp_tmp, &candidate_path_vector_tmp);


        //--------------------------------------------------
        // STEP 2-3: choose a final one in candidate_path_vector_tmp_tmp if possible
        //--------------------------------------------------
        // compare quality scores of candidate paths
        // if the number of paths in candidate_path_vector_tmp_tmp is larger than 1
        if (candidate_path_vector_tmp_tmp.size > 1) {
            // each path
            C_candidate_path* it_path;
            C_candidate_path* it_path_1st;
            C_candidate_path* it_path_2nd;

            uint32_t qs_1st = (INIT_MIN_QS);
            uint32_t qs_2nd = (INIT_MIN_QS);

            // find the first and second paths
            // each candidate path
            for (it_path = candidate_path_vector_tmp_tmp.c_path; it_path != candidate_path_vector_tmp_tmp.c_path + candidate_path_vector_tmp_tmp.size; it_path++) {
                // each modification
                for (uint32_t it_mod = 0; it_mod < (*it_path).size; it_mod++) {
                    // add quality scores of modified bases
                    if (sequence_modified[(*it_path).modified_bases[it_mod].first] != (*it_path).modified_bases[it_mod].second) {
                        (*it_path).sum_qs += ((unsigned short int)read_quality[(*it_path).modified_bases[it_mod].first] - arg.quality_score_offset);
                    }
                }

                // compare quality scores of each path
                if ((*it_path).sum_qs <= qs_1st) {
                    qs_2nd = qs_1st;
                    qs_1st = (*it_path).sum_qs;

                    it_path_2nd = it_path_1st;
                    it_path_1st = it_path;
                }
                else if ((*it_path).sum_qs <= qs_2nd) {
                    qs_2nd = (*it_path).sum_qs;

                    it_path_2nd = it_path;
                }
            }

            // check whether too many bases are modified in the first path, which implies indel may exist in the original read
            // if too many modifications exist, the read is just trimmedsequence_modification
            //bool keep_going(true);
            bool too_many_corrections(false);

            // at least over MAX_MODIFICATION times of modifications from the first modified position
            if (((read_length - (*it_path_1st).modified_bases[0].first) >= min_check_length) &&
                    ((*it_path_1st).size > MAX_MODIFICATION)) {
                // each modified base in a right range
                //std::uint32_t partial_qs((*it_path_1st).sum_qs);

                for (unsigned int it_mod_base = 0; it_mod_base < ((*it_path_1st).size - MAX_MODIFICATION); it_mod_base++) {
                    // at least MAX_MODIFICATION times of modifications within min_check_length bases
                    if (((*it_path_1st).modified_bases[it_mod_base + MAX_MODIFICATION].first - (*it_path_1st).modified_bases[it_mod_base].first) < min_check_length) {
                        //if (it_mod_base > 0) {
                        //   partial_qs -= ((unsigned short int)quality_score[(*it_path_1st).modified_bases[it_mod_base - 1].first] - quality_score_offset);
                        //}

                        // average quality score of the modified bases is too high
                        //if (1.0 * partial_qs / ((*it_path_1st).modified_bases.size() - it_mod_base) > quality_score_cutoff) {
                        // trim_5_end or trim_3_end
                        if ((read_length - (*it_path_1st).modified_bases[it_mod_base].first) <= max_trimmed_bases) {
                            trim_3_end = read_length - (*it_path_1st).modified_bases[it_mod_base].first;
                        }
                        else if (((*it_path_1st).modified_bases[it_mod_base].first + 1) <= max_trimmed_bases) {
                            trim_5_end = (*it_path_1st).modified_bases[it_mod_base].first + 1;
                        }

                        too_many_corrections = true;
                        break;
                        //}
                    }
                }
            }

            // use the 1st path if not too many corrections are made AND if the 1st path has a sufficiently low score
            // if an indel exists using quality scores is not a good way to choose the best path
            if ((too_many_corrections == false) && (qs_1st + MIN_QS_DIFF <= qs_2nd)) {
                for (unsigned int it_base = 0; it_base < (*it_path_1st).size; it_base++) {
                    // filter out the bases that are equal to the original ones
                    if (sequence_modified[(*it_path_1st).modified_bases[it_base].first] != (*it_path_1st).modified_bases[it_base].second) {
                        sequence_modification[(*it_path_1st).modified_bases[it_base].first] = (*it_path_1st).modified_bases[it_base].second;
                        sequence_modified[(*it_path_1st).modified_bases[it_base].first] = (*it_path_1st).modified_bases[it_base].second;

                        num_corrected_errors_local++;
                    }
                }
            }
            // hard to choose one path
            else {
                // A AND B
                vector_pair_uint_char base_vector_intersection;
                base_vector_intersection.size = 0;
                base_vector_intersection.sum_qs = 0;
                // A OR B
                vector_pair_uint_char base_vector_union;

                // temporary vectors
                vector_pair_uint_char base_vector_intersection_prev;
                memcpy(&base_vector_intersection_prev, &(candidate_path_vector_tmp_tmp.c_path[0]), sizeof(vector_pair_uint_char));
                vector_pair_uint_char base_vector_union_prev;
                //				(candidate_path_vector_tmp_tmp[0].modified_bases);
                memcpy(&base_vector_union_prev, &(candidate_path_vector_tmp_tmp.c_path[0]), sizeof(vector_pair_uint_char));

                // each candidate path
                for (uint32_t it_p = 1; it_p < candidate_path_vector_tmp_tmp.size; it_p++) {
                    clear(base_vector_intersection);
                    clear(base_vector_union);

                    base_intersection(base_vector_intersection_prev.uint_char, base_vector_intersection_prev.uint_char + base_vector_intersection_prev.size,
                            candidate_path_vector_tmp_tmp.c_path[it_p].modified_bases, candidate_path_vector_tmp_tmp.c_path[it_p].modified_bases + candidate_path_vector_tmp_tmp.c_path[it_p].size, base_vector_intersection);
                    base_union       (base_vector_union_prev.uint_char,        base_vector_union_prev.uint_char + base_vector_union_prev.size,
                            candidate_path_vector_tmp_tmp.c_path[it_p].modified_bases, candidate_path_vector_tmp_tmp.c_path[it_p].modified_bases+candidate_path_vector_tmp_tmp.c_path[it_p].size, base_vector_union);

                    base_vector_intersection_prev = base_vector_intersection;
                    base_vector_union_prev        = base_vector_union;
                }

                // A - B
                vector_pair_uint_char base_vector_difference;
                base_vector_difference.size = 0;
                base_vector_difference.sum_qs = 0;
                base_difference(base_vector_union.uint_char, base_vector_union.uint_char + base_vector_union.size, base_vector_intersection.uint_char,base_vector_intersection.uint_char + base_vector_intersection.size, base_vector_difference);

                // find trimmed region
                // correcting the 5'-end and 3'-end was not done
                // therefore the total number of trimmed bases is 0 yet
                if (base_vector_difference.size > 0) {
                    uint32_t vector_index_leftmost = (0);
                    uint32_t vector_index_rightmost = (base_vector_difference.size - 1);

                    bool keep_going(true);

                    while (keep_going) {
                        // # of trimmed bases at the 5' end: base_vector_difference[vector_index_leftmost].first + 1
                        // # of trimmed bases at the 3' end: read_length - base_vector_difference[vector_index_rightmost].first
                        // the 5'-end is smaller
                        if ((base_vector_difference.uint_char[vector_index_leftmost].first + 1) < (read_length - base_vector_difference.uint_char[vector_index_rightmost].first)) {
                            // check the total number of trimmed bases
                            if ((base_vector_difference.uint_char[vector_index_leftmost].first + 1 + trim_3_end) <= max_trimmed_bases) {
                                trim_5_end = base_vector_difference.uint_char[vector_index_leftmost].first + 1;

                                // two points are met
                                if (vector_index_leftmost == vector_index_rightmost) {
                                    keep_going = false;
                                }
                                else {
                                    vector_index_leftmost++;
                                }
                            }
                            // no need for more check
                            else {
                                keep_going = false;
                            }
                        }
                        // the 3'-end is smaller
                        else {
                            // check the total number of trimmed bases
                            if ((read_length - base_vector_difference.uint_char[vector_index_rightmost].first) <= max_trimmed_bases) {
                                trim_3_end = read_length - base_vector_difference.uint_char[vector_index_rightmost].first;

                                // two points are met
                                if (vector_index_leftmost == vector_index_rightmost) {
                                    keep_going = false;
                                }
                                else {
                                    vector_index_rightmost--;
                                }
                            }
                            // no need for more check
                            else {
                                keep_going = false;
                            }
                        }
                    }
                }

                // find consensus modifications
                for (uint32_t it_inter = 0; it_inter < base_vector_intersection.size; it_inter++) {
                    // check whether the base is not in the trimmed regions
                    if ((base_vector_intersection.uint_char[it_inter].first <  (read_length - trim_3_end)) &&
                            (base_vector_intersection.uint_char[it_inter].first >= trim_5_end)) {
                        // filter out the bases that are equal to the original ones
                        if (sequence_modified[base_vector_intersection.uint_char[it_inter].first] != base_vector_intersection.uint_char[it_inter].second) {
                            sequence_modification[base_vector_intersection.uint_char[it_inter].first] = base_vector_intersection.uint_char[it_inter].second;
                            sequence_modified[base_vector_intersection.uint_char[it_inter].first] = base_vector_intersection.uint_char[it_inter].second;

                            num_corrected_errors_local++;
                        }
                    }
                }
            }
        }
        // only one path
        // correction succeeds
        else if (candidate_path_vector_tmp_tmp.size == 1) {
            // check whether too many bases are modified in the first path, which implies indel may exist in the original read
            // if too many modifications exist, the read is just trimmed
            //bool keep_going(true);
            bool too_many_corrections(false);

            // at least over MAX_MODIFICATION times of modifications from the first modified position
            if (((read_length - candidate_path_vector_tmp_tmp.c_path[0].modified_bases[0].first) >= min_check_length) &&
                    (candidate_path_vector_tmp_tmp.c_path[0].size > MAX_MODIFICATION)) {
                // calculate sum_qs
                //for (unsigned int it_mod_base = 0; it_mod_base < candidate_path_vector_tmp_tmp[0].modified_bases.size(); it_mod_base++) {
                //   candidate_path_vector_tmp_tmp[0].sum_qs += ((unsigned short int)quality_score[candidate_path_vector_tmp_tmp[0].modified_bases[it_mod_base].first] - quality_score_offset);
                //}

                // each modified base in a right range
                //std::uint32_t partial_qs(candidate_path_vector_tmp_tmp[0].sum_qs);

                for (unsigned int it_mod_base = 0; it_mod_base < (candidate_path_vector_tmp_tmp.c_path[0].size - MAX_MODIFICATION); it_mod_base++) {
                    // at least MAX_MODIFICATION times of modifications within min_check_length bases
                    if ((candidate_path_vector_tmp_tmp.c_path[0].modified_bases[it_mod_base + MAX_MODIFICATION].first - candidate_path_vector_tmp_tmp.c_path[0].modified_bases[it_mod_base].first) < min_check_length) {
                        //if (it_mod_base > 0) {
                        //   partial_qs -= ((unsigned short int)quality_score[candidate_path_vector_tmp_tmp[0].modified_bases[it_mod_base - 1].first] - quality_score_offset);
                        //}

                        // average quality score of the modified bases is too high
                        //if (1.0 * partial_qs / (candidate_path_vector_tmp_tmp[0].modified_bases.size() - it_mod_base) > quality_score_cutoff) {
                        // trim_5_end or trim_3_end
                        if ((read_length - candidate_path_vector_tmp_tmp.c_path[0].modified_bases[it_mod_base].first) <= max_trimmed_bases) {
                            trim_3_end = read_length - candidate_path_vector_tmp_tmp.c_path[0].modified_bases[it_mod_base].first;

                            // update sequence_modification for the non-trimmed corrections
                            if (it_mod_base > 0) {
                                for (unsigned int it_base = 0; it_base < (it_mod_base - 1); it_base++) {
                                    sequence_modification[(candidate_path_vector_tmp_tmp.c_path[0]).modified_bases[it_base].first] = (candidate_path_vector_tmp_tmp.c_path[0]).modified_bases[it_base].second;
                                    sequence_modified[(candidate_path_vector_tmp_tmp.c_path[0]).modified_bases[it_base].first] = (candidate_path_vector_tmp_tmp.c_path[0]).modified_bases[it_base].second;
                                    num_corrected_errors_local++;
                                }
                            }
                        }
                        else if ((candidate_path_vector_tmp_tmp.c_path[0].modified_bases[it_mod_base].first + 1) <= max_trimmed_bases) {
                            trim_5_end = candidate_path_vector_tmp_tmp.c_path[0].modified_bases[it_mod_base].first + 1;

                            // update sequence_modification for the non-trimmed corrections
                            if (it_mod_base < (candidate_path_vector_tmp_tmp.c_path[0].size - 1)) {
                                for (unsigned int it_base = (it_mod_base + 1); it_base < candidate_path_vector_tmp_tmp.c_path[0].size; it_base++) {
                                    sequence_modification[(candidate_path_vector_tmp_tmp.c_path[0]).modified_bases[it_base].first] = (candidate_path_vector_tmp_tmp.c_path[0]).modified_bases[it_base].second;
                                    sequence_modified[(candidate_path_vector_tmp_tmp.c_path[0]).modified_bases[it_base].first] = (candidate_path_vector_tmp_tmp.c_path[0]).modified_bases[it_base].second;
                                    num_corrected_errors_local++;
                                }
                            }
                        }

                        too_many_corrections = true;
                        break;
                        //}
                    }
                }
            }

            // not too many modified bases
            if (too_many_corrections == false) {
                for (unsigned int it_base = 0; it_base < candidate_path_vector_tmp_tmp.c_path[0].size; it_base++) {
                    // filter out the bases that are equal to the original ones
                    if (sequence_modified[candidate_path_vector_tmp_tmp.c_path[0].modified_bases[it_base].first] != candidate_path_vector_tmp_tmp.c_path[0].modified_bases[it_base].second) {
                        sequence_modification[candidate_path_vector_tmp_tmp.c_path[0].modified_bases[it_base].first] = candidate_path_vector_tmp_tmp.c_path[0].modified_bases[it_base].second;
                        sequence_modified[candidate_path_vector_tmp_tmp.c_path[0].modified_bases[it_base].first] = candidate_path_vector_tmp_tmp.c_path[0].modified_bases[it_base].second;

                        num_corrected_errors_local++;
                    }
                }
            }
        }
        // no path
        else {
            // there are solid regions and they were ignored because the length among them is too short
            // okay, trust them
            if (solid_regions.size > 0) {
                trim_5_end = solid_regions.uint_uint[0].first;

                // trim as many bases as possible
                for (uint32_t it_solid_short = 0; it_solid_short < solid_regions.size; it_solid_short++) {
                    if ((trim_5_end + (read_length - arg.kmer_length - solid_regions.uint_uint[it_solid_short].second)) <= max_trimmed_bases) {
                        trim_3_end = read_length - arg.kmer_length - solid_regions.uint_uint[it_solid_short].second;
                        break;
                    }
                }
            }
            // the first k-mer has at least one error
            else if (arg.kmer_length <= max_trimmed_bases) {
                trim_5_end = arg.kmer_length;
            }
        }

        free(candidate_path_vector_tmp.c_path);
        free(candidate_path_vector_tmp_tmp.c_path);

    }

    //	free(sequence_modified);
    funcm.sequence_ind --;
    //*/
}



void correct_errors(init_args arg,  C_arg& c_inst_args,uint8_t* bit_vector) {
    //	unsigned int* hash_seed = (unsigned int*) malloc(
    //			arg.num_hash_func * sizeof(int));
    //	generate_hash_seed(arg.random_seed, arg.num_hash_func, hash_seed);

    char* read_name = (char*) malloc(BLOCK_READ_NUM * READ_MAX_LEN);
    char* read_plus = (char*) malloc(BLOCK_READ_NUM * READ_MAX_LEN);
    char* read_content = (char*) malloc(BLOCK_READ_NUM * READ_MAX_LEN);
    char* read_quality = (char*) malloc(BLOCK_READ_NUM * READ_MAX_LEN);
    char* read_modification = (char*) malloc(BLOCK_READ_NUM * READ_MAX_LEN);
    uint32_t* readns = (uint32_t*) malloc(BLOCK_READ_NUM * sizeof(uint32_t));
    memset(read_name, 0, BLOCK_READ_NUM * READ_MAX_LEN);
    memset(read_plus, 0, BLOCK_READ_NUM * READ_MAX_LEN);
    memset(read_content, 0, BLOCK_READ_NUM * READ_MAX_LEN);
    memset(read_quality, 0, BLOCK_READ_NUM * READ_MAX_LEN);
    memset(read_modification, 0, BLOCK_READ_NUM * READ_MAX_LEN);
    memset(readns, 0, BLOCK_READ_NUM * sizeof(uint32_t));

    //	char file_name[100];
    //	sprintf(file_name, "%s", );
    //	std::cout << c_inst_args.read_file_name1 << std::endl;
    FILE* pfile = fopen(c_inst_args.read_file_name1.c_str(), "r");
    std::cout << c_inst_args.read_file_name1 << std::endl;


#pragma offload target(mic:0)\
    nocopy(read_content:length(BLOCK_READ_NUM * READ_MAX_LEN) alloc_if(1) free_if(0))\
    nocopy(read_modification:length(BLOCK_READ_NUM * READ_MAX_LEN) alloc_if(1) free_if(0))\
    nocopy(read_quality:length(BLOCK_READ_NUM * READ_MAX_LEN) alloc_if(1) free_if(0))\
    nocopy(readns:length(BLOCK_READ_NUM) alloc_if(1) free_if(0))\
    in(bit_vector:length(arg.bit_vector_width_byte) alloc_if(1) free_if(0))
    {

    }
    uint32_t read_num = 1;
    uint32_t read_length = 101;

    int  bit_vector_width_byte = arg.bit_vector_width_byte;
   int total_time = 0;
    while (read_num) {
        read_fastq(pfile, read_num, read_name, read_plus, read_content,
                read_quality);
        std::cout << "read num == " << read_num << std::endl;

#pragma omp parallel for
        for (uint32_t i = 0; i < read_num; i++) {
            readns[i] = 0;
            for (uint32_t j = 0; j < read_length; j++) {
                read_content[i * READ_MAX_LEN + j] = toupper(
                        read_content[i * READ_MAX_LEN + j]);
                char c = read_content[i * READ_MAX_LEN + j];
                if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
                    read_content[i * READ_MAX_LEN + j] = 'A';
                    readns[i] += 1;
                }
            }
        }

        time_t starttime;
        time_t endtime;
        starttime = time(NULL);
        arg.read_num = read_num;
        arg.read_len = read_length;


#pragma offload target(mic:0)\
        in(read_content:length(BLOCK_READ_NUM * READ_MAX_LEN) alloc_if(0) free_if(0) )\
        in(read_quality:length(BLOCK_READ_NUM * READ_MAX_LEN) alloc_if(0) free_if(0) )\
        in(read_modification:length(BLOCK_READ_NUM * READ_MAX_LEN) alloc_if(0) free_if(0) )\
        in(readns:length(BLOCK_READ_NUM) alloc_if(0) free_if(0) )\
		in(bit_vector:length(arg.bit_vector_width_byte) alloc_if(0) free_if(0) )\
		in(arg:length(1))\
		in(read_num)
        {


            uint32_t hash_seed[6];
            hash_seed[0] = 1804289383;
            hash_seed[1] = 846930886;
            hash_seed[2] = 1681692777;
            hash_seed[3] = 1714636915;
            hash_seed[4] = 1957747793;
            hash_seed[5] = 424238335;
            printf("%d %d %d %d %d %d %d %d %d \n", arg.bit_vector_width, arg.num_hash_func, arg.kmer_length,
            		arg.read_num, arg.read_len, arg.quality_score_offset, arg.quality_score_cutoff, arg.extremely_low_quality_score, arg.max_extension);
            printf("sequence %s \n", read_content);
            printf("quality_score %s \n", read_quality);


            func_mem funcm[240];
            for (int i = 0; i < 240; i++) {
                funcm[i].kmer_length = arg.kmer_length + 5;
                funcm[i].kmer = (char*) _mm_malloc(
                        (2000  + 2)* (arg.kmer_length + 5) * sizeof(char), 64);
                funcm[i].kmer_size = 2000;
                funcm[i].sequence = (char*)malloc(100*READ_MAX_LEN * sizeof(char));
                funcm[i].sequence_size = 100;
//                funcm[i].candidate_path = (C_candidate_path*)malloc(200 * sizeof(C_candidate_path));
//                funcm[i].candidate_path_size = 200;
//                funcm[i].candidate_path_ind = 0;
                funcm[i].kmer_index = (uint32_t*)_mm_malloc(80 * sizeof(uint32_t), 64);
            }
            uint32_t num_trimmed_bases_tmp = 0;
            uint32_t num_corrected_errors_local = 0;
            uint32_t trim_5_end = 0;
            uint32_t trim_3_end = 0;
            uint32_t max_trimmed_bases = 0;
            uint32_t min_check_length = arg.read_len * CHECK_RANGE_RATIO;
            uint32_t max_allowed_ns = read_length * MAX_N_RATIO;
            bool too_many_errors;
            uint32_t num_corrected_errors_tmp = 0;
            uint32_t num_corrected_reads_tmp = 0;


#pragma omp parallel for schedule(dynamic, 10)  private(read_length, min_check_length, max_allowed_ns,\
        max_trimmed_bases, num_corrected_errors_local, trim_5_end, trim_3_end, too_many_errors,\
		num_corrected_errors_tmp, num_corrected_reads_tmp, num_trimmed_bases_tmp) num_threads(120)
            for (uint32_t i = 0; i < read_num; i++) {

                read_length = 101;
                min_check_length = read_length * CHECK_RANGE_RATIO;
                max_allowed_ns = read_length * MAX_N_RATIO;
                // forward
                // too short read: no trimming
                if ( read_length <= MIN_BASES_AFTER_TRIMMING) {
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
                funcm[ind].kmer_ind = 0;
                funcm[ind].sequence_ind = 0;
                //			&& i != 1415699

                //			printf("read_content %s \n", read_content + i * READ_MAX_LEN);
                //			printf("read_quality %s \n", read_quality + i * READ_MAX_LEN);
                if(readns[i] < max_allowed_ns){
                    mic_correct_errors(read_content + i * READ_MAX_LEN,
                            read_quality + i * READ_MAX_LEN,
                            read_modification + i * READ_MAX_LEN, hash_seed, bit_vector,
                            arg, num_corrected_errors_local, trim_5_end, trim_3_end,
                            max_trimmed_bases, min_check_length, read_length, i,
                            funcm[ind]);
                    //			printf("%d done\n", i);
                }
                //		printf("%d\n", i);
                //		if(i % 100000 == 0)
                //			printf("right %d \n", i);
                /*
                // no trim
                if (c_inst_args.notrim == true) {
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
                else if (num_corrected_errors_local > 0) {sequence_modification
                num_corrected_errors_tmp += num_corrected_errors_local;

                num_corrected_reads_tmp++;
                }
                else if (c_inst_args.notrim == false) {
                if ((trim_5_end > 0) || (trim_3_end > 0)) {
                num_corrected_reads_tmp++;
                }
                }

                // make a corrected read
                if (too_many_errors == false) {
                // apply modifications to the read
                for (unsigned int it_base = trim_5_end; it_base < (read_length - trim_3_end); it_base++) {
                if (read_modification[READ_MAX_LEN * i + it_base] != '0') {
                read_content[READ_MAX_LEN * i + it_base] = read_modification[READ_MAX_LEN * i + it_base];
                }
                }
                }

                // make a trimmed read
                if ((trim_5_end + trim_3_end) > 0) {
                memcpy(&read_content[READ_MAX_LEN * i], &read_content[READ_MAX_LEN * i  + trim_5_end], read_length);
                read_content[READ_MAX_LEN * i + read_length] ='\n';
                }
                else {
                read_content[READ_MAX_LEN * i + read_length] ='\n';
                }

                //		std::cout << "hhhhhhhhhhhhh" << std::endl;
                //		printf("%s", read_content);
                //	}
                //i
                */
        }
        for (int i = 0; i < 240; i++) {
            _mm_free(funcm[i].kmer);
            funcm[i].kmer_size = 0;
            free(funcm[i].sequence);
            funcm[i].sequence_size = 0;
            _mm_free(funcm[i].kmer_index);
        }


    }
    endtime = time(NULL);
    std::cout << "time == " << endtime - starttime << std::endl;
    total_time += endtime - starttime;
}

fclose(pfile);
printf("total time %d \n", total_time);
free(read_name);
free(read_quality);
free(read_plus);
free(read_content);
free(read_modification);
//	free(hash_seed);
}
