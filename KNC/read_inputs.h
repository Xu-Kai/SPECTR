#ifndef READ_INTPUTS_H_
#define READ_INPUTS_H_

//----------------------------------------------------------------------
// generate_hash_seed
//----------------------------------------------------------------------
void generate_hash_seed(const uint32_t random_seed, uint32_t num_hash_func,
        unsigned int* hash_seed); 
void read_fastq(FILE* pfile, uint32_t& read_num, char* read_name, char* read_plus,
        char* read_content, char* read_quality);
#endif
