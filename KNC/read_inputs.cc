#include "define.h"
#include "read_inputs.h"
#define READ_MAX_LEN 256
#define BLOCK_READ_NUM 2000000
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


