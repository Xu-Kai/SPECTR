#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>
#include<sys/stat.h>
void offset_data(char* buffer, uint64_t buffer_size, uint64_t& data_size, uint64_t remain_data_size){
    
    uint32_t n_num = 0;
    uint64_t offset_pos[8];
    if(remain_data_size <= buffer_size)
        data_size = remain_data_size;
    else if(remain_data_size > buffer_size){
            uint64_t i = buffer_size - 1;
                uint64_t j;
                uint32_t plus_num = 0;
                for(j = 0; j < 2048; j ++){
                    if(buffer[i - j] == '\n'){
                        offset_pos[n_num] = j;
                        n_num += 1;
                        if(buffer[i - j + 1] == '+'){
                           plus_num = n_num;
                        }
                        if(n_num - plus_num == 2 && buffer[i - j + 1] == '@'){ 
                           if(n_num >= 5)
                            break;
                        }
                    }
                }
    }
    data_size = buffer_size - offset_pos[n_num - 5];
    /* 
    for(int i = 1000; i > 0; i --){
        printf("%c", buffer[data_size - i]);

    }
    */
    printf("%c", buffer[data_size - 1]);
}
int main(){
    FILE* pfile = fopen("ERR1451602_1.fastq", "r");
    uint32_t iter_count = 0;
    uint64_t total_readed_size = 0;
    char* rbuffer_list[100];
    uint64_t rbuffer_size = 4L * 1024L * 1024L *1024L;
    uint64_t rbuffer_data_size[600];
    fseek(pfile, 0, SEEK_END);
    uint64_t file_len = ftell(pfile);
    fseek(pfile, 0, SEEK_SET);
    uint64_t device_data_size[600];
    uint32_t my_rank = 0;
    while(1){
            rbuffer_list[iter_count] = (char*)malloc(rbuffer_size * sizeof(char));
            rbuffer_data_size[iter_count] = fread(rbuffer_list[iter_count], sizeof(char), rbuffer_size,  pfile);
            total_readed_size += rbuffer_data_size[iter_count];
            iter_count += 1;
            if(total_readed_size >= file_len)
                 break;
    }

    fclose(pfile);
    
            uint64_t device_size;
            uint64_t rbuffer_index = 0;
            uint64_t total_compute_size = 0;
            iter_count = 0;
            uint32_t device_id = 0;
            uint32_t read_flag[600];
            uint64_t device_buffer_size[]={200*1024*1024, 200*1024*1024, 200*1024*1024, 60*1024*1024};
            uint32_t finish_flag = 0;
            uint64_t total_read_size = 0;
            srand(0);
            while(1){
                device_id = rand()%4;
                iter_count += 1;
                 
                if(total_compute_size + device_buffer_size[device_id] < rbuffer_data_size[rbuffer_index] + 1000){
                    uint64_t device_buffer_size_t = rbuffer_data_size[rbuffer_index] - total_compute_size;
                    if(device_buffer_size_t > device_buffer_size[device_id])
                           device_buffer_size_t =  device_buffer_size[device_id];
                    offset_data(rbuffer_list[rbuffer_index] + total_compute_size, device_buffer_size_t, device_data_size[device_id], file_len - total_read_size);
                    printf("device_buffer_size_t %lld\n", device_buffer_size_t);
                    total_compute_size += device_data_size[device_id]; 
                    total_read_size += device_data_size[device_id];
                }else{
                    uint64_t device_buffer_size_t = rbuffer_data_size[rbuffer_index] - total_compute_size;
                    device_data_size[device_id] = device_buffer_size_t;
                    free(rbuffer_list[rbuffer_index]); 
                    total_read_size += device_data_size[device_id];
                    if(total_read_size < file_len){
                        total_compute_size = 0;
                        rbuffer_index += 1;  
                        offset_data(rbuffer_list[rbuffer_index], device_buffer_size[device_id] - device_data_size[device_id], device_buffer_size_t, file_len - total_read_size);
                       total_read_size += device_buffer_size_t;
                       total_compute_size += device_buffer_size_t;
                    }
                }
                printf("%d rank readed size %lld\n", my_rank, total_read_size); 
                
                if(total_read_size >= file_len){
                    finish_flag = 1;
                }
                if(finish_flag == 1){
                    break;
                }
            }

        
    return 0;
}
