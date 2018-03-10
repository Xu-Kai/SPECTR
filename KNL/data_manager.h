#ifndef MIC_DATA_MANAGER_H_
#define MIC_DATA_MANAGER_H_
#include<stdio.h>
#include<stdlib.h>
#include "struct.h"
#include "define.h"
#include <stdint.h>
#include "read_inputs.h"
#include "parse_args.h"
void manage_data(init_args arg, C_time& c_inst_time, C_arg& c_inst_arg, uint8_t* bit_vector, int my_rank, struct stat& sb, FILE* pfile, FILE* f_corrected_read, double& mpi_total_time);
#endif
