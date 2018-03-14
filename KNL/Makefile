#ICC=mpicxx
ICC=mpiicpc
#MPICXX=mpic++
#MIC=
CFLAGS=-g -Wall -c -O3 -fopenmp -I ./kmc/kmc_api  -std=c++11 -xMIC-AVX512
LDFLAGS=-fopenmp -lpthread $(MIC)
SRCS=main_mic.cc mic_correct_errors.cc read_inputs.cc new_query_text.cc new_MurmurHash3.cc\
 parse_args.cc  MurmurHash3.cc check_inputs.cc\
  kmc/kmc_api/kmc_file.cpp kmc/kmc_api/kmer_api.cpp kmc/kmc_api/mmer.cpp count_solid_kmers.cc generate_bloom_filter.cc data_manager.cc 
OBJS=main_mic.o mic_correct_errors.o read_inputs.o new_query_text.o new_MurmurHash3.o\
parse_args.o MurmurHash3.o check_inputs.o\
 kmc/kmc_api/kmc_file.o kmc/kmc_api/kmer_api.o kmc/kmc_api/mmer.o count_solid_kmers.o generate_bloom_filter.o data_manager.o 

 
spectr-knl:$(OBJS)
	$(ICC) $(OBJS) $(LDFLAGS) -o spectr-knl
	
	
main_mic.o: main_mic.cc
	$(ICC) $(CFLAGS) main_mic.cc

new_query_text.o: new_query_text.cc new_query_text.h
	$(ICC) $(CFLAGS) new_query_text.cc

new_MurmurHash3.o: new_MurmurHash3.cc new_MurmurHash3.h
	$(ICC) $(CFLAGS) new_MurmurHash3.cc

read_inputs.o: read_inputs.cc read_inputs.h
	$(ICC) $(CFLAGS) read_inputs.cc


mic_correct_errors.o: mic_correct_errors.cc mic_correct_errors.h
	$(ICC) $(CFLAGS) mic_correct_errors.cc

parse_args.o: parse_args.cc  parse_args.h
	$(ICC) $(CFLAGS) parse_args.cc
	
#query_text.o: query_text.cc  query_text.h
#	$(ICC) $(CFLAGS) query_text.cc

#struct.o: struct.cc  struct.h
#	$(ICC) $(CFLAGS) struct.cc

#sdregion_query_text.o: sdregion_query_text.cc  sdregion_query_text.h
#	$(ICC) $(CFLAGS) sdregion_query_text.cc

MurmurHash3.o: MurmurHash3.cc MurmurHash3.h
	$(ICC) $(CFLAGS) MurmurHash3.cc
	
check_inputs.o: check_inputs.cc check_inputs.h
	$(ICC) $(CFLAGS) check_inputs.cc	
	
count_solid_kmers.o: count_solid_kmers.cc count_solid_kmers.h
	$(ICC) $(CFLAGS) count_solid_kmers.cc	

generate_bloom_filter.o: generate_bloom_filter.cc generate_bloom_filter.h
	$(ICC) $(CFLAGS) generate_bloom_filter.cc

data_manager.o: data_manager.cc data_manager.h
	$(ICC) $(CFLAGS) data_manager.cc
	
	#kmc/kmc_api/kmc_file.o: kmc/kmc_api/kmc_file.cpp kmc/kmc_api/kmc_file.h
#	$(ICC) $(CFLAGS) kmc/kmc_api/kmc_file.cpp	

#kmc/kmc_api/kmer_api.o: kmc/kmc_api/kmer_api.cpp kmc/kmc_api/kmer_api.h
#	$(ICC) $(CFLAGS) kmc/kmc_api/kmer_api.cpp	

#kmc/kmc_api/mmer.o: kmc/kmc_api/mmer.cpp kmc/kmc_api/mmer.h
#	$(ICC) $(CFLAGS) kmc/kmc_api/mmer.cpp	
							
#kmc/kmc_api/%.o: kmc/kmc_api/%.cpp
#	$(ICC) $(CFLAGS) -c $(DEF) $? -o $@


clean:
	rm *.o xbless-knl

