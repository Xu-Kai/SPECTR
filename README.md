# ParSHREC

We present a parallel error correction tool(XBLESS), which is designed to improve the throughput of DNA error correction for Illumina reads. Our design is based on the memory-efficient BLESS algorithm but is optimized towards AVX-512-based CPUs, Xeon Phi many-cores (both KNC and KNL), and heterogeneous compute clusters.

Currently, we presents three implementations for different platform: they are XBLESS-AVX512, XBLESS-KNC, XBLESS-KNL. We implemented them with Intel AVX512, KCI and Intel AVX512 instructuon set respectively. Besides, in order to support efficient processing of large-scale NGS datasets, we have further developed a distributed version of XBLESS using MPI. 

Our experimental results show that XBLESS achieves a speedup of 2.8, 5.2 and 9.3 on a CPU(Xeon W-2123), a KNC-based Xeon Phi(31S1P), and a KNL-based Xeon Phi(7210), respectively. Furthermore, when executed on same hardware, XBLESS achieves a speedup of 1.73, 2.4 and 6.37, compared to the state-of-the-art error correction tools Lighter, RECKONER and Musket. Our MPI cluster version exhibits strong scalability with an efficiency of around 86% when executed on 32 nodes of Tianhe-2 supercomputer.

The techniques presented in XBLESS can also be adapted to map similar applications exhibiting the massively probing of look-up data structures, such as searching of large-scale transcriptomic sequencing databases using sequence Bloom tress or de novo genome assembly onto heterogeneous many-core cluster architectures.


# Getting start

To build XBLESS, after you choose one in three implementation, based on your hardware platform, then simply

do: make

# Running a test 

To run a test using XBLESS, simply do:

$  ./xbless-knl -max_mem 32 -prefix fixed -readlist ecc.list -kmerlength 31

The meaning of augments:

1: -max_mem, the maxmimum memory resource for KMC software (GB)
2: prefix, the prefix for the corrected read file
3: -readlist, the file list of input raw reads
4: kmerlength, the length of kmer used in error correction tool

After the execution is finished, a summary of the run is printed out, the corrected reads will be write in a file according to the prefix.
