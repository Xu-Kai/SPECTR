set -x
#./vbf -max_mem 64 -prefix fixed/fragfix -readlist cluster.list -kmerlength 31  -load fixed/fragfix
i=31
./xbless-knl -max_mem 32 -prefix fixed/fixsrr${i} -readlist ecc.list -kmerlength ${i}
