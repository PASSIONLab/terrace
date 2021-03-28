#!/bin/sh

set -x

make clean
make CILK=1 graph_bm_no_update
./graph_bm_no_update -src 212 -f $1 > gc.out

diff bc.out ligra_out/bc_ligra_lj_212.out
diff bfs.out ligra_out/bfs_ligra_lj_212.out

