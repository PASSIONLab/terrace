#!/bin/bash
#set -x
ROUNDS=10

function ligra_run {
#  printf "ligra, $ROUNDS, $1, $2, $3, 1, "
#  bc <<<"scale=5;$(taskset -c 0 ./ligra/apps/$1 -s  -r $2  -rounds $ROUNDS $3 | awk '{print $NF}' | paste -s -d+ - | bc)/$ROUNDS"
#  printf "ligra, $ROUNDS, $1, $2, $3, 2, "
#  bc <<<"scale=5;$(taskset -c 0-1 ./ligra/apps/$1 -s  -r $2  -rounds $ROUNDS $3 | awk '{print $NF}' | paste -s -d+ - | bc)/$ROUNDS"
#  printf "ligra, $ROUNDS, $1, $2, $3, 4, "
#  bc <<<"scale=5;$(taskset -c 0-3 ./ligra/apps/$1 -s  -r $2  -rounds $ROUNDS $3 | awk '{print $NF}' | paste -s -d+ - | bc)/$ROUNDS"
#  printf "ligra, $ROUNDS, $1, $2, $3, 8, "
#  bc <<<"scale=5;$(taskset -c 0-7 ./ligra/apps/$1 -s  -r $2  -rounds $ROUNDS $3 | awk '{print $NF}' | paste -s -d+ - | bc)/$ROUNDS"
#  printf "ligra, $ROUNDS, $1, $2, $3, 16, "
#  bc <<<"scale=5;$(taskset -c 0-15 ./ligra/apps/$1 -s  -r $2  -rounds $ROUNDS $3 | awk '{print $NF}' | paste -s -d+ - | bc)/$ROUNDS"
#  printf "ligra, $ROUNDS, $1, $2, $3, 24, "
#  bc <<<"scale=5;$(taskset -c 0-23 ./ligra/apps/$1 -s  -r $2  -rounds $ROUNDS $3 | awk '{print $NF}' | paste -s -d+ - | bc)/$ROUNDS"
#  printf "ligra, $ROUNDS, $1, $2, $3, 48h, "
#  bc <<<"scale=5;$(taskset -c 0-23,48-71 ./ligra/apps/$1 -s  -r $2  -rounds $ROUNDS $3 | awk '{print $NF}' | paste -s -d+ - | bc)/$ROUNDS"
#  printf "ligra, $ROUNDS, $1, $2, $3, 48n, "
#  bc <<<"scale=5;$(taskset -c 0-47 ./ligra/apps/$1 -s  -r $2  -rounds $ROUNDS $3 | awk '{print $NF}' | paste -s -d+ - | bc)/$ROUNDS"
# Running command below - prashant
  #printf "ligra, $ROUNDS, $1, $2, $3, 48, "
  #bc <<<"scale=5;$(taskset -c 0-23,48-71 ./ligra/apps/$1 -s  -r $2  -rounds $ROUNDS $3 | grep "Running" | awk '{print $NF}' | paste -s -d+ - | bc)/$ROUNDS"
  printf "ligra, $ROUNDS, $1, $2, $3, 1, "
  bc <<<"scale=5;$(taskset -c 1 ./ligra/apps/$1 -s  -r $2  -rounds $ROUNDS $3 | grep "Running" | awk '{print $NF}' | paste -s -d+ - | bc)/$ROUNDS"
}
function ligra_run_tc {
 # printf "ligra, $ROUNDS, $1, $2, $3, 1, "
 # bc <<<"scale=5;$(taskset -c 0 ./ligra/apps/$1 -s  -r $2  -rounds $ROUNDS $3 | grep time | awk '{print $NF}' | paste -s -d+ - | bc)/$ROUNDS"
 # printf "ligra, $ROUNDS, $1, $2, $3, 2, "
 # bc <<<"scale=5;$(taskset -c 0-1 ./ligra/apps/$1 -s  -r $2  -rounds $ROUNDS $3 | grep time | awk '{print $NF}' | paste -s -d+ - | bc)/$ROUNDS"
 # printf "ligra, $ROUNDS, $1, $2, $3, 4, "
 # bc <<<"scale=5;$(taskset -c 0-3 ./ligra/apps/$1 -s  -r $2  -rounds $ROUNDS $3 | grep time | awk '{print $NF}' | paste -s -d+ - | bc)/$ROUNDS"
 # printf "ligra, $ROUNDS, $1, $2, $3, 8, "
 # bc <<<"scale=5;$(taskset -c 0-7 ./ligra/apps/$1 -s  -r $2  -rounds $ROUNDS $3 | grep time | awk '{print $NF}' | paste -s -d+ - | bc)/$ROUNDS"
 # printf "ligra, $ROUNDS, $1, $2, $3, 16, "
 # bc <<<"scale=5;$(taskset -c 0-15 ./ligra/apps/$1 -s  -r $2  -rounds $ROUNDS $3 | grep time | awk '{print $NF}' | paste -s -d+ - | bc)/$ROUNDS"
 # printf "ligra, $ROUNDS, $1, $2, $3, 24, "
 # bc <<<"scale=5;$(taskset -c 0-23 ./ligra/apps/$1 -s  -r $2  -rounds $ROUNDS $3 | grep time | awk '{print $NF}' | paste -s -d+ - | bc)/$ROUNDS"
 # printf "ligra, $ROUNDS, $1, $2, $3, 48, "
 # bc <<<"scale=5;$(taskset -c 0-23,48-71 ./ligra/apps/$1 -s  -r $2  -rounds $ROUNDS $3 | grep time | awk '{print $NF}' | paste -s -d+ - | bc)/$ROUNDS"
# Running command below - prashant
  #printf "ligra, $ROUNDS, $1, $2, $3, 48, "
  #bc <<<"scale=5;$(taskset -c 0-23,48-71 ./ligra/apps/$1 -s  -r $2  -rounds 10 $3 | grep time | awk '{print $NF}' | paste -s -d+ - | bc)/$ROUNDS"
  printf "ligra, $ROUNDS, $1, $2, $3, 1, "
  bc <<<"scale=5;$(taskset -c 1 ./ligra/apps/$1 -s  -r $2  -rounds 10 $3 | grep time | awk '{print $NF}' | paste -s -d+ - | bc)/$ROUNDS"
}
function ligra_plus_run {
  printf "ligra+, $ROUNDS, $1, $2, $3, 1, "
  bc <<<"scale=5;$(taskset -c 0 ./ligra/apps/$1  -c -s -r $2  -rounds $ROUNDS $3 | tail -n $ROUNDS | awk '{print $NF}' | paste -s -d+ - | bc)/$ROUNDS"
  printf "ligra+, $ROUNDS, $1, $2, $3, 2, "
  bc <<<"scale=5;$(taskset -c 0-1 ./ligra/apps/$1 -c -s -r $2  -rounds $ROUNDS $3 | tail -n $ROUNDS | awk '{print $NF}' | paste -s -d+ - | bc)/$ROUNDS"
  printf "ligra+, $ROUNDS, $1, $2, $3, 4, "
  bc <<<"scale=5;$(taskset -c 0-3 ./ligra/apps/$1  -c -s -r $2  -rounds $ROUNDS $3 | tail -n $ROUNDS | awk '{print $NF}' | paste -s -d+ - | bc)/$ROUNDS"
  printf "ligra+, $ROUNDS, $1, $2, $3, 8, "
  bc <<<"scale=5;$(taskset -c 0-7 ./ligra/apps/$1 -c -s -r $2  -rounds $ROUNDS $3 | tail -n $ROUNDS | awk '{print $NF}' | paste -s -d+ - | bc)/$ROUNDS"
  printf "ligra+, $ROUNDS, $1, $2, $3, 16, "
  bc <<<"scale=5;$(taskset -c 0-15 ./ligra/apps/$1 -c -s -r $2  -rounds $ROUNDS $3 | tail -n $ROUNDS | awk '{print $NF}' | paste -s -d+ - | bc)/$ROUNDS"
  printf "ligra+, $ROUNDS, $1, $2, $3, 24, "
  bc <<<"scale=5;$(taskset -c 0-23 ./ligra/apps/$1 -c -s -r $2  -rounds $ROUNDS $3 | tail -n $ROUNDS | awk '{print $NF}' | paste -s -d+ - | bc)/$ROUNDS"
  printf "ligra+, $ROUNDS, $1, $2, $3, 48h, "
  bc <<<"scale=5;$(taskset -c 0-23,48-71 ./ligra/apps/$1  -c -s -r $2  -rounds $ROUNDS $3 | tail -n $ROUNDS | awk '{print $NF}' | paste -s -d+ - | bc)/$ROUNDS"
  printf "ligra+, $ROUNDS, $1, $2, $3, 48n, "
  bc <<<"scale=5;$(taskset -c 0-47 ./ligra/apps/$1  -c -s -r $2  -rounds $ROUNDS $3 | tail -n $ROUNDS | awk '{print $NF}' | paste -s -d+ - | bc)/$ROUNDS"
  printf "ligra+, $ROUNDS, $1, $2, $3, 96, "
  bc <<<"scale=5;$(numactl -i all  ./ligra/apps/$1  -r $2 -c -s -rounds $ROUNDS $3 | tail -n $ROUNDS | awk '{print $NF}' | paste -s -d+ - | bc)/$ROUNDS"
}
function aspen_run {
 # printf "aspen, $ROUNDS, $1, $2, $3, 1, "
 # taskset -c 0 ./aspen/code/run_static_algorithm -s -t $1 -src $2 -f $3 -rounds $ROUNDS | grep "RESULT (AVG)" | awk '{print $5}'
 # printf "aspen, $ROUNDS, $1, $2, $3, 2, "
 # taskset -c 0-1 ./aspen/code/run_static_algorithm -s -t $1 -src $2 -f $3 -rounds $ROUNDS | grep "RESULT (AVG)" | awk '{print $5}'
 # printf "aspen, $ROUNDS, $1, $2, $3, 4, "
 # taskset -c 0-3 ./aspen/code/run_static_algorithm -s -t $1 -src $2 -f $3 -rounds $ROUNDS | grep "RESULT (AVG)" | awk '{print $5}'
 # printf "aspen, $ROUNDS, $1, $2, $3, 8, "
 # taskset -c 0-7 ./aspen/code/run_static_algorithm -s -t $1 -src $2 -f $3 -rounds $ROUNDS | grep "RESULT (AVG)" | awk '{print $5}'
 # printf "aspen, $ROUNDS, $1, $2, $3, 16, "
 # taskset -c 0-15 ./aspen/code/run_static_algorithm -s -t $1 -src $2 -f $3 -rounds $ROUNDS | grep "RESULT (AVG)" | awk '{print $5}'
 # printf "aspen, $ROUNDS, $1, $2, $3, 24, "
 # taskset -c 0-23 ./aspen/code/run_static_algorithm -s -t $1 -src $2 -f $3 -rounds $ROUNDS | grep "RESULT (AVG)" | awk '{print $5}'
 # printf "aspen, $ROUNDS, $1, $2, $3, 48h, "
 # taskset -c 0-23,48-71 ./aspen/code/run_static_algorithm -s -t $1 -src $2 -f $3 -rounds $ROUNDS | grep "RESULT (AVG)" | awk '{print $5}'
 # printf "aspen, $ROUNDS, $1, $2, $3, 48n, "
 # taskset -c 0-47 ./aspen/code/run_static_algorithm -s -t $1 -src $2 -f $3 -rounds $ROUNDS | grep "RESULT (AVG)" | awk '{print $5}'
# Running command below - prashant
  #printf "aspen, $ROUNDS, $1, $2, $3, 48, "
  #taskset -c 0-23,48-71 ./aspen/code/run_static_algorithm -s -t $1 -src $2 -f $3 -rounds $ROUNDS | grep "RESULT (AVG)" | awk '{print $5}'
  printf "aspen, $ROUNDS, $1, $2, $3, 1, "
  taskset -c 1 ./aspen/code/run_static_algorithm -s -t $1 -src $2 -f $3 -rounds $ROUNDS | grep "RESULT (AVG)" | awk '{print $5}'
}


ligra_run BFS 9 /graphs/soc-LiveJournal_shuf.adj
ligra_run BFS 28 /graphs/com-orkut.ungraph.adj.shuf
ligra_run BFS 19372 /graphs/rmat_ligra.adj.shuf
#ligra_run BFS 0 ~/graphs/soc-LiveJournal1_sym.adj
#ligra_run BFS 1000 ~/graphs/com-orkut.ungraph.adj
#ligra_run BFS 0 ~/graphs/rmat_ligra.adj
#ligra_run BFS 12 ~/graphs/twitter.adj
ligra_run BFS 6662945 /graphs/twitter.adj.shuf
ligra_run BFS 35 /graphs/protein.adj

#ligra_plus_run BFS 9 ~/graphs/soc-LiveJournal_shuf.adj.compressed
#ligra_plus_run BFS 28 ~/graphs/com-orkut.ungraph.adj.shuf.compressed
#ligra_plus_run BFS 19372 ~/graphs/rmat_ligra.adj.shuf.compressed
#ligra_plus_run BFS 0 ~/graphs/soc-LiveJournal1_sym.adj.compressed
#ligra_plus_run BFS 1000 ~/graphs/com-orkut.ungraph.adj.compressed
#ligra_plus_run BFS 0 ~/graphs/rmat_ligra.adj.compressed
##ligra_plus_run BFS 12 ~/graphs/twitter.adj.compressed
##ligra_plus_run BFS 6662945 ~/graphs/twitter.adj.shuf.compressed
#ligra_plus_run BFS 35 ~/graphs/protein.adj.compressed

#aspen_run BFS 9 /graphs/soc-LiveJournal_shuf.adj
#aspen_run BFS 28 /graphs/com-orkut.ungraph.adj.shuf
#aspen_run BFS 19372 /graphs/rmat_ligra.adj.shuf
#aspen_run BFS 0 ~/graphs/soc-LiveJournal1_sym.adj
#aspen_run BFS 1000 ~/graphs/com-orkut.ungraph.adj
#aspen_run BFS 0 ~/graphs/rmat_ligra.adj
#aspen_run BFS 12 ~/graphs/twitter.adj
#aspen_run BFS 6662945 /graphs/twitter.adj.shuf
#aspen_run BFS 35 /graphs/protein.adj

ligra_run BC 9 /graphs/soc-LiveJournal_shuf.adj
ligra_run BC 28 /graphs/com-orkut.ungraph.adj.shuf
ligra_run BC 19372 /graphs/rmat_ligra.adj.shuf
#ligra_run BC 0 ~/graphs/soc-LiveJournal1_sym.adj
#ligra_run BC 1000 ~/graphs/com-orkut.ungraph.adj
#ligra_run BC 0 ~/graphs/rmat_ligra.adj
#ligra_run BC 12 ~/graphs/twitter.adj
ligra_run BC 6662945 /graphs/twitter.adj.shuf
ligra_run BC 35 /graphs/protein.adj

#ligra_plus_run BC 9 ~/graphs/soc-LiveJournal_shuf.adj.compressed
#ligra_plus_run BC 28 ~/graphs/com-orkut.ungraph.adj.shuf.compressed
#ligra_plus_run BC 19372 ~/graphs/rmat_ligra.adj.shuf.compressed
#ligra_plus_run BC 0 ~/graphs/soc-LiveJournal1_sym.adj.compressed
#ligra_plus_run BC 1000 ~/graphs/com-orkut.ungraph.adj.compressed
#ligra_plus_run BC 0 ~/graphs/rmat_ligra.adj.compressed
##ligra_plus_run BC 12 ~/graphs/twitter.adj.compressed
##ligra_plus_run BC 6662945 ~/graphs/twitter.adj.shuf.compressed
#ligra_plus_run BC 35 ~/graphs/protein.adj.comressed

#aspen_run BC 9 /graphs/soc-LiveJournal_shuf.adj
#aspen_run BC 28 /graphs/com-orkut.ungraph.adj.shuf
#aspen_run BC 19372 /graphs/rmat_ligra.adj.shuf
#aspen_run BC 0 ~/graphs/soc-LiveJournal1_sym.adj
#aspen_run BC 1000 ~/graphs/com-orkut.ungraph.adj
#aspen_run BC 0 ~/graphs/rmat_ligra.adj
#aspen_run BC 12 ~/graphs/twitter.adj
#aspen_run BC 6662945 /graphs/twitter.adj.shuf
#aspen_run BC 35 /graphs/protein.adj

ligra_run Components 9 /graphs/soc-LiveJournal_shuf.adj
ligra_run Components 28 /graphs/com-orkut.ungraph.adj.shuf
ligra_run Components 19372 /graphs/rmat_ligra.adj.shuf
#ligra_run Components 0 ~/graphs/soc-LiveJournal1_sym.adj
#ligra_run Components 1000 ~/graphs/com-orkut.ungraph.adj
#ligra_run Components 0 ~/graphs/rmat_ligra.adj
#ligra_run Components 12 ~/graphs/twitter.adj
ligra_run Components 6662945 /graphs/twitter.adj.shuf
ligra_run Components 35 /graphs/protein.adj

#ligra_plus_run Components 9 ~/graphs/soc-LiveJournal_shuf.adj.compressed
#ligra_plus_run Components 28 ~/graphs/com-orkut.ungraph.adj.shuf.compressed
#ligra_plus_run Components 19372 ~/graphs/rmat_ligra.adj.shuf.compressed
#ligra_plus_run Components 0 ~/graphs/soc-LiveJournal1_sym.adj.compressed
#ligra_plus_run Components 1000 ~/graphs/com-orkut.ungraph.adj.compressed
#ligra_plus_run Components 0 ~/graphs/rmat_ligra.adj.compressed
##ligra_plus_run Components 12 ~/graphs/twitter.adj.compressed
##ligra_plus_run Components 6662945 ~/graphs/twitter.adj.shuf.compressed
#ligra_plus_run Components 35 ~/graphs/protein.adj.compressed

#aspen_run CC 9 /graphs/soc-LiveJournal_shuf.adj
#aspen_run CC 28 /graphs/com-orkut.ungraph.adj.shuf
#aspen_run CC 19372 /graphs/rmat_ligra.adj.shuf
#aspen_run CC 0 ~/graphs/soc-LiveJournal1_sym.adj
#aspen_run CC 1000 ~/graphs/com-orkut.ungraph.adj
#aspen_run CC 0 ~/graphs/rmat_ligra.adj
#aspen_run CC 12 ~/graphs/twitter.adj
#aspen_run CC 6662945 /graphs/twitter.adj.shuf
#aspen_run CC 35 /graphs/protein.adj

ligra_run_tc Triangle 9 /graphs/soc-LiveJournal_shuf.adj
ligra_run_tc Triangle 28 /graphs/com-orkut.ungraph.adj.shuf
ligra_run_tc Triangle 19372 /graphs/rmat_ligra.adj.shuf
#ligra_run_tc Triangle 0 ~/graphs/soc-LiveJournal1_sym.adj
#ligra_run_tc Triangle 1000 ~/graphs/com-orkut.ungraph.adj
#ligra_run Triangle 0 ~/graphs/rmat_ligra.adj

#aspen_run TC 9 ~/graphs/soc-LiveJournal_shuf.adj
#aspen_run TC 28 ~/graphs/com-orkut.ungraph.adj.shuf
#aspen_run TC 19372 ~/graphs/rmat_ligra.adj.shuf
#aspen_run TC 0 ~/graphs/soc-LiveJournal1_sym.adj
#aspen_run TC 1000 ~/graphs/com-orkut.ungraph.adj
#aspen_run TC 0 ~/graphs/rmat_ligra.adj

ligra_run PageRank 9 /graphs/soc-LiveJournal_shuf.adj
ligra_run PageRank 28 /graphs/com-orkut.ungraph.adj.shuf
ligra_run PageRank 19372 /graphs/rmat_ligra.adj.shuf
#ligra_run PageRank 0 ~/graphs/soc-LiveJournal1_sym.adj
#ligra_run PageRank 1000 ~/graphs/com-orkut.ungraph.adj
#ligra_run PageRank 0 ~/graphs/rmat_ligra.adj
#ligra_run PageRank 12 ~/graphs/twitter.adj
ligra_run PageRank 6662945 /graphs/twitter.adj.shuf
ligra_run PageRank 35 /graphs/protein.adj

#ligra_plus_run PageRank 9 ~/graphs/soc-LiveJournal_shuf.adj.compressed
#ligra_plus_run PageRank 28 ~/graphs/com-orkut.ungraph.adj.shuf.compressed
#ligra_plus_run PageRank 19372 ~/graphs/rmat_ligra.adj.shuf.compressed
#ligra_plus_run PageRank 0 ~/graphs/soc-LiveJournal1_sym.adj.compressed
#ligra_plus_run PageRank 1000 ~/graphs/com-orkut.ungraph.adj.compressed
#ligra_plus_run PageRank 0 ~/graphs/rmat_ligra.adj.compressed
##ligra_plus_run PageRank 12 ~/graphs/twitter.adj.compressed
##ligra_plus_run PageRank 6662945 ~/graphs/twitter.adj.shuf.compressed
#ligra_plus_run PageRank 35 ~/graphs/protein.adj.comressed

#aspen_run PR 9 /graphs/soc-LiveJournal_shuf.adj
#aspen_run PR 28 /graphs/com-orkut.ungraph.adj.shuf
#aspen_run PR 19372 /graphs/rmat_ligra.adj.shuf
#aspen_run PR 0 ~/graphs/soc-LiveJournal1_sym.adj
#aspen_run PR 1000 ~/graphs/com-orkut.ungraph.adj
#aspen_run PR 0 ~/graphs/rmat_ligra.adj
#aspen_run PR 12 ~/graphs/twitter.adj
#aspen_run PR 6662945 /graphs/twitter.adj.shuf
#aspen_run PR 35 /graphs/protein.adj

ligra_run BellmanFord 9 /graphs/rand-out-soc-LiveJournal1.txt.shuf.weighted.adj
ligra_run BellmanFord 28 /graphs/com-orkut.sym.txt.shuf.weighted.adj
ligra_run BellmanFord 19372 /graphs/rmat_ligra.el.shuf.weighted.adj
