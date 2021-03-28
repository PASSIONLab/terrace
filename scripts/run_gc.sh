#!/bin/bash
set -x
ROUNDS=10

GC_QUERIES=gc_queries
GC_UPDATES=gc_updates

# run single-core and not single-core
function gc_run_lj {
  # output file is in the format graph_numcores_srcvtx
#numactl -i all ./graphstore/graph_bm -TC 1 -f $1 -rounds $ROUNDS -src $2 >> $3\_96\_$2 
  taskset -c 0-23,48-71 ./graphstore/graph_bm_no_update -TC 1 -f $1 -rounds $ROUNDS -src $2 >> $3\_48\_$2 
  #taskset -c 0 ./graphstore/graph_bm_no_update -TC 1 -f $1 -rounds $ROUNDS -src $2 >> $3\_1\_$2 
  #taskset -c 0-23 ./graphstore/graph_bm -TC 1 -f $1 -rounds $ROUNDS -src $2 >> $3\_24\_$2 
  #taskset -c 0 ./graphstore/graph_bm_no_update -TC 1 -f $1 -rounds $ROUNDS -src $2 >> $3\_1\_$2 
#python post_gc.py $3\_96\_$2 $GC_QUERIES $GC_UPDATES
  python post_gc.py $3\_48\_$2 $GC_QUERIES $GC_UPDATES
  #python post_gc.py $3\_1\_$2 $GC_QUERIES $GC_UPDATES
  #python post_gc.py $3\_24\_$2 $GC_QUERIES $GC_UPDATES
  #python post_gc.py $3\_1\_$2 $GC_QUERIES $GC_UPDATES
}

function gc_run_updates_no_tc {
  # output file is in the format graph_numcores_srcvtx
#numactl -i all ./graphstore/graph_bm -f $1 -rounds $ROUNDS -src $2 >> $3\_96\_$2 
  #taskset -c 0-23,48-71 ./graphstore/graph_bm -f $1 -rounds $ROUNDS -src $2 >> $3\_48\_$2 
  taskset -c 0 ./graphstore/graph_bm -f $1 -rounds $ROUNDS -src $2 >> $3\_1\_$2 
  #taskset -c 0-23 ./graphstore/graph_bm -f $1 -rounds $ROUNDS -src $2 >> $3\_24\_$2 
  #taskset -c 0 ./graphstore/graph_bm_no_update -f $1 -rounds $ROUNDS -src $2 >> $3\_1\_$2 
#python post_gc.py $3\_96\_$2 $GC_QUERIES $GC_UPDATES
  #python post_gc.py $3\_48\_$2 $GC_QUERIES $GC_UPDATES
  python post_gc.py $3\_1\_$2 $GC_QUERIES $GC_UPDATES
  #python post_gc.py $3\_24\_$2 $GC_QUERIES $GC_UPDATES
  #python post_gc.py $3\_1\_$2 $GC_QUERIES $GC_UPDATES
}

function gc_run_no_updates_with_tc {
  # output file is in the format graph_numcores_srcvtx
#numactl -i all ./graphstore/graph_bm_no_update -TC 1 -f $1 -rounds $ROUNDS -src $2 >> $3\_96\_$2 
  taskset -c 0-23,48-71 ./graphstore/graph_bm_no_update -TC 1 -f $1 -rounds $ROUNDS -src $2 >> $3\_48\_$2 
  #taskset -c 0 ./graphstore/graph_bm_no_update -TC 1 -f $1 -rounds $ROUNDS -src $2 >> $3\_1\_$2 
  #taskset -c 0-23 ./graphstore/graph_bm_no_update -TC 1 -f $1 -rounds $ROUNDS -src $2 >> $3\_24\_$2 
  #taskset -c 0 ./graphstore/graph_bm_no_update -TC 1 -f $1 -rounds $ROUNDS -src $2 >> $3\_1\_$2 
#python post_gc.py $3\_96\_$2 $GC_QUERIES $GC_UPDATES
  python post_gc.py $3\_48\_$2 $GC_QUERIES $GC_UPDATES
  #python post_gc.py $3\_1\_$2 $GC_QUERIES $GC_UPDATES
  #python post_gc.py $3\_24\_$2 $GC_QUERIES $GC_UPDATES
  #python post_gc.py $3\_1\_$2 $GC_QUERIES $GC_UPDATES
}

function gc_run_no_updates_no_tc {
  # output file is in the format graph_numcores_srcvtx
#numactl -i all ./graphstore/graph_bm_no_update -f $1 -rounds $ROUNDS -src $2 >> $3\_96\_$2 
  taskset -c 0-23,48-71 ./graphstore/graph_bm_no_update -f $1 -rounds $ROUNDS -src $2 >> $3\_48\_$2 
  #taskset -c 0 ./graphstore/graph_bm_no_update -f $1 -rounds $ROUNDS -src $2 >> $3\_1\_$2 
  #taskset -c 0-23 ./graphstore/graph_bm_no_update -f $1 -rounds $ROUNDS -src $2 >> $3\_24\_$2 
  #taskset -c 0 ./graphstore/graph_bm_no_update -f $1 -rounds $ROUNDS -src $2 >> $3\_1\_$2 
#python post_gc.py $3\_96\_$2 $GC_QUERIES $GC_UPDATES
  python post_gc.py $3\_48\_$2 $GC_QUERIES $GC_UPDATES
  #python post_gc.py $3\_1\_$2 $GC_QUERIES $GC_UPDATES
  #python post_gc.py $3\_24\_$2 $GC_QUERIES $GC_UPDATES
  #python post_gc.py $3\_1\_$2 $GC_QUERIES $GC_UPDATES
}

function gc_run_weighted {
  # output file is in the format graph_numcores_srcvtx
#numactl -i all ./graphstore/graph_bm_weighted -f $1 -rounds $ROUNDS -src $2 >> $3\_96\_$2 
  #taskset -c 0-23,48-71 ./graphstore/graph_bm_weighted -f $1 -rounds $ROUNDS -src $2 >> $3\_48\_weighted\_$2 
  taskset -c 0 ./graphstore/graph_bm_weighted -f $1 -rounds $ROUNDS -src $2 >> $3\_1\_weighted\_$2 
  #taskset -c 0-23 ./graphstore/graph_bm_weighted -f $1 -rounds $ROUNDS -src $2 >> $3\_24\_$2 
  #taskset -c 0 ./graphstore/graph_bm_weighted -f $1 -rounds $ROUNDS -src $2 >> $3\_1\_$2 

  # post-process log into times
#python post_gc.py $3\_96\_$2 $GC_QUERIES $GC_UPDATES
  #python post_gc.py $3\_48\_weighted\_$2 $GC_QUERIES $GC_UPDATES
  python post_gc.py $3\_1\_weighted\_$2 $GC_QUERIES $GC_UPDATES
  #python post_gc.py $3\_24\_$2 $GC_QUERIES $GC_UPDATES
  #python post_gc.py $3\_1\_$2 $GC_QUERIES $GC_UPDATES
}

# set up output file
echo "structure,rounds,algorithm,src,graph,threads,time" >> $GC_QUERIES
echo "structure,update type,batch size,graph,throughput" >> $GC_UPDATES

: <<'END'
# test all subroutines
gc_run_lj ~/../wheatman/graphs/slashdot.adj 9 slashdot
gc_run_updates_no_tc ~/../wheatman/graphs/slashdot.adj 9 slashdot
gc_run_no_updates_with_tc ~/../wheatman/graphs/slashdot.adj 9 slashdot
gc_run_no_updates_no_tc ~/../wheatman/graphs/slashdot.adj 9 slashdot
END

#TODO: update paths according to local setup
# args are <path to graph> <src vtx> <graph name>
# outputs into file graph_numcores_srcvtx for cores 1 and 24


# yes queries, yes TC, yes update: LJ
gc_run_lj /graphs/soc-LiveJournal_shuf.adj 9 lj*
#gc_run_lj ~/graphs/soc-LiveJournal1_sym.adj 0 lj

# yes queries, yes TC, no update: orkut, rmat
gc_run_no_updates_with_tc /graphs/com-orkut.ungraph.adj.shuf 28 orkut*
gc_run_no_updates_with_tc /graphs/rmat_ligra.adj.shuf 19372 rmat*
gc_run_no_updates_no_tc /graphs/protein.adj 35 protein
gc_run_no_updates_no_tc /graphs/twitter.adj.shuf 6662945 twitter*

# yes queries, no TC, yes update: ER
# gc_run_no_updates_with_tc /graphs/com-orkut.ungraph.adj 1000 orkut
# gc_run_updates_no_tc /graphs/er_graph.adj 0 ER
# yes queries, no TC, no update: twitter, protein
# gc_run_no_updates_no_tc /graphs/twitter.adj 12 twitter
# gc_run_no_updates_with_tc /graphs/rmat_ligra.adj 0 rmat

'''
# sssp 
gc_run_weighted /graphs/rand-out-soc-LiveJournal1.txt.shuf.weighted.adj 9 LJ
gc_run_weighted /graphs/com-orkut.sym.txt.shuf.weighted.adj 28 Orkut
gc_run_weighted /graphs/rmat_ligra.el.shuf.weighted.adj 19372 rMAT
'''
