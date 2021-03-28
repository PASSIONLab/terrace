
# replace with your local paths
LJ_PATH=~/graphs/soc-LiveJournal_shuf.adj
ORKUT_PATH=~/graphs/com-orkut.ungraph.adj.shuf
RMAT_PATH=~/graphs/rmat_ligra.adj.shuf
TWITTER_PATH=~/graphs/twitter.adj.shuf

LJ_SRC=9
ORKUT_SRC=28
RMAT_SRC=19372
TWITTER_SRC=6662945

# all
numactl -i all ./graph_bm  -src $LJ_SRC -f $LJ_PATH -rounds 10 
numactl -i all ./graph_bm  -src $ORKUT_SRC -f $ORKUT_PATH -rounds 10
numactl -i all ./graph_bm  -src $RMAT_SRC -f $RMAT_PATH  -rounds 10
numactl -i all ./graph_bm  -src $TWITTER_SRC -f $TWITTER_PATH  -rounds 10

# bfs
#numactl -i all ./graph_bm -t BFS -src $LJ_SRC -f $LJ_PATH -rounds 10 
#numactl -i all ./graph_bm -t BFS -src $ORKUT_SRC -f $ORKUT_PATH -rounds 10
#numactl -i all ./graph_bm -t BFS -src $RMAT_SRC -f $RMAT_PATH  -rounds 10
#numactl -i all ./graph_bm -t BFS -src $TWITTER_SRC -f $TWITTER_PATH  -rounds 10

# PR
#numactl -i all ./graph_bm -t PR -f $LJ_PATH -rounds 10
#numactl -i all ./graph_bm -t PR -f $ORKUT_PATH -rounds 10
#numactl -i all ./graph_bm -t PR -f $RMAT_PATH -rounds 10
#numactl -i all ./graph_bm -t PR -f $TWITTER_PATH -rounds 10

# BC
#numactl -i all ./graph_bm -t BC -src $LJ_SRC -f $LJ_PATH -rounds 10
#numactl -i all ./graph_bm -t BC -src $ORKUT_SRC -f $ORKUT_PATH -rounds 10
#numactl -i all ./graph_bm -t BC -src $RMAT_SRC -f $RMAT_PATH -rounds 10
#numactl -i all ./graph_bm -t BC -src $TWITTER_SRC -f $TWITTER_PATH -rounds 10

# CC
#numactl -i all ./graph_bm -t CC -f $LJ_PATH -rounds 10
#numactl -i all ./graph_bm -t CC -f $ORKUT_PATH -rounds 10
#numactl -i all ./graph_bm -t CC -f $RMAT_PATH -rounds 10
#numactl -i all ./graph_bm -t CC -f $TWITTER_PATH -rounds 10

# TC
#numactl -i all ./graph_bm -t TC -f $LJ_PATH -rounds 10
#numactl -i all ./graph_bm -t TC -f $ORKUT_PATH -rounds 10
