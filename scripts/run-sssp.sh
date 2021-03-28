
# replace with your local paths
LJ_PATH=~/../wheatman/extended-csr/graphs/rand-out-soc-LiveJournal1.txt.shuf
ORKUT_PATH=/efs/home/hjxu/aspen/inputs/com-orkut.ungraph.txt.shuf
RMAT_PATH=~/../wheatman/extended-csr/rmat_ligra.el.shuf

# weighted paths
LJ_WGH_PATH=data/rand-out-soc-LiveJournal1.txt.shuf.weighted
ORKUT_WGH_PATH=data/com-orkut.ungraph.txt.shuf.weighted
RMAT_WGH_PATH=data/rmat_ligra.el.shuf.weighted

LJ_SRC=9
ORKUT_SRC=28
RMAT_SRC=19372

# SSSP BF
numactl -i all ./graph_bm_weighted -t SSSP_BF -src $LJ_SRC -f $LJ_WGH_PATH 
numactl -i all ./graph_bm_weighted -t SSSP_BF -src $ORKUT_SRC -f $ORKUT_WGH_PATH 
numactl -i all ./graph_bm_weighted -t SSSP_BF -src $RMAT_SRC -f $RMAT_WGH_PATH 

