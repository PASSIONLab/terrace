# Terrace
Terrace: A Hierarchical Graph Container for Skewed Dynamic Graphs

This work appeared at SIGMOD 2021. If you use this software please cite us:

```
@inproceedings{PandeyWHB21,
  author    = {Prashant Pandey and
               Brian Wheatman and
               Helen Xu and
               Aydin Buluc},
  title     = {Terrace: A Hierarchical Graph Container for Skewed Dynamic Graphs},
  booktitle={Proceedings of the 2021 ACM international conference on Management of Data},
  year      = {2021},
}
```

Overview
--------
Terrace is a system for streaming graphs that uses
a hierarchical data structure design to store a vertex’s neighbors
in different data structures depending on the degree of the vertex.
This multi-level structure enables Terrace to dynamically partition
vertices based on their degrees and adapt to skewness in the
underlying graph. 

Build
-------
We implemented Terrace as a C++ library parallelized using Cilk Plus and the
Tapir branch of the LLVM compiler.

We recommend compiling Terrace with Clang++. We downloaded OpenCilk from
[here](https://github.com/OpenCilk/opencilk-project/releases). For recent
updates regarding OpenCilk please check [here](https://cilk.mit.edu/).

```bash
make graph_bm_no_update
./graph_bm_no_update -src 9 -maxiters 5 -f data/slashdot.adj 
```

Contributing
------------
Contributions via GitHub pull requests are welcome.


Authors
-------
- Prashant Pandey <ppandey@berkeley.edu>
- Brian Wheatman <wheatman@cs.jhu.edu>
- Helen Xu <hjxu@mit.edu>
- Aydin Buluc <abuluc@lbl.gov>
