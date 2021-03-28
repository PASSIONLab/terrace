#pragma once
#include "graph.h"
#include "VertexSubset.cpp"

using namespace graphstore;

template<class F>
struct EDGE_MAP_SPARSE {
  Graph &G;
  VertexSubset &output_vs;
  F f;
  bool output;
  EDGE_MAP_SPARSE(Graph &G_, VertexSubset &output_vs_, F f_, bool output_) : 
    G(G_), output_vs(output_vs_), f(f_), output(output_) {}
  inline bool update(uint32_t val) {
    G.map_sparse(f, output_vs, val, output);
    return false;
  }
};

template <class F>
VertexSubset EdgeMapSparse(Graph &G, VertexSubset &vs, F f, bool output) {
  // printf("edge map sparse\n");
  vs.convert_to_sparse();
  VertexSubset output_vs = VertexSubset(vs, false);
  struct EDGE_MAP_SPARSE<F> v(G, output_vs, f, output);
  vs.map(v);
  output_vs.finalize();
  return output_vs;
}

template <class F>
VertexSubset EdgeMapDense(Graph &G, VertexSubset &vs, F f, bool output) {
  // printf("edge map dense\n");
  vs.convert_to_dense();
  VertexSubset output_vs = VertexSubset(vs, false);
  // needs a grainsize of at least 512 
  // so writes to the bitvector storing the next vertex set are going to different cache lines
  if (vs.all) {
      parallel_for(uint64_t i_ = 0; i_ < G.get_num_vertices(); i_+=512) {
        uint64_t end = std::min(i_+512, (uint64_t) G.get_num_vertices());
        for (uint64_t i = i_; i < end; i++) {
          if (f.cond(i) == 1) {
            //printf("processing row %lu\n", i);
            G.map_dense_vs_all(f, vs, output_vs, i, output);
          }
        }
    }
  } else {
			parallel_for(uint64_t i_ = 0; i_ < G.get_num_vertices(); i_+=512) {
				uint64_t end = std::min(i_+512, (uint64_t) G.get_num_vertices());
				for (uint64_t i = i_; i < end; i++) {
			//for(uint64_t i = 0; i < G.get_num_vertices(); i++) {
          if (f.cond(i) == 1) {
            //printf("processing row %lu\n", i);
            G.map_dense_vs_not_all(f, vs, output_vs, i, output);
          }
				}
		}
  }
      return output_vs;

}

template <class F>
VertexSubset edgeMap(Graph &G, VertexSubset &vs, F f, bool output = true, uint32_t threshold = 20) {
  //vs.print();
  //printf("%u, %u, %u\n", G.rows, threshold, vs.get_n());
  if (G.get_num_vertices()/threshold <= vs.get_n()) {
    return EdgeMapDense(G, vs, f, output);
  } else {
    return EdgeMapSparse(G, vs, f, output);
  }
}

template<class F, bool output>
struct VERTEX_MAP {
  VertexSubset &output_vs;
  F f;
  VERTEX_MAP(VertexSubset &output_vs_, F f_) : output_vs(output_vs_), f(f_) {}
  inline bool update(uint32_t val) {
    if constexpr (output) {
      if (f(val) == 1) {
        output_vs.insert(val);
      }
    } else {
      f(val);
    }
    return false;
  }
};

template <class F>
inline VertexSubset vertexMap(VertexSubset &vs, F f, bool output = true) {
  //TODO the compilier should have been able to do this itself
  if (output) {
    VertexSubset output_vs = VertexSubset(vs, false);
    struct VERTEX_MAP<F, true> v(output_vs, f);
    vs.map(v);
    output_vs.finalize();
    return output_vs;
  } else {
    VertexSubset null_vs = VertexSubset();
    struct VERTEX_MAP<F, false> v(null_vs, f);
    vs.map(v);
    return null_vs;
  }
}
