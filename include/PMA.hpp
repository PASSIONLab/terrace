#pragma once

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <vector>
#include <tuple>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <immintrin.h>
#include <atomic>


#include <stdint.h>
#include <queue>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "parallel.h"
#include "PMA_Lock.hpp"

#include "integerSort/blockRadixSort/blockRadixSort.h"

typedef struct _edge {
  uint32_t value; // length in array
  uint32_t dest;
} edge_t;

typedef struct _pair_double {
  double x;
  double y;
} pair_double;

#define EDGELONG 1
#if EDGELONG == 1
typedef uint64_t uint_t;
typedef int64_t int_t;
#define REDISTRIBUTE_PAR_SIZE (UINT64_MAX)
#else
typedef uint32_t uint_t;
typedef int32_t int_t;
#define REDISTRIBUTE_PAR_SIZE (UINT32_MAX)
#endif


typedef struct _node {
  // beginning and end of the associated region in the edge list
  uint_t beginning;     // deleted = max int
  uint_t end;           // end pointer is exclusive
  uint32_t num_neighbors; // number of edgess with this node as source
#if ENABLE_PMA_LOCK == 1
  PMA_Lock lock;
#endif
} node_t;


#define NULL_VAL (UINT32_MAX)
#define SENT_VAL (UINT32_MAX -1)
//#define REDISTRIBUTE_PAR_SIZE (4096)

typedef struct edge_list {
  uint_t N;
  uint32_t H;
  uint32_t logN;
  uint32_t loglogN;
  uint32_t mask_for_leaf;
  uint32_t * vals;
  uint32_t * dests;

  // Lock list_lock;

  double density_limit;
} edge_list_t;

class PMA {
public:
  // data members
  edge_list_t edges;
  std::vector<node_t> nodes;
#if ENABLE_PMA_LOCK == 1
  PMA_Lock node_lock;
#endif
  uint64_t next_task_id;

  double upper_density_bound[32];
  double lower_density_bound[32];
  // graph_t g;

  // function headings
  PMA(uint32_t init_n = 16);
  PMA(PMA &other);
  ~PMA();
  void double_list(uint64_t task_id, std::vector<uint_t> &sub_counts, uint_t num_elements);
  void half_list(uint64_t task_id, std::vector<uint_t> &sub_counts, uint_t num_elements);
  //void half_list();
  void slide_right(uint_t index, uint32_t *vals, uint32_t *dests);
  void slide_left(uint_t index, uint32_t *vals, uint32_t *dests);
  void redistribute(uint_t index, uint64_t len);
  void redistribute_par(uint_t index, uint64_t len, std::vector<uint_t> &sub_counts, uint_t num_elements, bool for_double = false);
  uint_t fix_sentinel(uint32_t node_index, uint_t in);
  void print_array(uint64_t worker_num = 0);
  uint32_t find_value(uint32_t src, uint32_t dest);
  void print_graph();
  void add_node();

  //assumes the edge is not already in the graph
  void add_edge(uint32_t src, uint32_t dest, uint32_t value);
  bool add_edge_update(uint32_t src, uint32_t dest, uint32_t value);


  bool add_edge_update_fast(uint32_t src, uint32_t dest, uint32_t value, uint64_t task_id);
  
  void add_edge_batch_update(uint32_t *srcs, uint32_t *dests, uint32_t *values, uint_t edge_count);
  void add_edge_batch_update_no_val(uint32_t *srcs, uint32_t *dests, uint_t edge_count);

  // merge functions in original PMA with no val
  void add_edge_batch_update_no_val_parallel(pair_uint *es, uint64_t edge_count);
  void add_edge_batch_wrapper(pair_uint *es, uint64_t edge_count, int64_t threshold = -1);
  void remove_edge_batch_update_no_val_parallel(pair_uint *es, uint64_t edge_count);
  void remove_edge_batch_wrapper(pair_uint *es, uint64_t edge_count, int64_t threshold = -1);

  bool remove_edge(uint32_t src, uint32_t dest);
  void remove_edge_fast(uint32_t src, uint32_t dest, uint64_t task_id);
  void remove_edge_batch(uint32_t *srcs, uint32_t *dests, uint_t edge_count);


  void insert(uint64_t task_id, uint_t index, uint32_t elem_dest, uint32_t elem_value, uint32_t src
  	#if ENABLE_PMA_LOCK == 1
  	, pair_int held_locks
  	#endif
  	);
  void remove(uint64_t task_id, uint_t index, uint32_t elem_dest, uint32_t src
    #if ENABLE_PMA_LOCK == 1
    , pair_int held_locks
    #endif
    );
  uint64_t get_size();
  uint64_t get_n();
  std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> get_edges();
  void clear();
  bool check_no_locks();
  bool check_no_locks_for_me(uint64_t task_id);
  bool check_no_node_locks_for_me(uint64_t task_id);
  bool check_no_node_locks_held_by_me(uint64_t task_id);
  bool grab_all_locks(uint64_t task_id, bool exclusive, REASONS reason = GENERAL);
  void release_all_locks(uint64_t task_id, bool exclusive, REASONS reason = GENERAL);
  pair_int grab_locks_in_range_with_resets(uint64_t task_id, uint_t index, uint_t len, REASONS reason, uint32_t guess);
  pair_int grab_locks_for_leaf_with_resets(uint64_t task_id, uint32_t src, REASONS reason = GENERAL);
  void release_locks_in_range(uint64_t task_id, pair_int locks, REASONS reason = GENERAL);
  uint32_t find_contaning_node(uint_t index);

  pair_int which_locks_in_range(uint_t index, uint_t len, uint32_t guess);
  pair_int which_locks_for_leaf(uint32_t src);
  bool check_every_lock_in_leaf(uint64_t task_id, uint64_t index);
  bool check_every_lock_in_node(uint64_t task_id, uint64_t index, uint64_t len);
  
  uint32_t num_neighbors(uint32_t node) {
    return nodes[node].num_neighbors;
  }
  uint64_t num_edges() {
    uint64_t num = 0;
    for (uint64_t i = 0; i < get_n(); i++) {
      num += num_neighbors(i);
    }
    return num;
  }

  class iterator {
  public:
    uint_t place;
    uint_t end;
    uint32_t * vals;
    uint32_t * dests;
    uint8_t loglogN;
    iterator(const PMA *G, uint32_t node, bool start) {
      if (!start) {
        place = G->nodes[node].end;
        return;
      }
      place = G->nodes[node].beginning + 1;
      end = G->nodes[node].end;
      vals = (uint32_t *)G->edges.vals;
      dests = (uint32_t *)G->edges.dests;
      loglogN = G->edges.loglogN;
      while ((place < end) && (dests[place] == NULL_VAL)) {
        place = ((place >> loglogN) + 1) << (loglogN);
      }
      if (place > end) {
        place = end;
      }
      return;
    }
		iterator () {}
		iterator(const iterator & other) {
			place = other.place;
			end = other.end;
			vals = other.vals;
			dests = other.dests;
			loglogN = other.loglogN;
		}
    bool operator==(const iterator& other) const {
      return (place == other.place);
    }
    bool operator!=(const iterator& other) const {
      return (place != other.place);
    }
    iterator& operator++() {
      place += 1;
      while ((place < end) && (dests[place] == NULL_VAL)) {
        place = ((place >> loglogN) + 1) << (loglogN);
      }
      if (place > end) {
        place = end;
      }
      return *this;
    }
    edge_t operator*() const {
      return {vals[place], dests[place]};
    }
  };
  iterator begin(uint32_t node) const {
    return iterator(this, node, true);
  }
  iterator end(uint32_t node) const {
    return iterator(this, node, false);
  }
};

// same as find_leaf, but does it for any level in the tree
// index: index in array
// len: length of sub-level.
static inline uint_t find_node(uint_t index, uint_t len) { return (index / len) * len; }
