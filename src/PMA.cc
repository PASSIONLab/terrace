using namespace std;

#include "PMA.hpp"
#include <immintrin.h>

# define assert(x)

// given index, return the starting index of the leaf it is in
//TODO this could be aster if we store a mask and just do a single and
static uint_t find_leaf(edge_list_t *list, uint_t index) {
  return index & list->mask_for_leaf;
}


static uint_t find_prev_valid(uint32_t volatile  * volatile dests, uint_t start) {
  while (dests[start] == NULL_VAL) {
    start--;
  }
  return start;
}

/*
// same as find_leaf, but does it for any level in the tree
// index: index in array
// len: length of sub-level.
static inline int find_node(int index, int len) { return (index / len) * len; }
*/

bool PMA::check_no_locks() {
  bool ret = true;
#if ENABLE_PMA_LOCK == 1
  ret = ret && node_lock.check_unlocked();
  assert(node_lock.check_unlocked());
  // ret = ret && edges.list_lock.check_unlocked();
  // assert(edges.list_lock.check_unlocked());
  for (uint32_t i = 0; i < nodes.size(); i++) {
    ret = ret && nodes[i].lock.check_unlocked();
    if (!ret) {
      printf("found lock on iter %d\n", i);
    }
    ret = ret && (nodes[i].lock.reason == GENERAL);
    if (!ret) {
      printf("found lock on iter %d with %d reason\n", i, nodes[i].lock.reason);
    }
    assert(nodes[i].lock.check_unlocked());
    assert(nodes[i].lock.reason == GENERAL);
  }
#endif
  return ret;
}

bool PMA::check_no_locks_for_me(uint64_t task_id) {

  bool ret = true;
#if ENABLE_PMA_LOCK == 1
  ret = ret && !node_lock.i_own_lock(task_id);
  assert(!node_lock.i_own_lock(task_id));
  // ret = ret && edges.list_lock.check_unlocked();
  // assert(edges.list_lock.check_unlocked());
  for (uint32_t i = 0; i < nodes.size(); i++) {
    ret = ret && !nodes[i].lock.i_own_lock(task_id);
    if (!ret) {
      printf("found lock on iter %d\n", i);
    }
    ret = ret && (nodes[i].lock.reason_set_by != task_id);
    if (!ret) {
      printf("found lock on iter %d that worker %u set the reason with %d\n", i,task_id, nodes[i].lock.reason);
    }
    assert(!nodes[i].lock.i_own_lock(task_id));
    assert(nodes[i].lock.reason_set_by != task_id);
  }
#endif
  return ret;
}  

bool PMA::check_no_node_locks_for_me(uint64_t task_id) {
#if ENABLE_PMA_LOCK == 0
  return true;
#else
bool ret = true;
  for (uint32_t i = 0; i < nodes.size(); i++) {
    ret = ret && !nodes[i].lock.i_own_lock(task_id);
    if (!ret) {
      printf("found lock on iter %d\n", i);
    }
    ret = ret && (nodes[i].lock.reason_set_by != task_id);
    if (!ret) {
      printf("found lock on iter %d that worker %u set the reason with %d\n", i,task_id, nodes[i].lock.reason);
    }
    assert(!nodes[i].lock.i_own_lock(task_id));
    assert(nodes[i].lock.reason_set_by != task_id);
  }
  return ret;
#endif
}  

bool PMA::check_no_node_locks_held_by_me(uint64_t task_id) {
  #if ENABLE_PMA_LOCK != 1
      return true;
  #else
  bool ret = true;
  for (uint32_t i = 0; i < nodes.size(); i++) {
    ret = ret && !nodes[i].lock.i_own_lock(task_id);
    if (!ret) {
      printf("found lock on iter %d\n", i);
    }
    assert(!nodes[i].lock.i_own_lock(task_id));
  }
  return ret;
#endif
}  


uint_t next_leaf(uint_t index, int loglogN) {
  return ((index >> loglogN) + 1) << (loglogN);
}


bool PMA::grab_all_locks(uint64_t task_id, bool exclusive, REASONS reason) {
  #if ENABLE_PMA_LOCK != 1
      return true;
  #else
  for (uint32_t i = 0; i < nodes.size(); i++) {
    if (exclusive) {
      if (!nodes[i].lock.lock(task_id, reason)) {
        return false;
      }
    } else  {
      if (!nodes[i].lock.lock_shared(task_id)) {
        return false;
      }
    }
  }
  return true;
#endif
}
void PMA::release_all_locks(uint64_t task_id, bool exclusive, REASONS reason) {
  #if ENABLE_PMA_LOCK != 1
      return;
  #else
  parallel_for (uint32_t i = 0; i < nodes.size(); i++) {
    if (exclusive) {
      nodes[i].lock.unlock(task_id, reason);
    } else  {
      nodes[i].lock.unlock_shared(task_id);
    }
  }
#endif
}

void PMA::clear() {
  printf("clear called\n");
  uint_t n = 0;
  uint64_t task_id = __sync_fetch_and_add(&next_task_id, 2);
  grab_all_locks(task_id, true, GENERAL);
  
  free((void*)edges.vals);
  free((void*)edges.dests);
  edges.N = 2UL << bsr_word(n);
  // printf("%d\n", bsf_word(list->N));
  edges.loglogN = bsr_word(bsr_word(edges.N) + 1);
  edges.logN = (1 << edges.loglogN);
  edges.H = bsr_word(edges.N / edges.logN);
}


// TODO jump to next leaf
vector<tuple<uint32_t, uint32_t, uint32_t>> PMA::get_edges() {
  // TODO grab locks in the lock list
  // for now, grabs the global lock
#if ENABLE_PMA_LOCK == 1
  node_lock.lock_shared(); // lock node array
#endif
  // edges.list_lock.lock();
  uint64_t n = get_n();
  vector<tuple<uint32_t, uint32_t, uint32_t>> output;

  for (uint_t i = 0; i < n; i++) {
    uint_t start = nodes[i].beginning;
    uint_t end = nodes[i].end;
#if ENABLE_PMA_LOCK == 1
    nodes[i].lock.lock_shared();
#endif
    for (uint_t j = start + 1; j < end; j++) {
      if (edges.dests[j]!=NULL_VAL) {
        output.push_back(
            make_tuple(i, edges.dests[j], edges.vals[j]));
      }
    }
#if ENABLE_PMA_LOCK == 1
    nodes[i].lock.unlock_shared();
#endif
    }
  // edges.list_lock.unlock();
#if ENABLE_PMA_LOCK == 1
  node_lock.unlock_shared(); // lock node array
#endif
  return output;
}


uint64_t PMA::get_n() {
#if ENABLE_PMA_LOCK == 1
  node_lock.lock_shared();
#endif
  uint64_t size = nodes.size();
#if ENABLE_PMA_LOCK == 1
  node_lock.unlock_shared();
#endif
  return size;
}

uint64_t PMA::get_size() {
#if ENABLE_PMA_LOCK == 1
  node_lock.lock_shared();
#endif
  uint64_t size = nodes.capacity() * sizeof(node_t);
  size += sizeof(PMA);
  size += (uint64_t)edges.N * sizeof(edge_t);
#if ENABLE_PMA_LOCK == 1
  node_lock.unlock_shared();
#endif
  return size;
}


void print_array(edge_list_t *edges) {
  printf("N = %d, logN = %d\n", edges->N, edges->logN);
  for (uint_t i = 0; i < edges->N; i++) {
    if (edges->dests[i]==NULL_VAL) {
      printf("%d-x ", i);
    } else if ((edges->dests[i]==SENT_VAL) || i == 0) {
      uint32_t value = edges->vals[i];
      if (value == NULL_VAL) {
        value = 0;
      }
      printf("\n%d-s(%u):(?, ?) ", i, value);
    } else {
      printf("%d-(%d, %u) ", i, edges->dests[i], edges->vals[i]);
    }
  }
  printf("\n\n");
}

void PMA::print_array(uint64_t worker_num) {
  printf("worker num: %lu, N = %d, logN = %d, density_limit = %f\n", worker_num, edges.N, edges.logN, edges.density_limit);
  for (uint_t i = 0; i < edges.N; i++) {
    if (edges.dests[i]==NULL_VAL) {
      printf("%d-x ", i);
    } else if ((edges.dests[i] == SENT_VAL) || i == 0) {
      uint32_t value = edges.vals[i];
      if (value == NULL_VAL) {
        value = 0;
      }
      printf("\n worker num: %lu, %d-s(%u):(%d, %d)(%d) ", worker_num, i, value, nodes[value].beginning,
             nodes[value].end, nodes[value].num_neighbors);
      #ifndef NDEBUG
      #if ENABLE_PMA_LOCK == 1
        printf("(%s, %d by %u)", nodes[value].lock.check_unlocked() ? "Unlocked": "Locked", nodes[value].lock.reason, nodes[value].lock.owner);
      #endif
      #else
        printf(")");
      #endif
    } else {
      printf("%d-(%d, %u) ", i, edges.dests[i], edges.vals[i]);
    }
  }
  printf("\n\n");
}


uint_t get_density_count(edge_list_t *list, uint_t index, uint_t len) {
 /* 
  uint32_t full = 0;
  uint32_t i = index;
  while (i < index + len) {
    if (!is_null(list->items[i].e)) {
      full++;
      i++;
    } else {
      i = next_leaf(i, list->logN);
    }
  }
  return full;
*/
  /*
  cilk::reducer< cilk::op_add<uint32_t> > full;
  parallel_for (uint32_t i = index; i < index+len; i++) {
    if (!is_null(list->items[i].e)) {
      (*full)++;
    }
  }
  return full.get_value();
  */
  // fater without paralleliszation since it gets properly vectorized
  uint32_t * dests = (uint32_t *) list->dests;
  uint_t full = 0;
#ifdef  __AVX__
  if (len >= 8) {
    uint_t null_count = 0;
    for (uint_t i = index; i < index+len; i+=8) {
      //TODO if we keep things aligned this could be faster, but then we cant use realloc in double
      __m256i a = _mm256_loadu_si256((__m256i *)& dests[i]);
      __m256i b =  _mm256_set1_epi32(NULL_VAL);
      uint32_t add = __builtin_popcount(_mm256_movemask_ps((__m256)_mm256_cmpeq_epi32(a, b)));
      null_count += add;
    }
    full = len - null_count;
    return full;
  }
#endif
  for (uint_t i = index; i < index+len; i+=4) {
#ifdef __SSE2__
      __m128i a = _mm_load_si128((__m128i *)& dests[i]);
      __m128i b =  _mm_set1_epi32(NULL_VAL);
      uint32_t add = 4-__builtin_popcount(_mm_movemask_ps((__m128) _mm_cmpeq_epi32(a, b)));
#else
      uint32_t add = (dests[i]!=NULL_VAL) + (dests[i+1]!=NULL_VAL) + (dests[i+2]!=NULL_VAL) + (dests[i+3]!=NULL_VAL);
#endif
      //__sync_fetch_and_add(&full, add);
      full+=add;
  }
  return full;
  /*
  cilk::reducer< cilk::op_add<uint32_t> > full;
  parallel_for (uint32_t i = index; i < index+len; i+=4) {
      uint32_t add = !is_null(list->items[i].e) + !is_null(list->items[i+1].e) + !is_null(list->items[i+2].e) + !is_null(list->items[i+3].e);
      (*full)+=add;
  }
  return full.get_value();
  */
  
}

uint64_t get_density_count_par(edge_list_t *list, uint_t index, uint_t len, std::vector<uint_t> &sub_counts) {
  std::vector<uint64_t> worker_counts(getWorkers()*8);
  uint32_t volatile  * volatile dests = list->dests;
  parallel_for(uint_t j = index; j < index+len; j+= REDISTRIBUTE_PAR_SIZE) {
    uint_t full = 0;
    for (uint_t i = j; i < j+REDISTRIBUTE_PAR_SIZE; i+=4) {
        uint32_t add = (dests[i]!=NULL_VAL) + (dests[i+1]!=NULL_VAL) + (dests[i+2]!=NULL_VAL) + (dests[i+3]!=NULL_VAL);
        full+=add;
    }
    worker_counts[getWorkerNum()]+=full;
    sub_counts[(j-index)/REDISTRIBUTE_PAR_SIZE] = full;
  }

  uint64_t total = 0;  
  for (uint64_t i = 0; i < worker_counts.size(); i+=8) {
    total += worker_counts[i];
  }
  return total;  
}

// get density of a node
// should already be locked if you are calling get density
double get_density(edge_list_t *list, uint_t index, uint_t len) {
  double full_d = (double)get_density_count(list, index, len);
  return full_d / len;
}

bool check_no_full_leaves(edge_list_t *list, uint_t index, uint_t len) {
  for (uint_t i = index; i < index + len; i+= list->logN) {
    bool full = true;
    for (uint_t j = i; j < i + list->logN; j++) {
       if (list->dests[j]==NULL_VAL) {
        full = false;
      }
    }
    if (full) {
      return false;
    }
  }
  return true;
}

// height of this node in the tree
int get_depth(edge_list_t *list, uint_t len) { return bsr_word(list->N / len); }

// when adjusting the list size, make sure you're still in the
// density bound
pair_double density_bound(edge_list_t *list, int depth) {
  pair_double pair;

  // between 1/4 and 1/2
  // pair.x = 1.0/2.0 - (( .25*depth)/list->H);
  // between 1/8 and 1/4
  pair.x = 1.0 / 4.0 - ((.125 * depth) / list->H);
  pair.y = 3.0 / 4.0 + ((.25 * depth) / list->H);
  if (pair.y > list->density_limit) {
    pair.y = list->density_limit-.001;
  }
  return pair;
}

//TODO make it so the first element is known to always be the first element and don't special case it so much
//assumes the node_lock is held
//returns where me have to start checking for sentinals again
// ths is just the start + the degree so we can run some fasts paths
uint_t PMA::fix_sentinel(uint32_t node_index, uint_t in) {
  // we know the first sentinal will never move so we just ignore it
  assert(node_index > 0);
  //node_lock.lock_shared();
  nodes[node_index - 1].end = in;

  nodes[node_index].beginning = in;
  if (node_index == nodes.size() - 1) {
    nodes[node_index].end = edges.N - 1;
  }
  //node_lock.unlock_shared();
  return nodes[node_index].beginning + nodes[node_index].num_neighbors;
}


// Evenly redistribute elements in the ofm, given a range to look into
// index: starting position in ofm structure
// len: area to redistribute
// should already be locked
void PMA::redistribute(uint_t index, uint64_t len) {
  //printf("len = %u\n", len);
  assert(find_leaf(&edges, index) == index);
  
  //printf("REDISTRIBUTE START: index:%u, len %u, worker %lu\n", index, len, get_worker_num());
  //print_array(get_worker_num());
  // std::vector<edge_t> space(len); //
  //TODO if its small use the stack
  // for small cases put on the stack
  uint32_t *space_vals;
  uint32_t *space_dests;
  uint32_t * vals = (uint32_t *) edges.vals;
  uint32_t  * dests = (uint32_t *) edges.dests;
  uint_t j = 0;
  if (len == edges.logN) {
    return;
  } else {
    space_vals = (uint32_t *)malloc(len * sizeof(*(edges.vals)));
    if (space_vals == 0) {

      printf("bad malloc, len = %lu\n", len);
      while (1){}
    }
    space_dests = (uint32_t *)malloc(len * sizeof(*(edges.dests)));
  }

  // move all items in ofm in the range into
  // a temp array
  /*
  int i = index;
  while (i < index + len) {
    if (!is_null(edges.items[i])) {
      space[j] = edges.items[i];
      edges.items[i].value = 0;
      edges.items[i].dest = 0;
      i++;
      j++;
    } else {
      i = next_leaf(i, edges.logN);
    }
  }
  */
  //TODO could parralize if get_density_count gave us more data, but doesn't seem to be a bottle neck
  // could get better cache behavior if I go back and forth with reading and writing
  assert(len >= 8);
  // AVX code seems to have a bug, but I can't find it s leaving it in the iffalse
#ifdef false & __AVX__
  for (uint_t i = index; i < index + len; i+=8) {
    __m256i dests_vec = _mm256_loadu_si256((__m256i *)& dests[i]);
    __m256i vals_vec = _mm256_loadu_si256((__m256i *)& vals[i]);
    _mm256_storeu_si256((__m256i *) &space_dests[j], dests_vec);
    _mm256_storeu_si256((__m256i *) &space_vals[j], vals_vec);
    __m256i b =  _mm256_set1_epi32(NULL_VAL);
    j += 8-__builtin_popcount(_mm256_movemask_ps((__m256)_mm256_cmpeq_epi32(dests_vec, b)));
    _mm256_storeu_si256((__m256i *)& dests[i], b);
    _mm256_storeu_si256((__m256i *)& vals[i], _mm256_setzero_si256());
  }
#else
  for (uint_t i = index; i < index + len; i+=8) {
    for (uint_t k = i; k < i+8; k++) {
      space_vals[j] = vals[k];
      space_dests[j] = dests[k];
      // counting non-null edges
      j += (space_dests[j]!=NULL_VAL);
      // setting section to null
      //vals[i] = 0;
      //dests[i] = NULL_VAL;
    }
    memset (__builtin_assume_aligned((void*)&vals[i], 32), 0, 8*sizeof(uint32_t));
    //setting by byte, but NULL_VAL is all ones so it is fine
    memset (__builtin_assume_aligned((void*)&dests[i], 32), NULL_VAL, 8*sizeof(uint32_t));
  }
#endif
  
  /*
  if (((double)j)/len > ((double)(edges.logN-1)/edges.logN)) {
    printf("too dense in redistribute, j = %u, len = %u, index = %u for worker %lu\n",j, len, index, get_worker_num() );
    print_array(get_worker_num());
  }*/
  assert( ((double)j)/len <= ((double)(edges.logN-1)/edges.logN));

  uint_t num_leaves = len >> edges.loglogN;
  uint_t count_per_leaf = j / num_leaves;
  uint_t extra = j % num_leaves;
  __builtin_prefetch ((void *)&nodes, 0, 3);

  // parallizing does not make it faster
  uint_t end_sentinel = 0;
  for (uint_t i = 0; i < num_leaves; i++) {
    uint_t count_for_leaf = count_per_leaf + (i < extra);
    uint_t in = index + (edges.logN * (i));
    uint_t j2 = count_per_leaf*i + min(i,extra);
    //TODO could be parallized, but normally only up to size 32
    uint_t j3 = j2;
    memcpy(__builtin_assume_aligned((void*)&vals[in],16), (void*)&space_vals[j2], count_for_leaf*sizeof(uint32_t));
    /*
    for (uint32_t k = in; k < count_for_leaf+in; k++) {
      vals[k] = space_vals[j2];
      j2++;
    }
    */
    if (end_sentinel > in + count_for_leaf) {
      memcpy(__builtin_assume_aligned((void*)&dests[in], 16), (void*)&space_dests[j2], count_for_leaf*sizeof(uint32_t));
    } else { 
      for (uint_t k = in; k < count_for_leaf+in; k++) {
        dests[k] = space_dests[j3];
        if (dests[k]==SENT_VAL) {
          // fixing pointer of node that goes to this sentinel
          uint32_t node_index = vals[k];
          end_sentinel = fix_sentinel(node_index, k);
        }
        j3++;
      }
  }

  }
  free(space_dests);
  free(space_vals);
  /*
  if (!check_no_full_leaves(&edges, index, len)) {
    printf("some leaves are full, current density is %f, index = %u, len = %u\n", get_density(&edges, index, len), index, len);
    print_array(get_worker_num());
  }*/
  assert(check_no_full_leaves(&edges, index, len));
  //printf("REDISTRIBUTE END: index:%u, len %u, worker %lu\n", index, len, get_worker_num());
}

void PMA::redistribute_par(uint_t index, uint64_t len, std::vector<uint_t> &sub_counts, uint_t num_elements, bool for_double) {
  assert(find_leaf(&edges, index) == index);
  //printf("par len = %u\n", len);

  uint32_t *space_vals = (uint32_t *)aligned_alloc(64, len * sizeof(*(edges.vals)));
  uint32_t *space_dests = (uint32_t *)aligned_alloc(64, len * sizeof(*(edges.dests)));
  uint_t j = 0;

  // move all items in ofm in the range into
  // a temp array


  //TODO parallel prefix sum
  for (size_t i = 1; i < sub_counts.size(); i++) {
    sub_counts[i] += sub_counts[i-1];
    //printf("sub_counts[%d] = %u\n", i, sub_counts[i]);
  }

  // could get better cache behavior if I go back and forth with reading and writing
  uint32_t volatile  * volatile vals = edges.vals;
  uint32_t volatile  * volatile dests = edges.dests;
 /* 
  for (uint32_t i = index; i < index + len; i++) {
    j += (dests[i]!=NULL_VAL);
    printf("start = %u, i = %u\n", j, i);
  }
  */
  j = num_elements;
  uint_t end = index+len;
  if (for_double) {
    end = end/2;
  }
  parallel_for(uint_t k = index; k < end; k+=REDISTRIBUTE_PAR_SIZE) {
    //printf("index = %u, len = %u, k = %u REDISTRIBUTE_PAR_SIZE = %u\n", index, len, k , REDISTRIBUTE_PAR_SIZE);
      uint_t start;
      if (k == index) {
        start = 0;
      } else {
        start = sub_counts[(k-index)/REDISTRIBUTE_PAR_SIZE-1];
      }

      for (uint_t i = k; i < k+REDISTRIBUTE_PAR_SIZE; i+=8) {
        //printf("start = %u, i = %u\n", start, i);
        for (uint_t j = i; j < i+8; j++) {
          if (dests[j] != NULL_VAL) {
            space_vals[start] = vals[j];
            space_dests[start] = dests[j];
            start += 1;
          }
          //TODO use avx below when possible
          vals[i] = 0;
          dests[i] = NULL_VAL;
        }
        // setting section to null
        //_mm256_storeu_si256((__m256i *) (&edges.vals[i]), _mm256_setzero_si256());
        //_mm256_storeu_si256((__m256i *) (&edges.dests[i]), _mm256_set_epi32(NULL_VAL, NULL_VAL, NULL_VAL, NULL_VAL, NULL_VAL, NULL_VAL, NULL_VAL, NULL_VAL));  
    }

  }
/*  
  j = 0;
  for(uint32_t k = index; k < index+len; k+=REDISTRIBUTE_PAR_SIZE) {
    for (uint32_t i = k; i < k + REDISTRIBUTE_PAR_SIZE; i++) {
      space_vals[j] = vals[i];
      space_dests[j] = dests[i];
      // counting non-null edges
      j += (space_dests[j]!=NULL_VAL);
      // setting section to null
      vals[i] = 0;
      dests[i] = NULL_VAL;
    }
  }
*/
  assert( ((double)j)/len <= ((double)(edges.logN-1)/edges.logN));

  uint_t num_leaves = len >> edges.loglogN;
  uint_t count_per_leaf = j / num_leaves;
  uint_t extra = j % num_leaves;
  __builtin_prefetch ((void *)&nodes, 0, 3);
  // parallizing does not make it faster
  parallel_for (uint_t i = 0; i < num_leaves; i++) {
    uint_t count_for_leaf = count_per_leaf + (i < extra);
    uint_t in = index + ((i) << edges.loglogN);
    uint_t j2 = count_per_leaf*i + min(i,extra);
    //TODO could be parallized, but normally only up to size 32
    uint_t j3 = j2;
    assert(j3 < len);
    
    for (uint_t k = in; k < count_for_leaf+in; k++) {
      vals[k] = space_vals[j2];
      j2++;
    }
    for (uint_t k = in; k < count_for_leaf+in; k++) {
      dests[k] = space_dests[j3];
      if (dests[k]==SENT_VAL) {
        // fixing pointer of node that goes to this sentinel
        uint32_t node_index = vals[k];
         
        fix_sentinel(node_index, k);
      }
      j3++;
      assert(j3 < len);
    }

  }
  free(space_dests);
  free(space_vals);

  assert(check_no_full_leaves(&edges, index, len));
}


//TODO pass in subcounts and do redistibute_par when big
void PMA::double_list(uint64_t task_id, std::vector<uint_t> &sub_counts, uint_t num_elements) {
  //printf("doubling list by worker %lu\n", get_worker_num());
  grab_all_locks(task_id, true, DOUBLE);
  uint64_t new_N = edges.N * 2;
  edges.loglogN = bsr_word(bsr_word(edges.N) + 1);
  edges.logN = (1 << edges.loglogN);
  edges.mask_for_leaf = ~(edges.logN - 1);
  assert(edges.logN > 0);
  edges.density_limit = ((double) edges.logN - 1)/edges.logN;
  edges.H = bsr_word(new_N / edges.logN);
  for (uint32_t i = 0; i <= edges.H; i++) {
    upper_density_bound[i] = density_bound(&edges, i).y;
    lower_density_bound[i] = density_bound(&edges, i).x;
  }

  edges.N = new_N;


  /* 
  uint32_t *new_dests = (uint32_t *)malloc(new_N * sizeof(*(edges.dests)));
  uint32_t *new_vals = (uint32_t *)malloc(new_N * sizeof(*(edges.vals)));
  parallel_for (uint32_t i = 0; i < edges.N / 2; i++) {
    new_vals[i] = vals[i]; // setting second half to null
    new_dests[i] = dests[i]; // setting second half to null
  }
  edges.dests = new_dests;
  edges.vals = new_vals;
  vals = edges.vals;
  dests = edges.dests;
  */
  //edges.dests = (uint32_t *)realloc((void*)edges.dests, new_N * sizeof(*(edges.dests)));
  //edges.vals = (uint32_t *)realloc((void*)edges.vals, new_N * sizeof(*(edges.vals)));
  uint32_t *new_dests = (uint32_t *)aligned_alloc(32, new_N * sizeof(*(edges.dests)));
  uint32_t *new_vals = (uint32_t *)aligned_alloc(32, new_N * sizeof(*(edges.vals)));
  uint32_t * vals = (uint32_t *)edges.vals;
  uint32_t * dests = (uint32_t *)edges.dests;
  parallel_for (uint_t i = 0; i < new_N / 2; i++) {
    new_vals[i] = vals[i];
    new_dests[i] = dests[i];
  }
  parallel_for (uint_t i = new_N / 2; i < new_N; i++) {
    new_vals[i] = 0; // setting second half to null
    new_dests[i] = NULL_VAL; // setting second half to null
  }
  free((void*)edges.vals);
  edges.vals = new_vals;
  free((void*)edges.dests);
  edges.dests = new_dests;

  assert(edges.dests != NULL && edges.vals != NULL);
  //memset((void*)(edges.items+(new_N / 2)), 0, sizeof(edge_t) * (new_N / 2) );
  
  
  //printf("List doubled: N - %u, logN = %u, H = %u\n", edges.N, edges.logN, edges.H);
  if (sub_counts.size() == 0) {
    redistribute(0, edges.N);
  } else {
    redistribute_par(0, edges.N, sub_counts, num_elements, true);
  }
  release_all_locks(task_id, true, GENERAL);

}

void PMA::half_list(uint64_t task_id, std::vector<uint_t> &sub_counts, uint_t num_elements) {
  //printf("doubling list by worker %lu\n", get_worker_num());
  grab_all_locks(task_id, true, DOUBLE);
  uint64_t new_N = edges.N / 2;
  edges.loglogN = bsr_word(bsr_word(edges.N) - 1);
  edges.logN = (1 << edges.loglogN);
  edges.mask_for_leaf = ~(edges.logN - 1);
  assert(edges.logN > 0);
  edges.density_limit = ((double) edges.logN - 1)/edges.logN;
  edges.H = bsr_word(new_N / edges.logN);
  for (uint32_t i = 0; i <= edges.H; i++) {
    upper_density_bound[i] = density_bound(&edges, i).y;
    lower_density_bound[i] = density_bound(&edges, i).x;
  }


  /* 
  uint32_t *new_dests = (uint32_t *)malloc(new_N * sizeof(*(edges.dests)));
  uint32_t *new_vals = (uint32_t *)malloc(new_N * sizeof(*(edges.vals)));
  cilk_for (uint32_t i = 0; i < edges.N / 2; i++) {
    new_vals[i] = vals[i]; // setting second half to null
    new_dests[i] = dests[i]; // setting second half to null
  }
  edges.dests = new_dests;
  edges.vals = new_vals;
  vals = edges.vals;
  dests = edges.dests;
  */
  //uint32_t *new_dests = (uint32_t *)malloc(new_N * sizeof(*(edges.dests)));
  //uint32_t *new_vals = (uint32_t *)malloc(new_N * sizeof(*(edges.vals)));
  uint32_t *new_dests = (uint32_t *)aligned_alloc(32, new_N * sizeof(*(edges.dests)));
  uint32_t *new_vals = (uint32_t *)aligned_alloc(32, new_N * sizeof(*(edges.vals)));
	//memset(new_dests, 0, new_N * sizeof(*(edges.dests)));
	//memset(new_vals, 0, new_N * sizeof(*(edges.vals)));
  uint32_t * vals = (uint32_t *)edges.vals;
  uint32_t * dests = (uint32_t *)edges.dests;

  assert(edges.dests != NULL && edges.vals != NULL);
  
  
  //printf("List doubled: N - %u, logN = %u, H = %u\n", edges.N, edges.logN, edges.H);
  if (sub_counts.size() == 0) {
    uint_t start = 0;
    for (uint_t i = 0; i < edges.N; i++) {
      if (dests[i] != NULL_VAL) {
        new_vals[start] = vals[i];
        new_dests[start] = dests[i];
        start += 1;
      }
    }
		for (uint_t i = start; i < new_N; i++) {
			new_vals[i] = 0; new_dests[i] = NULL_VAL;
		}
    free(vals);
    free(dests);
    edges.vals = new_vals;
    edges.dests = new_dests;
    edges.N = new_N;
    redistribute(0, new_N);
  } else {
    parallel_for(size_t k = 0; k < sub_counts.size(); k+=1) {
      //printf("index = %u, len = %u, k = %u REDISTRIBUTE_PAR_SIZE = %u\n", index, len, k , REDISTRIBUTE_PAR_SIZE);
      uint_t start;
      if (k == 0) {
        start = 0;
      } else {
        start = sub_counts[k-1];
      }

      for (uint_t i = REDISTRIBUTE_PAR_SIZE*k; i < REDISTRIBUTE_PAR_SIZE*(k+1); i+=8) {
        //printf("start = %u, i = %u\n", start, i);
        for (uint_t j = i; j < i+8; j++) {
          if (dests[j] != NULL_VAL) {
            new_vals[start] = vals[j];
            new_dests[start] = dests[j];
            start += 1;
          }
        }
      }
    }
		for (uint_t i = sub_counts[sub_counts.size()-1]; i < new_N; i++) {
			new_vals[i] = 0; new_dests[i] = NULL_VAL;
		}
    free(vals);
    free(dests);
    edges.vals = new_vals;
    edges.dests = new_dests;
    edges.N = new_N;
    redistribute(0, new_N);
  }
  release_all_locks(task_id, true, GENERAL);
}


// index is the beginning of the sequence that you want to slide right.
// we wil always hold locks to the end of the leaf so we don't need to lock here
void PMA::slide_right(uint_t index, uint32_t *vals, uint32_t *dests) {
  // edges.list_lock.lock_shared();
  uint32_t el_val = vals[index];
  uint32_t el_dest = dests[index];
  dests[index] = NULL_VAL;
  vals[index] = 0;
  // uint32_t original_index = index;
  //printf("start of slide right, original_index: %d, worker number: %lu, \n", original_index, get_worker_num());
  // edges.lock_array[current_lock].print();
  index++;
  // uint32_t leaf = find_leaf(&edges, index);
	//uint32_t cnt{0};
  while (index < edges.N && (dests[index]!=NULL_VAL)) {
    // assert(find_leaf(&edges, index) == leaf);
    uint32_t temp_val = vals[index];
    uint32_t temp_dest = dests[index];
    vals[index] = el_val;
    dests[index] = el_dest;
    if (el_dest == SENT_VAL) {
      // fixing pointer of node that goes to this sentinel
      uint32_t node_index = el_val;
      fix_sentinel(node_index, index);
    }
    el_val = temp_val;
    el_dest = temp_dest;
    index++;
		//cnt++;
  }
	//std::cout << "right: " << cnt << '\n';
  if (el_dest == SENT_VAL) {
    // fixing pointer of node that goes to this sentinel
    uint32_t node_index = el_val;
    fix_sentinel(node_index, index);
  }
  //printf("middle of slide right, original_index: %d, worker number: %lu, current lock = %d\n", original_index, get_worker_num(), current_lock);

  // TODO There might be an issue with this going of the end sometimes
  assert(index != edges.N);

  vals[index] = el_val;
  dests[index] = el_dest;
  // assert(find_leaf(&edges, index) == leaf);
  // assert(check_no_full_leaves(&edges, original_index, edges.logN));
}


// index is the beginning of the sequence that you want to slide left.
// the element we start at will be deleted
// we wil always hold locks to the end of the leaf so we don't need to lock here
void PMA::slide_left(uint_t index, uint32_t *vals, uint32_t *dests) {
  // edges.list_lock.lock_shared();
  // uint32_t original_index = index;
  //printf("start of slide right, original_index: %d, worker number: %lu, \n", original_index, get_worker_num());
  // edges.lock_array[current_lock].print();
  // uint32_t leaf = find_leaf(&edges, index);
	//uint32_t cnt{0};
  while (index+1 < edges.N) {
		//cnt++;
    // assert(find_leaf(&edges, index) == leaf);
    uint32_t temp_val = vals[index+1];
    uint32_t temp_dest = dests[index+1];
    vals[index] = temp_val;
    dests[index] = temp_dest;
    if (temp_dest == SENT_VAL) {
      // fixing pointer of node that goes to this sentinel
      uint32_t node_index = temp_val;
      fix_sentinel(node_index, index);
    }
    if (dests[index] == NULL_VAL) {
      break;
    }
    index++;
  }
	// std::cout << "left: " << cnt << '\n';

  assert(index != edges.N);

  // assert(find_leaf(&edges, index) == leaf);
  // assert(check_no_full_leaves(&edges, original_index, edges.logN));
}

// important: make sure start, end don't include sentinels
// returns the index of the smallest element bigger than you in the range
// [start, end)
// if no such element is found, returns end (because insert shifts everything to
// the right)
// assumes we already hold the list_lock and the relevant lock_array locks
uint_t binary_search(const edge_list_t *list, uint32_t elem_dest, uint32_t elem_val, uint_t start,
                       uint_t end) {
  uint32_t *dests = (uint32_t *) list->dests;

  uint_t mid = (start + end) / 2;
  // print_array(list);
  while (start + 1 < end) {
    
    __builtin_prefetch ((void *)&dests[(mid+end)/2], 0, 3);
    __builtin_prefetch ((void *)&dests[(start + mid)/2], 0, 3);
    // printf("start = %d, end = %d, dest = %d, mid = %d, val =%u\n", start, end, elem_dest, mid, elem_val);
    /*
    if (mid % list->logN > list->logN / 2) {
      // if beginning of next leaf is before end of binary search
      uint32_t temp = next_leaf(mid, list->loglogN);
      if(temp < end) {
        mid = temp;
      }
    }
    */
    uint32_t item_dest = dests[mid];

    //if is_null
    if (item_dest==NULL_VAL) {
      // first check the next leaf
      uint_t check = next_leaf(mid, list->loglogN);
      //TODO deal with check is null
      if (check > end) {
        end = mid;
        mid = (start + end) / 2;
        __builtin_prefetch ((void *)&dests[mid], 0, 3);
        continue;
      }
      // if is_null
      if (dests[check]==NULL_VAL) {
        uint_t early_check = find_prev_valid(dests, mid);
        // if we found the sentinel, go right after it
        if (dests[early_check] == SENT_VAL) {
          return early_check + 1;
        }
        if (early_check < start) {
          start = mid;
          mid = (start + end) / 2;
          __builtin_prefetch ((void *)&dests[mid], 0, 3);
          continue;
        }  
        check = early_check;
      } 
      // printf("check = %d\n", check);
      uint32_t dest = dests[check];
      if (elem_dest == dest) {
        // cleanup before return
        return check;
      } else if (elem_dest < dest) {
        end = find_prev_valid(dests, mid) + 1;

      } else {
        if (check == start) {
          start = check + 1;
        } else {
          start = check;
        }
        // otherwise, searched for item is more than current and we set start
      }
      mid = (start + end) / 2;
      __builtin_prefetch ((void *)&dests[mid], 0, 3);
      continue;
    }

    if (elem_dest < item_dest) {
      end = mid; // if the searched for item is less than current item, set end
      mid = (start + end) / 2;
    } else if (elem_dest > item_dest) {
      start = mid;
      mid = (start + end) / 2;
      // otherwise, sesarched for item is more than current and we set start
    } else if (elem_dest == item_dest) {  // if we found it, return
      // cleanup before return
      return mid;
    }
  }
  if (end < start) {
    start = end;
  }
  assert(start >= 0);
  //tbassert(end < list->N, "end: %u, list->N: %u\n", end, list->N);

  //trying to encourage the packed left property so if they are both null go to the left
  if ((dests[start]==NULL_VAL) && (dests[end]==NULL_VAL)) {
    end = start;
  }

  // handling the case where there is one element left
  // if you are leq, return start (index where elt is)
  // otherwise, return end (no element greater than you in the range)
  // printf("start = %d, end = %d, n = %d\n", start,end, list->N);
  if (elem_dest <= dests[start] && (dests[start]!=NULL_VAL)) {
    end = start;
  }
  // cleanup before return

  return end;
}

uint32_t PMA::find_value(uint32_t src, uint32_t dest) {

#if ENABLE_PMA_LOCK == 1
  uint64_t task_id = __sync_fetch_and_add(&next_task_id, 2);
  node_lock.lock_shared(task_id);
#endif
  uint32_t e_value = 0;
  uint32_t e_dest = dest;
  //printf("src: %d, dest: %d\n", src, dest);
  //printf("beginning: %d, end:%d\n", nodes[src].beginning+1, nodes[src].end);

#if ENABLE_PMA_LOCK == 1
  //the lock has been deleted
  // not sure why this has to be exclusive
  if (!nodes[src].lock.lock(task_id, GENERAL)) {
    node_lock.unlock_shared(task_id);
    return find_value(src,dest);
  }
#endif
  uint_t loc =
      binary_search(&edges, e_dest, e_value, nodes[src].beginning + 1, nodes[src].end);
  //printf("loc = %d, looking for %u, %u, found %u, %u\n",loc, src, dest, edges.dests[loc], edges.vals[loc]);
  e_dest = edges.dests[loc];
  e_value = edges.vals[loc];
#if ENABLE_PMA_LOCK == 1
  nodes[src].lock.unlock(task_id, GENERAL);
  // edges.list_lock.unlock_shared();
  node_lock.unlock_shared(task_id);
#endif
  // print_array();
  // printf("loc: %d\n", loc);
  //TODO probably don't need the first check since we will never look for null
  if ((e_dest != NULL_VAL) && e_dest == dest) {
    return e_value;
  } else {
    return 0;
  }
}

//assumes node_lock is held
uint32_t PMA::find_contaning_node(uint_t index) {
  uint32_t start = 0; 
  uint32_t end = nodes.size()-1;
  while (end - start > 1) {
    uint32_t middle = (end + start) / 2;
    uint_t node_start = nodes[middle].beginning;
    uint_t node_end = nodes[middle].end;
    if ( index >= node_start && index < node_end){
      return middle;
    } else if (index < node_start) {
      end = middle;
    } else if (index >= node_end) {
      start = middle;
    } else {
      printf("should not happen\n");
      assert(false);
    }
  }
  if ( index >= nodes[start].beginning && index < nodes[start].end){
    return start;
  } else if ( index >= nodes[end].beginning && index < nodes[end].end) {
    return end;
  } else if (index >= nodes[nodes.size() - 1].end) {
      return nodes.size() - 1;
  } else {
    //printf("no containing node trying again\n");
    return find_contaning_node(index);
  }
}

pair_int PMA::which_locks_in_range(uint_t index, uint_t len, uint32_t guess) {
  uint32_t start;
which_locks_in_range_start:

  if (index >= nodes[guess].beginning && index < nodes[guess].end) {
      start = guess;
  } else {
    bool found = false;
    uint64_t range = 50;
    uint64_t st = (guess > range) ? guess - range : 0;
    uint64_t end = (guess + range < nodes.size()) ? guess + range : nodes.size() - 1;
    if (index >= nodes[st].beginning && index < nodes[end].end) {
      for (uint64_t est = st; est < end; est++) {
        if (index >= nodes[est].beginning && index < nodes[est].end) {
          start = est;
          found = true;
          break;
        }
      }
    }
    if (!found) {
      start = find_contaning_node(index);
    }
  }
  if (index == nodes[start].beginning) {
    if (start > 0) {
      start--;
    }
  }

  uint64_t end_index = next_leaf(index + len, edges.loglogN) - 1;
  uint32_t end;
  if (end_index < nodes[start].end) {
    end = start;
  } else {
    bool found = false;
    if (len <= 512) {
      uint64_t range = 200;
      uint64_t end_range = (start + range < nodes.size()) ? start + range : nodes.size() - 1;
      for (uint64_t est = start; est < end_range; est++) {
        if (end_index >= nodes[est].beginning && end_index < nodes[est].end) {
          end = est;
          found = true;
          break;
        }
      }
    }
    if (!found) {
      end = find_contaning_node(end_index);
    }
  }
  //printf("grabbing locks %d through %d, by %lu with reason %d\n", start, end, get_worker_num(), reason)
  if (end < start) {
    goto which_locks_in_range_start;
  }
  return {start, end};
}

pair_int PMA::grab_locks_in_range_with_resets(uint64_t task_id, uint_t index, uint_t len, REASONS reason, uint32_t guess) {
  uint_t orig_n = edges.N;
#if ENABLE_PMA_LOCK == 1
  start_grab_locks_in_range_with_resets:
#endif
  if (orig_n != edges.N) {
    return {0xFFFFFFFF,0};
  }
  pair_int ends = which_locks_in_range(index, len, guess);
#if ENABLE_PMA_LOCK == 1
  //printf("grabbing locks %d through %d, by %lu with reason %d\n", ends.x, ends.y, get_worker_num(), reason);
  for (uint32_t i = ends.x; i <= ends.y; i++) {
    if (!nodes[i].lock.try_lock(task_id, reason)) {
      if ( i > 0) {
        release_locks_in_range(task_id, {ends.x, i - 1}, SAME);
      }
      //usleep(1);
      //return grab_locks_in_range_with_resets(task_id, index, len, reason);
      goto start_grab_locks_in_range_with_resets;
    }
  }
  //printf("grabed locks %d through %d, by %lu with reason %d\n", ends.x, ends.y, get_worker_num(), reason);
  // make sure nothing changed before we managed to grab the locks
  if (!((nodes[ends.x].beginning <= index) && ((index+len <= nodes[ends.y].end) || (ends.y == nodes.size() -1)))) {
    release_locks_in_range(task_id, ends, SAME);
    goto start_grab_locks_in_range_with_resets;
  }
#endif
  return ends;
}

void PMA::release_locks_in_range(uint64_t task_id, pair_int locks, REASONS reason) {
#if ENABLE_PMA_LOCK == 1
uint32_t start = locks.x;
  uint32_t end = locks.y;
  //printf("releasing locks %d through %d, by worker %lu with reason %d\n", start, end, get_worker_num(), reason);
  for (uint32_t i = start; i <= end; i++) {
    nodes[i].lock.unlock(task_id, reason);
  }
#endif
}

pair_int PMA::which_locks_for_leaf(uint32_t src) {
  uint_t start_index = find_leaf(&edges, nodes[src].beginning);
  uint_t end_index = next_leaf(nodes[src].end, edges.loglogN);
  uint32_t first_node = src;
  while (nodes[first_node].beginning > start_index) {
    first_node--;
  }
  uint32_t last_node = src;
  while (nodes[last_node].end < end_index && last_node < nodes.size() -1) {
    last_node++;
  }
  return {first_node, last_node};
}

pair_int PMA::grab_locks_for_leaf_with_resets(uint64_t task_id, uint32_t src, REASONS reason) {
#if ENABLE_PMA_LOCK == 1
start_grab_locks_for_leaf_with_resets:
#endif
  pair_int ends = which_locks_for_leaf(src);
#if ENABLE_PMA_LOCK == 1
  //printf("grabbing locks %d through %d, by %lu with reason %d\n", ends.x, ends.y, get_worker_num(), reason);
  for (uint32_t i = ends.x; i <= ends.y; i++) {
    if (!nodes[i].lock.try_lock(task_id, reason)) {
      if ( i > 0) {
        release_locks_in_range(task_id, {ends.x, i - 1}, SAME);
      }
      usleep(1);
     //sched_yield();
      goto start_grab_locks_for_leaf_with_resets;
    }
  }
  //printf("grabed locks %d through %d, by %lu with reason %d\n", ends.x, ends.y, get_worker_num(), reason);
  pair_int ends_check = which_locks_for_leaf(src);
  if (ends.x != ends_check.x || ends.y != ends_check.y) {
    release_locks_in_range(task_id, ends);
    goto start_grab_locks_for_leaf_with_resets;
  }
#endif
  return ends;
}

bool PMA::check_every_lock_in_leaf(uint64_t task_id, uint_t index) {
  uint_t start = find_leaf(&edges, index);
  uint_t end = next_leaf(index, edges.loglogN);
  for (uint_t i = start; i < end; i++) {
#if ENABLE_PMA_LOCK == 1
    uint32_t node = find_contaning_node(i);
    assert(nodes[node].lock.i_own_lock(task_id));
    if (!nodes[node].lock.i_own_lock(task_id)) {
      return false;
    }
#endif
  }
  return true;
}

bool PMA::check_every_lock_in_node(uint64_t task_id, uint_t index, uint_t len) {
#if ENABLE_PMA_LOCK == 1
  uint_t start = find_leaf(&edges, index);
  for (uint_t i = start; i < start + len; i++) {
    uint32_t node = find_contaning_node(i);
    if (!nodes[node].lock.i_own_lock(task_id)) {
      return false;
    }
    assert(nodes[node].lock.i_own_lock(task_id));
  }
#endif
  return true;
}


// insert elem at index
// assumes the lock on index is held by the parent
// and releases it when it is done with it
void PMA::insert(uint64_t task_id, uint_t index, uint32_t elem_dest, uint32_t elem_val, uint32_t src
  #if ENABLE_PMA_LOCK == 1
  , pair_int held_locks
  #endif
  ) {
  #if ENABLE_PMA_LOCK == 1
  assert(check_every_lock_in_leaf(task_id, index));
  #endif
  assert(check_no_full_leaves(&edges, find_leaf(&edges, index), edges.logN));
  
  /*
  if (index > 0 && edges.dests[index - 1] == NULL_VAL && find_leaf(&edges, index) != index) {
    printf("insert conditional || index = %d src = %u elem dest = %u, elem value = %u, next leaf: %d\n", index, src, elem_dest, elem_val, next_leaf(index, edges.loglogN));
    print_array(get_worker_num());
  }*/


  //printf("index = %d src = %u elem dest = %u, elem value = %u, worker %lu\n", index, src, elem.dest, elem.value, get_worker_num());
  //print_array(get_worker_num());
  uint_t orig_n = edges.N;
  
  if ((elem_dest != SENT_VAL) && (elem_val != 0)) {
    // -1 for the case that the location is the spot of the next sentinal
    //assert(src == find_contaning_node(index-1));
  }
  // edges.list_lock.lock_shared();
  #if ENABLE_PMA_LOCK == 1
  assert(nodes[src].lock.i_own_lock(task_id));
  #endif
  int level = edges.H;
  uint_t len = edges.logN;

  uint32_t * vals = (uint32_t *) edges.vals;
  uint32_t * dests = (uint32_t *) edges.dests;
  // always deposit on the left
  if (dests[index]==NULL_VAL) {
    // printf("added to empty\n");
    vals[index] = elem_val;
    dests[index] = elem_dest;
    #if ENABLE_PMA_LOCK == 1
    assert(check_every_lock_in_leaf(task_id, index));
    #endif

  } else {
    assert(index < edges.N - 1);
    // if the edge already exists in the graph, update its value
    // do not make another edge
    if ((elem_dest != SENT_VAL) && dests[index] == elem_dest) {
      vals[index] = elem_val;
      #if ENABLE_PMA_LOCK == 1
      assert(check_every_lock_in_leaf(task_id, index));
      #endif
      // edges.list_lock.unlock_shared();
      assert(check_no_full_leaves(&edges, find_leaf(&edges, index), edges.logN));
      #if ENABLE_PMA_LOCK == 1
      release_locks_in_range(task_id, held_locks);
      assert(check_no_node_locks_for_me(task_id));
      #endif
      return;
    } else {
      // slide right assumes we hold the lock till the end of the leaf
      #if ENABLE_PMA_LOCK == 1
      assert(check_every_lock_in_leaf(task_id, index));
      #endif
      slide_right(index, vals, dests);
      // printf("after sliding, index = %d\n", index);
      // print_array();
      vals[index] = elem_val;
      dests[index] = elem_dest;
      #if ENABLE_PMA_LOCK == 1
      assert(check_every_lock_in_leaf(task_id, index));
      #endif
      // print_array();
    }
  }
  
  //assert(vals[index]!=0);
  //print_array();
  uint_t node_index = find_leaf(&edges, index);
  double density = get_density(&edges, node_index, len);
  //printf("density = %f, %d\n", density, density == 1);

  // spill over into next level up, node is completely full.
  if (density == 1) {
    //printf("first rebalence\n");
    //TODO continue work here
    len*=2;
    level--;
    node_index = find_node(node_index, len);
  } else {
    redistribute(node_index, len);
  }
  assert(edges.N == orig_n);
  #if ENABLE_PMA_LOCK == 1
  assert(check_every_lock_in_leaf(task_id, index));
  // being fancy here since some of the locks are done with and others might be used for the rebalence
  pair_int locks_for_rebalance = which_locks_in_range(node_index, len, src);
  //printf("relasing locks %d through %d with reason %d\n", locks_for_leaf.x, locks_for_rebalance.x -1, GENERAL);
  //printf("before with worker %lu\n", get_worker_num());
  //printf("worker %lu unlocking %d through %d with GENERAL \n", get_worker_num(), locks_for_leaf.x, locks_for_rebalance.x);
  parallel_for (uint32_t i = held_locks.x; i < locks_for_rebalance.x; i++) {
    nodes[i].lock.unlock(task_id);
  }
  //printf("worker %lu relasing locks %d through %d with reason REBALANCE\n",get_worker_num(), max(locks_for_rebalance.x, locks_for_leaf.x), locks_for_leaf.y);
  for (uint32_t i = max(locks_for_rebalance.x, held_locks.x) ; i <= min(held_locks.y, locks_for_rebalance.y); i++) {
    nodes[i].lock.unlock(task_id, REBALANCE);
  }
  //printf("worker %lu unlocking %d through %d with GENERAL \n", get_worker_num(), locks_for_rebalance.y +1,locks_for_leaf.y);

  for (uint32_t i = locks_for_rebalance.y +1; i <= held_locks.y; i++) {
    nodes[i].lock.unlock(task_id);
  }

  //printf("unlocking %d by worker %lu\n",src, get_worker_num());
  //printf("grabbing %d, %d\n", node_index, len);
  //printf("about to grab\n");
  pair_int lock_span = {max(locks_for_rebalance.x, held_locks.x),  max(locks_for_rebalance.y, held_locks.y)};
  assert(check_no_node_locks_held_by_me(task_id));


  held_locks = grab_locks_in_range_with_resets(task_id, node_index, len, REBALANCE, src);
  lock_span.x = min(held_locks.x, lock_span.x);
  lock_span.y = max(held_locks.y, lock_span.y);
  #endif
  // if somebody doubled and we let them with a reset
  if (edges.N != orig_n) {
    #if ENABLE_PMA_LOCK == 1
    //printf("worker %lu: somebody else doubled on us so we assume it has been distributed properly 1\n", get_worker_num());
    //printf("worker %lu: edges.N = %u, orig_n = %u\n",get_worker_num(), edges.N, orig_n);

    release_locks_in_range(task_id, held_locks, GENERAL);
    for (uint32_t i = lock_span.x; i <= lock_span.y; i++) {
      if (i < held_locks.x || i > held_locks.y) {
        nodes[i].lock.lock(task_id, REBALANCE);
        nodes[i].lock.unlock(task_id);
      }
    }
    assert(check_no_node_locks_for_me(task_id));
    #endif
    return;
  }
  // printf("node_index3 = %d\n", node_index);
  // print_array();

  // get density of the leaf you are in
  double density_b = upper_density_bound[level];
  uint_t density_count = get_density_count(&edges, node_index, len);
  density = ((double)density_count)/len;
  // printf("density %f, upperbound %f, len = %d, N = %d, logN = %d\n", density,
  // density_b.y, len, list->N, list->logN);

  // while density too high, go up the implicit tree
  // go up to the biggest node above the density bound
  //printf("node_index = %d, desnsity = %f, density bound = %f\n", node_index, density, density_b.y);
  std::vector<uint_t> sub_counts(0);
  // while (density >= density_b) {
  while(density >= density_b) {
  // density >= (((double) edges.logN - 1) / edges.logN)) {

    //printf("node_index = %d, desnsity = %f, density bound = %f, len = %d, worker = %lu\n", node_index, density, density_b.y, len, get_worker_num());
    len *= 2;
    if (len <= edges.N) {
      level--;
      uint_t new_node_index = find_node(node_index, len);
      #if ENABLE_PMA_LOCK == 1
      release_locks_in_range(task_id, held_locks, REBALANCE);

      held_locks = grab_locks_in_range_with_resets(task_id, new_node_index, len, REBALANCE, src);
      if (held_locks.x < lock_span.x) {
        lock_span.x = held_locks.x;
      }
      if (held_locks.y > lock_span.y) {
        lock_span.y = held_locks.y;
      }
      #endif
      if (edges.N != orig_n) {
        //printf("worker %lu: somebody else doubled on us so we assume it has been distributed properly 2\n", get_worker_num());
        //printf("worker %lu: edges.N = %u, orig_n = %u\n",get_worker_num(), edges.N, orig_n);
        #if ENABLE_PMA_LOCK == 1
        release_locks_in_range(task_id, held_locks, GENERAL);
        for (uint32_t i = lock_span.x; i <= lock_span.y; i++) {
          if (i < held_locks.x || i > held_locks.y) {
            nodes[i].lock.lock(task_id, REBALANCE);
            nodes[i].lock.unlock(task_id);
          }
        }
        assert(check_no_node_locks_for_me(task_id));
        #endif
        return;
      }
      if (len <= REDISTRIBUTE_PAR_SIZE) {
        density_count = get_density_count(&edges, new_node_index, len);
      } else {
        sub_counts.resize(len/REDISTRIBUTE_PAR_SIZE);
        density_count = get_density_count_par(&edges, new_node_index, len, sub_counts);
      }
      // to help prevent double doubling by knowing how big it was on the last count
      orig_n = edges.N;
      node_index = new_node_index;
      density_b = upper_density_bound[level];
      density = ((double) density_count )/len;
    } else {
      // if you reach the root, double the list
      #if ENABLE_PMA_LOCK == 1
      release_locks_in_range(task_id, held_locks, DOUBLE);
      //if (edges.N == orig_n) {
        // for (int i = 0; i < lock_count; i++) {
        //  edges.list_lock.unlock_shared();
        //}
        //TODO don't double double
        //printf("second double due to worker %lu\n", get_worker_num());
        //print_array(get_worker_num());
        assert(check_no_node_locks_held_by_me(task_id));
        #endif
        double_list(task_id, sub_counts, density_count);
        // -1 for the lock at the begining of the function
        // we want to leave with the same number as we entered with
        //for (int i = 0; i < lock_count-1; i++) {
        //  edges.list_lock.lock_shared();
        //}
     // }
        #if ENABLE_PMA_LOCK == 1
        assert(check_no_node_locks_for_me(task_id));
        #endif
      return;
    }
    // printf("density %f, upperbound %f, len = %d, N = %d, logN = %d\n",
    // density, density_b.y, len, list->N, list->logN);
  }
  assert(((double)get_density_count(&edges, node_index, len))/ len == density);
  assert(density < density_b);
  assert(density <= (((double) edges.logN - 1) / edges.logN));
  //print_array(get_worker_num());
  if(len > edges.logN) {
     if (len <= REDISTRIBUTE_PAR_SIZE) {
      redistribute(node_index, len);
     } else {
      redistribute_par(node_index, len, sub_counts, density_count);
     }
    
  }
  //printf("this relase? %d, %d\n", node_index, len);
  assert(check_no_full_leaves(&edges, node_index, len));
  #if ENABLE_PMA_LOCK == 1
  release_locks_in_range(task_id, held_locks, GENERAL);

  for (uint32_t i = lock_span.x; i < held_locks.x; i++) {
    nodes[i].lock.lock(task_id, REBALANCE);
    nodes[i].lock.unlock(task_id);
  }
  for (uint32_t i = held_locks.y+1; i <= lock_span.y; i++) {
    nodes[i].lock.lock(task_id, REBALANCE);
    nodes[i].lock.unlock(task_id);
  }
  // printf("node_index5 = %d\n", node_index);
  // print_array();
  // edges.list_lock.unlock_shared();
    assert(check_no_node_locks_for_me(task_id));
  #endif
  return;
}



// remove elem at index
// assumes the lock on index is held by the parent
// and releases it when it is done with it
void PMA::remove(uint64_t task_id, uint_t index, uint32_t elem_dest, uint32_t src
  #if ENABLE_PMA_LOCK == 1
  , pair_int held_locks
  #endif
  ) {
  #if ENABLE_PMA_LOCK == 1
  assert(check_every_lock_in_leaf(task_id, index));
  assert(check_no_full_leaves(&edges, find_leaf(&edges, index), edges.logN));
  #endif
  
  /*
  if (index > 0 && edges.dests[index - 1] == NULL_VAL && find_leaf(&edges, index) != index) {
    printf("insert conditional || index = %d src = %u elem dest = %u, elem value = %u, next leaf: %d\n", index, src, elem_dest, elem_val, next_leaf(index, edges.loglogN));
    print_array(get_worker_num());
  }*/


  //printf("reoving: index = %d src = %u elem dest = %u\n", index, src, elem_dest);
  //print_array(get_worker_num());
  uint_t orig_n = edges.N;
  
  if ((elem_dest != SENT_VAL) ) {
    // -1 for the case that the location is the spot of the next sentinal
    assert(src == find_contaning_node(index-1));
  }
  // edges.list_lock.lock_shared();
  #if ENABLE_PMA_LOCK == 1
  assert(nodes[src].lock.i_own_lock(task_id));
  #endif
  int level = edges.H;
  uint_t len = edges.logN;

  uint32_t * vals = (uint32_t *) edges.vals;
  uint32_t * dests = (uint32_t *) edges.dests;
  // always deposit on the left
  assert(index < edges.N - 1);
  // if the edge already exists in the graph, update its value
  // do not make another edge
  // slide right assumes we hold the lock till the end of the leaf
  #if ENABLE_PMA_LOCK == 1
  assert(check_every_lock_in_leaf(task_id, index));
  #endif
  //printf("before sliding, index = %d\n", index);
  //print_array();
  slide_left(index, vals, dests);
  //printf("after sliding, index = %d\n", index);
  //print_array();

  
  //assert(vals[index]!=0);
  //print_array();
  uint_t node_index = find_leaf(&edges, index);
  double density = get_density(&edges, node_index, len);
  //printf("density = %f, %d\n", density, density == 1);

  // spill over into next level up, node is completely full.
  if (density == 0) {
    //printf("first rebalence\n");
    len*=2;
    level--;
    node_index = find_node(node_index, len);
  } else {
    redistribute(node_index, len);
  }
  assert(edges.N == orig_n);
  #if ENABLE_PMA_LOCK == 1
  assert(check_every_lock_in_leaf(task_id, index));
  // being fancy here since some of the locks are done with and others might be used for the rebalence
  pair_int locks_for_rebalance = which_locks_in_range(node_index, len, src);
  //printf("relasing locks %d through %d with reason %d\n", locks_for_leaf.x, locks_for_rebalance.x -1, GENERAL);
  //printf("before with worker %lu\n", get_worker_num());
  //printf("worker %lu unlocking %d through %d with GENERAL \n", get_worker_num(), locks_for_leaf.x, locks_for_rebalance.x);
  parallel_for (int i = held_locks.x; i < locks_for_rebalance.x; i++) {
    nodes[i].lock.unlock(task_id);
  }
  //printf("worker %lu relasing locks %d through %d with reason REBALANCE\n",get_worker_num(), max(locks_for_rebalance.x, locks_for_leaf.x), locks_for_leaf.y);
  for (int i = max(locks_for_rebalance.x, held_locks.x) ; i <= min(held_locks.y, locks_for_rebalance.y); i++) {
    nodes[i].lock.unlock(task_id, REBALANCE);
  }
  //printf("worker %lu unlocking %d through %d with GENERAL \n", get_worker_num(), locks_for_rebalance.y +1,locks_for_leaf.y);

  for (int i = locks_for_rebalance.y +1; i <= held_locks.y; i++) {
    nodes[i].lock.unlock(task_id);
  }

  //printf("unlocking %d by worker %lu\n",src, get_worker_num());
  //printf("grabbing %d, %d\n", node_index, len);
  //printf("about to grab\n");
  pair_int lock_span = {max(locks_for_rebalance.x, held_locks.x),  max(locks_for_rebalance.y, held_locks.y)};
  assert(check_no_node_locks_held_by_me(task_id));

  held_locks = grab_locks_in_range_with_resets(task_id, node_index, len, REBALANCE, src);
  lock_span.x = min(held_locks.x, lock_span.x);
  lock_span.y = max(held_locks.y, lock_span.y);
  #endif

  // if somebody doubled and we let them with a reset
  if (edges.N != orig_n) {
    //printf("worker %lu: somebody else doubled on us so we assume it has been distributed properly 1\n", get_worker_num());
    //printf("worker %lu: edges.N = %u, orig_n = %u\n",get_worker_num(), edges.N, orig_n);
    #if ENABLE_PMA_LOCK == 1
    release_locks_in_range(task_id, held_locks, GENERAL);
    for (int i = lock_span.x; i <= lock_span.y; i++) {
      if (i < held_locks.x || i > held_locks.y) {
        nodes[i].lock.lock(task_id, REBALANCE);
        nodes[i].lock.unlock(task_id);
      }
    }
    assert(check_no_node_locks_for_me(task_id));
    #endif
    return;
  }
  // printf("node_index3 = %d\n", node_index);
  // print_array();

  // get density of the leaf you are in
  double density_b = lower_density_bound[level];
  uint_t density_count = get_density_count(&edges, node_index, len);
  density = ((double)density_count)/len;
  // printf("density %f, upperbound %f, len = %d, N = %d, logN = %d\n", density,
  // density_b.y, len, list->N, list->logN);

  // while density too high, go up the implicit tree
  // go up to the biggest node above the density bound
  //printf("node_index = %d, desnsity = %f, density bound = %f\n", node_index, density, density_b.y);
  std::vector<uint_t> sub_counts(0);
  while (density <= density_b) {
    //printf("node_index = %d, desnsity = %f, density bound = %f, len = %d, worker = %lu\n", node_index, density, density_b.y, len, get_worker_num());
    len *= 2;
    if (len <= edges.N) {
      level--;
      uint_t new_node_index = find_node(node_index, len);
      
      
      #if ENABLE_PMA_LOCK == 1
      release_locks_in_range(task_id, held_locks, REBALANCE);
      held_locks = grab_locks_in_range_with_resets(task_id, new_node_index, len, REBALANCE, src);
      if (held_locks.x < lock_span.x) {
        lock_span.x = held_locks.x;
      }
      if (held_locks.y > lock_span.y) {
        lock_span.y = held_locks.y;
      }
      #endif
      if (edges.N != orig_n) {
        #if ENABLE_PMA_LOCK == 1
        //printf("worker %lu: somebody else doubled on us so we assume it has been distributed properly 2\n", get_worker_num());
        //printf("worker %lu: edges.N = %u, orig_n = %u\n",get_worker_num(), edges.N, orig_n);
        release_locks_in_range(task_id, held_locks, GENERAL);
        for (int i = lock_span.x; i <= lock_span.y; i++) {
          if (i < held_locks.x || i > held_locks.y) {
            nodes[i].lock.lock(task_id, REBALANCE);
            nodes[i].lock.unlock(task_id);
          }
        }
        assert(check_no_node_locks_for_me(task_id));
        #endif
        return;
      }
      if (len <= REDISTRIBUTE_PAR_SIZE) {
        density_count = get_density_count(&edges, new_node_index, len);
      } else {
        sub_counts.resize(len/REDISTRIBUTE_PAR_SIZE);
        density_count = get_density_count_par(&edges, new_node_index, len, sub_counts);
      }
      // to help prevent double doubling by knowing how big it was on the last count
      orig_n = edges.N;
      node_index = new_node_index;
      density_b = lower_density_bound[level];
      density = ((double) density_count )/len;
    } else {
      // if you reach the root, double the list
      #if ENABLE_PMA_LOCK == 1
      release_locks_in_range(task_id, held_locks, DOUBLE);
      #endif
      //if (edges.N == orig_n) {
        // for (int i = 0; i < lock_count; i++) {
        //  edges.list_lock.unlock_shared();
        //}
        //TODO don't double double
        //printf("second double due to worker %lu\n", get_worker_num());
        //print_array(get_worker_num());
        assert(check_no_node_locks_held_by_me(task_id));
        half_list(task_id, sub_counts, density_count);
        // -1 for the lock at the begining of the function
        // we want to leave with the same number as we entered with
        //for (int i = 0; i < lock_count-1; i++) {
        //  edges.list_lock.lock_shared();
        //}
     // }
        assert(check_no_node_locks_for_me(task_id));
      return;
    }
    // printf("density %f, upperbound %f, len = %d, N = %d, logN = %d\n",
    // density, density_b.y, len, list->N, list->logN);
  }
  assert(((double)get_density_count(&edges, node_index, len))/ len == density);
  assert(density > density_b);
  //print_array(get_worker_num());
  if(len > edges.logN) {
     if (len <= REDISTRIBUTE_PAR_SIZE) {
      redistribute(node_index, len);
     } else {
      redistribute_par(node_index, len, sub_counts, density_count);
     }
    
  }
  //printf("this relase? %d, %d\n", node_index, len);
  assert(check_no_full_leaves(&edges, node_index, len));
  #if ENABLE_PMA_LOCK == 1
  release_locks_in_range(task_id, held_locks, GENERAL);

  for (uint32_t i = lock_span.x; i < held_locks.x; i++) {
    nodes[i].lock.lock(task_id, REBALANCE);
    nodes[i].lock.unlock(task_id);
  }
  for (uint32_t i = held_locks.y+1; i <= lock_span.y; i++) {
    nodes[i].lock.lock(task_id, REBALANCE);
    nodes[i].lock.unlock(task_id);
  }
  // printf("node_index5 = %d\n", node_index);
  // print_array();
  // edges.list_lock.unlock_shared();
    assert(check_no_node_locks_for_me(task_id));
    #endif
  return;
}


#include <immintrin.h>
#include <iostream>
#include <iomanip>    

template<class T> inline void Log(const __m256i & value) {
    const size_t n = sizeof(__m256i) / sizeof(T);
    T buffer[n];
    _mm256_storeu_si256((__m256i*)buffer, value);
    for (int i = 0; i < n; i++)
        std::cout << buffer[i] << " ";
    printf("\n");
}


/*
//I don't think this is quite correct, but it did work at some point
std::vector<uint32_t>
PMA::sparse_matrix_vector_multiplication(std::vector<uint32_t> const &v) {
  std::vector<uint32_t> result(nodes.size(), 0);

  uint32_t num_vertices = nodes.size();
  
  bool vector = true;

  parallel_for (uint32_t i = 0; i < num_vertices; i++) {
    nodes[i].lock.lock_shared();
    // printf("i: %d\n", i);
    // +1 to avoid sentinel
    uint32_t j = nodes[i].beginning + 1;
    __m256i temp_sum = _mm256_setzero_si256();
    while (j < nodes[i].end) {
       if (!vector || (nodes[i].end - j < 4 || j >= next_leaf(j, edges.logN)-4)) { 
        if (!is_null(edges.items[j].e)) {
          result[i] += edges.items[j].e.value * v[edges.items[j].e.dest];
          j++;
        } else {
          j = next_leaf(j, edges.logN);
        }
      
      } else {
        __m256i edgegroup = _mm256_loadu_si256((__m256i*) &edges.items[j]);
        __m256i values_group = _mm256_srli_epi64(edgegroup, 32);
        __m256i mask = _mm256_setr_epi64x(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF);
        __m256i dest_group_targets = _mm256_and_si256(edgegroup, mask);
        __m256i dest_group = _mm256_i32gather_epi32((int*)v.data(), dest_group_targets, sizeof(v[0]));
        __m256i partial_sums = _mm256_mul_epi32(dest_group, values_group);
        temp_sum = _mm256_add_epi64(temp_sum, partial_sums);
        //temp_sums = _mm256_and_si256(temp_sums, mask);    


        __m256i mask2 = _mm256_setr_epi64x(0xFFFFFFFFFFFFFFFFUL, 0xFFFFFFFFFFFFFFFFUL, 0xFFFFFFFFFFFFFFFFUL, 0xFFFFFFFF00000000UL);
        if (_mm256_testc_si256(mask2, edgegroup)) {
          j = next_leaf(j, edges.logN);
        } else {
          j+=4;
        }
      }
    }
    if (vector) {
      result[i] += _mm256_extract_epi32(temp_sum, 0);
      result[i] += _mm256_extract_epi32(temp_sum, 2);
      result[i] += _mm256_extract_epi32(temp_sum, 4);
      result[i] += _mm256_extract_epi32(temp_sum, 6);
    } 
  } 
  return result;
}
*/

void PMA::print_graph() {
  uint32_t num_vertices = nodes.size();
  for (uint32_t i = 0; i < num_vertices; i++) {
    // printf("i: %d\n", i);
    // +1 to avoid sentinel
    uint32_t matrix_index = 0;

    for (uint_t j = nodes[i].beginning + 1; j < nodes[i].end; j++) {
      if (edges.dests[j]!=NULL_VAL) {
        while (matrix_index < edges.dests[j]) {
          printf("000 ");
          matrix_index++;
        }
        printf("%03d ", edges.vals[j]);
        matrix_index++;
      }
    }
    for (uint32_t j = matrix_index; j < num_vertices; j++) {
      printf("000 ");
    }
    printf("\n");
  }
}

// add a node to the graph
void PMA::add_node() {
  uint64_t task_id = __sync_fetch_and_add(&next_task_id, 2);
#if ENABLE_PMA_LOCK == 1
  node_lock.lock(task_id);
#endif
  node_t node;
  uint32_t len = nodes.size();
  edge_t sentinel;
  sentinel.dest = SENT_VAL; // placeholder
  sentinel.value = len;       // back pointer

  if (len > 0) {
    node.beginning = nodes[len - 1].end;
    if(node.beginning == edges.N - 1) {
      uint32_t volatile  * volatile dests = edges.dests;
      uint_t leaf = find_leaf(&edges, node.beginning);
      // moving the beginning of the node you are inserting left
      //TODO jump back to the last leaf and look forward from there
      if (nodes[len-1].num_neighbors == 0) {
        node.beginning = nodes[len - 1].beginning + 1;
      } else {
        while(dests[node.beginning - 1] == NULL_VAL && node.beginning != leaf /* && next_leaf(node.beginning, edges.loglogN) != node.beginning*/) {
          node.beginning -= 1;
        }
      }
    }
    node.end = node.beginning + 1;
    // fix previous node to set its end to your beginning
    nodes[len - 1].end = node.beginning;
  } else {
    node.beginning = 0;
    node.end = 1;
    // be special to the first one since we know it never moves do it doesn't need to look like a sentinal since it doesn't have to be fixed ever
    sentinel.value = NULL_VAL;
    sentinel.dest = 0;
  }
  // printf("sentinel dest: %u, value: %u\n", sentinel.dest, sentinel.value);
  node.num_neighbors = 0;
#if ENABLE_PMA_LOCK == 1
  node.lock.name("nodeLock");
  node.lock.number(nodes.size());
#endif
  nodes.push_back(node);
  #if ENABLE_PMA_LOCK == 1
  pair_int held_locks = grab_locks_for_leaf_with_resets(task_id, nodes.size() - 1);
  // edges.list_lock.lock_shared();
  node_lock.make_shared(task_id);
  #endif
  uint_t loc = node.beginning;
  //tbassert(loc < edges.N, "loc: %d, edges.N: %d\n", loc, edges.N);
  insert(task_id, loc, sentinel.dest, sentinel.value, nodes.size() - 1
  #if ENABLE_PMA_LOCK == 1
    , held_locks
  #endif
  );
  // printf("end of insert reason: %d\n", edges.lock_array[find_lock_index(&edges, loc)].reason);
  // edges.list_lock.unlock_shared();
#if ENABLE_PMA_LOCK == 1
  node_lock.unlock_shared(task_id);
#endif
  // print_array(len);
}

//TODO deal with the case that multiple threads do the binary search and try and make the same region exclusive
void PMA::add_edge(uint32_t src, uint32_t dest, uint32_t value) {
  // printf("src = %d, dest = %d, val = %d\n", src,dest,value);
  if (value != 0) {
    uint64_t task_id = __sync_fetch_and_add(&next_task_id, 2);
    #if ENABLE_PMA_LOCK == 1
    assert(check_no_locks_for_me(task_id));
    node_lock.lock_shared(task_id);
    pair_int held_locks = grab_locks_for_leaf_with_resets(task_id, src);
    assert(nodes[src].lock.i_own_lock(task_id));
    #endif
    node_t node = nodes[src];


    edge_t e;
    e.dest = dest;
    e.value = value;
    
    // edges.list_lock.lock_shared();
    //uint32_t node_begin = node.beginning;
    //uint32_t node_end = node.end;
    __sync_fetch_and_add(&nodes[src].num_neighbors,1);
    // printf("looking in region %u, %u by worker %lu\n", node.beginning + 1, node.end, get_worker_num());
    uint_t loc_to_add =
        binary_search(&edges, e.dest, e.value, node.beginning + 1, node.end);
    __builtin_prefetch ((void *)&edges.vals[loc_to_add], 1, 3);
#if ENABLE_PMA_LOCK == 1
    pair_int needed_locks = which_locks_in_range(find_leaf(&edges, loc_to_add), edges.logN, src);
    for (uint32_t i = held_locks.x; i < needed_locks.x; i++) {
      nodes[i].lock.unlock(task_id, GENERAL);
      held_locks.x = i+1;
    }
    for (uint32_t i = needed_locks.y+1; i <= held_locks.y; i++) {
      nodes[i].lock.unlock(task_id, GENERAL);
      
    }
    if (needed_locks.y+1 <= held_locks.y) {
      held_locks.y = needed_locks.y;
    }
    // print_array();
    //printf("loc_to_add: %u by worker %lu\n", loc_to_add, get_worker_num());
    // printf("src: %d, dest: %d by worker %lu\n", src, dest, get_worker_num());
    // print_array();

    
    assert(nodes[src].lock.i_own_lock(task_id));
    assert(check_every_lock_in_leaf(task_id, loc_to_add));
    /*
    if (!check_no_full_leaves(&edges, find_leaf(&edges, loc_to_add), edges.logN)) {
      printf("worker %lu is trying to insert into a full leaf\n", get_worker_num());
      print_array(get_worker_num());
    }
    */
#endif
    assert(check_no_full_leaves(&edges, find_leaf(&edges, loc_to_add), edges.logN));
    insert(task_id, loc_to_add, e.dest, e.value, src
      #if ENABLE_PMA_LOCK == 1
      , held_locks
      #endif
     );
    //printf("worker %lu done\n", get_worker_num());
    // edges.list_lock.unlock_shared();
    //print_array();
    #if ENABLE_PMA_LOCK == 1
    node_lock.unlock_shared(task_id);
    assert(check_no_locks_for_me(task_id));
    #endif
  }
}

bool PMA::add_edge_update(uint32_t src, uint32_t dest, uint32_t value) {
  if (value != 0) {
#if ENABLE_PMA_LOCK == 1
    uint64_t task_id = __sync_fetch_and_add(&next_task_id, 2);
    node_lock.lock_shared(task_id);
    pair_int held_locks = grab_locks_for_leaf_with_resets(task_id, src);
    assert(nodes[src].lock.i_own_lock(task_id));
#else
    uint64_t task_id = next_task_id;
    next_task_id += 2;
#endif
    node_t node = nodes[src];


    edge_t e;
    e.dest = dest;
    e.value = value;
    
    // edges.list_lock.lock_shared();
    //uint32_t node_begin = node.beginning;
    //uint32_t node_end = node.end;
    // printf("looking in region %u, %u by worker %lu\n", node.beginning + 1, node.end, get_worker_num());
    uint_t loc_to_add =
        binary_search(&edges, e.dest, e.value, node.beginning + 1, node.end);
    //__builtin_prefetch ((void *)&edges.vals[loc_to_add], 1, 3);


    if (edges.dests[loc_to_add] == dest) {
      edges.vals[loc_to_add] = value;
      #if ENABLE_PMA_LOCK == 1
      release_locks_in_range(task_id, held_locks);
      node_lock.unlock_shared(task_id);
      #endif
      return false;
    }
    #if ENABLE_PMA_LOCK == 1
    pair_int needed_locks = which_locks_in_range(find_leaf(&edges, loc_to_add), edges.logN, src);
    for (uint32_t i = held_locks.x; i < needed_locks.x; i++) {
      nodes[i].lock.unlock(task_id, GENERAL);
      held_locks.x = i+1;
    }
    for (uint32_t i = needed_locks.y+1; i <= held_locks.y; i++) {
      nodes[i].lock.unlock(task_id, GENERAL);
      
    }
    if (needed_locks.y+1 <= held_locks.y) {
      held_locks.y = needed_locks.y;
    }
    #endif
    // printf("loc_to_add: %d\n", loc_to_add);
    // printf("src: %d, dest: %d\n", src, dest);
#if ENABLE_PMA_LOCK == 1
    __sync_fetch_and_add(&nodes[src].num_neighbors,1);
#else
    nodes[src].num_neighbors += 1;
#endif
    // print_array();
    insert(task_id, loc_to_add, e.dest, e.value, src
      #if ENABLE_PMA_LOCK == 1
      , held_locks
      #endif 
      );
    // edges.list_lock.unlock_shared();
    // print_array();
#if ENABLE_PMA_LOCK == 1
    node_lock.unlock_shared(task_id);
#endif
    //assert(check_no_locks_for_me(task_id));
    return true;
  }
  // if the value was 0
  return false;
}

bool __attribute__ ((noinline)) PMA::add_edge_update_fast(uint32_t src, uint32_t dest, uint32_t value, uint64_t task_id) {
    __builtin_prefetch(&nodes[src]);
    __builtin_prefetch(&edges);
    if (value != 0) {
      assert(check_no_locks_for_me(task_id));
      pair_int held_locks = grab_locks_for_leaf_with_resets(task_id, src);
      //assert(nodes[src].lock.i_own_lock(task_id));
      node_t node = nodes[src];


      edge_t e;
      e.dest = dest;
      e.value = value;
      
      // edges.list_lock.lock_shared();
      //uint32_t node_begin = node.beginning;
      //uint32_t node_end = node.end;
      // printf("looking in region %u, %u by worker %lu\n", node.beginning + 1, node.end, get_worker_num());
      uint_t loc_to_add =
          binary_search(&edges, e.dest, e.value, node.beginning + 1, node.end);
      __builtin_prefetch ((void *)&edges.vals[loc_to_add], 1, 3);


      if (edges.dests[loc_to_add] == dest) {
        edges.vals[loc_to_add] = value;
        release_locks_in_range(task_id, held_locks);
        return false;
      }
#if ENABLE_PMA_LOCK == 1
      pair_int needed_locks = which_locks_in_range(find_leaf(&edges, loc_to_add), edges.logN, src);
      for (uint32_t i = held_locks.x; i < needed_locks.x; i++) {
        nodes[i].lock.unlock(task_id, GENERAL);
        held_locks.x = i+1;
      }
      for (uint32_t i = needed_locks.y+1; i <= held_locks.y; i++) {
        nodes[i].lock.unlock(task_id, GENERAL);
        
      }
      if (needed_locks.y+1 <= held_locks.y) {
        held_locks.y = needed_locks.y;
      }
#endif
      // printf("loc_to_add: %d\n", loc_to_add);
      // printf("src: %d, dest: %d\n", src, dest);
      __sync_fetch_and_add(&nodes[src].num_neighbors,1);
      // print_array();
      insert(task_id, loc_to_add, e.dest, e.value, src
        #if ENABLE_PMA_LOCK == 1
        , held_locks
        #endif
        );
      // edges.list_lock.unlock_shared();
      // print_array();
      assert(check_no_locks_for_me(task_id));
      return true;
    }
  // if the value was 0
  return false;
}

void PMA::add_edge_batch_update_no_val(uint32_t *srcs, uint32_t *dests, uint_t edge_count) {
  uint64_t task_id = __sync_fetch_and_add(&next_task_id, 2*edge_count);
#if ENABLE_PMA_LOCK == 1
  node_lock.lock_shared(0);
  parallel_for (uint_t i = 0; i < edge_count; i++) {
#else
  for (uint_t i = 0; i < edge_count; i++) {
#endif
    uint32_t src = srcs[i];
    uint32_t dest = dests[i];
    add_edge_update_fast(src,dest,1,task_id+2*i);
  }
#if ENABLE_PMA_LOCK == 1
  node_lock.unlock_shared(0);
#endif
}
void PMA::add_edge_batch_update(uint32_t *srcs, uint32_t *dests, uint32_t *values, uint_t edge_count) {
  uint64_t task_id = __sync_fetch_and_add(&next_task_id, 2*edge_count);
#if ENABLE_PMA_LOCK == 1
  node_lock.lock_shared(0);
  parallel_for (uint_t i = 0; i < edge_count; i++) {
#else
  for (uint_t i = 0; i < edge_count; i++) {
#endif
    uint32_t src = srcs[i];
    uint32_t dest = dests[i];
    uint32_t value = values[i];
    add_edge_update_fast(src,dest,value,task_id+2*i);
  }
#if ENABLE_PMA_LOCK == 1
  node_lock.unlock_shared(0);
#endif
}

// merge a dense batch with PMA into a specified region of the giant output
// array
// return nnz in the range of the merge
uint64_t do_merge(uint32_t* new_dests, uint32_t* new_vals, uint64_t
output_start, pair_uint *es, uint64_t batch_start, uint64_t batch_end,
uint32_t* dests, uint32_t* vals, uint64_t pma_start, uint64_t pma_end, uint32_t loglogN)
{
  // dests = pma_dests, vals = pma_vals

  uint64_t es_idx = batch_start;
  uint64_t pma_idx = pma_start;
  uint64_t result_idx = output_start;
  // TODO: look around (could be a binary search if its too slow)
  uint32_t current_src_pma = 0;
  while(pma_idx > 0 && dests[pma_idx] != SENT_VAL) {
    pma_idx--;
  }
  if (pma_idx != 0) {
    current_src_pma = vals[pma_idx];
  }

  pma_idx = pma_start;
  if (pma_idx == 0) {
    new_dests[0] = dests[0];
    new_vals[0] = vals[0];
    result_idx++;
    pma_idx++;
  }

  // merge
  // printf("edges.N %lu, edge_count %u, max_elts %lu, new_size %lu, nnzs %lu\n", edges.N, edge_count, m
  while(es_idx < batch_end || pma_idx < pma_end) {
    // merge the rest of the PMA in
    if (es_idx == batch_end) {
      while(pma_idx < pma_end) {
        new_dests[result_idx] = dests[pma_idx];
        new_vals[result_idx] = vals[pma_idx];
        //printf("end, %u, %u\n", new_dests[result_idx], new_vals[result_idx]);
        result_idx+=(dests[pma_idx] != NULL_VAL);
        pma_idx++;
      }
      break;
    }

    uint32_t current_src_batch = es[es_idx].x;
    uint32_t current_dest_batch = es[es_idx].y;
    // skip over duplicates
    while(es_idx < batch_end && current_src_batch == es[es_idx+1].x && current_dest_batch == es[es_idx+1].y) {
      es_idx++;
    }
    assert(es_idx < batch_end);
    uint32_t current_val_batch = 1;

    assert(pma_idx < pma_end);

    // TODO: skip to next leaf if there is a performance issue
    while (pma_idx < pma_end && dests[pma_idx] == NULL_VAL) {
      pma_idx = next_leaf(pma_idx, loglogN);
      /*
      uint32_t logN = 1 << loglogN;
      if (current_src_pma < num_nodes) {
        while (dests[pma_idx] == NULL_VAL) {
          pma_idx += logN;
        }
      } else {
        while (pma_idx < pma_end && dests[pma_idx] == NULL_VAL) {
          pma_idx += logN;
        }
      }
      */
    }
    // at the last node
    if (pma_idx == pma_end) {
      // merge the rest of the batch in
      for(;es_idx < batch_end; es_idx++) {
        // TODO: remove later after correctness
        /*if (es[es_idx].x != current_src_pma) {
          printf("something is wrong\n");
        }*/
        // TODO: if you make edges SoA, you can do a memcpy here
        new_dests[result_idx] = es[es_idx].y;
        new_vals[result_idx] = current_val_batch;
        result_idx++;
      }
      break;
    }
    uint32_t current_dest_pma = dests[pma_idx];
    assert(current_dest_pma != NULL_VAL);
    // check if pma elt is a sentinel
    if (current_dest_pma == SENT_VAL) {
      // TODO: remove after correctness tests
      /*if (!(current_src_pma == 0 || current_src_pma == vals[pma_idx] - 1 || current_src_pma == vals[pm
        printf("something wrong with pma sentinals giving us a bad pma source, %u, %u\n", current_src_pm
      }*/
      current_src_pma = vals[pma_idx];
    }
    if (current_src_pma < current_src_batch) {
      //printf("1. %u, %u\n", current_src_pma, dests[pma_idx]);
      new_dests[result_idx] = dests[pma_idx];
      new_vals[result_idx] = vals[pma_idx];
      pma_idx++;
    } else if (current_src_pma == current_src_batch) {
      if (current_dest_pma < current_dest_batch || current_dest_pma == SENT_VAL) {
        //printf("2. %u, %u\n", current_src_pma, dests[pma_idx]);
        new_dests[result_idx] = dests[pma_idx];
        new_vals[result_idx] = vals[pma_idx];
        pma_idx++;
      } else if (current_dest_pma == current_dest_batch) {
        //printf("3. %u, %u\n", current_src_pma, current_dest_batch);
        new_dests[result_idx] = current_dest_batch;
        new_vals[result_idx] = current_val_batch;
        pma_idx++;
        es_idx++;
      } else {
        //printf("4. %u, %u\n", current_src_batch, current_dest_batch);
        new_dests[result_idx] = current_dest_batch;
        new_vals[result_idx] = current_val_batch;
        es_idx++;
      }
    } else { // current_src_pma > current_src_batch
      //printf("5. %u, %u\n", current_src_batch, current_dest_batch);
      new_dests[result_idx] = current_dest_batch;
      new_vals[result_idx] = current_val_batch;
      es_idx++;
    }
    result_idx++;
  }
  return result_idx - output_start;
}

// large merge
void PMA::add_edge_batch_update_no_val_parallel(pair_uint *es, uint64_t edge_count) {
  //printf("starting add_edge_batch_update_no_val_parallel\n");
  uint64_t task_id = __sync_fetch_and_add(&next_task_id, 2*edge_count);
  node_lock.lock_shared(0);

  //printf("\nlarge merge in batch\n");
  grab_all_locks(task_id, true, GENERAL);
  //std::sort(es, es + edge_count, sort_helper_x);
  //printf("starting sort\n");
  integerSort_y((pair_els*)es, edge_count, get_n());
  integerSort_x((pair_els*)es, edge_count, get_n());
  //printf("finished sort\n");
  //for (int i = 0; i < edge_count; i++) {
  //  printf("%u, %u\n", es[i].x, es[i].y);
  //}

  uint32_t num_workers = getWorkers()*5;
  //printf("edges.N = %lu, edge_count = %lu\n",edges.N, edge_count);
  uint64_t max_elts = edges.N + edge_count;
  //uint64_t new_size = 1 << bsr_word(max_elts);
  uint64_t new_size = max_elts - 1;
  new_size |= new_size >> 1;
  new_size |= new_size >> 2;
  new_size |= new_size >> 4;
  new_size |= new_size >> 8;
  new_size |= new_size >> 16;
  new_size |= new_size >> 32;
  new_size+=1;
  if (new_size < max_elts) { new_size *= 2; }

  uint32_t* dests = (uint32_t*) edges.dests;
  uint32_t* vals = (uint32_t*) edges.vals;
  uint32_t* new_vals = (uint32_t*) malloc(new_size * sizeof(uint32_t));
  if (new_vals == NULL) {
    printf("bad alloc for new_vals, trying to alloc something of size %lu\n", new_size);
    exit(-1);
  }
  uint32_t* new_dests = (uint32_t*) malloc(new_size * sizeof(uint32_t));
  if (new_dests == NULL) {
    printf("bad alloc for new_dests, trying to alloc something of size %lu\n", new_size);
    exit(-1);
  }

  std::vector<uint64_t> counts(num_workers+1);
  counts[0] = 0;
  std::vector<uint64_t> output_starts(num_workers);
  uint32_t increment_batch = edge_count / num_workers;
  //printf("starting merge\n");
  parallel_for(uint32_t i = 0; i < num_workers; i++) {

    uint64_t batch_start = increment_batch * i;
    uint64_t batch_end = std::min(batch_start + increment_batch, (uint64_t) edge_count);
    if (i == num_workers-1) {
      batch_end = edge_count;
    }

    uint64_t pma_start = 0;
    if (i > 0) {
      pma_start = binary_search(&edges, es[batch_start].y, 0, nodes[es[batch_start].x].beginning + 1,  nodes[es[batch_start].x].end);
    }
    uint64_t pma_end = edges.N;
    if (i < num_workers - 1) {
      pma_end = binary_search(&edges, es[batch_end].y, 0, nodes[es[batch_end].x].beginning + 1, nodes[es[batch_end].x].end);
    }

    uint64_t output_start = batch_start + pma_start;
    output_starts[i] = output_start;
    counts[i+1] = do_merge(new_dests, new_vals, output_start, es, batch_start, batch_end, dests, vals,  pma_start, pma_end, edges.loglogN);
  }

  // prefix sum
  for(uint32_t i = 1; i < num_workers + 1; i++) {
    counts[i] += counts[i-1];
  }

  // should be nnz
  uint64_t result_idx = counts[num_workers];

  uint64_t final_size = new_size;
  double cur_density = (double) result_idx / final_size;
  while (cur_density < .25) {
    final_size /= 2;
    cur_density = (double) result_idx / final_size;
  }
  while (cur_density > .75) {
    final_size *= 2;
    cur_density = (double) result_idx / final_size;
  }
  assert(cur_density >= .25 && cur_density <= .75);
  // printf("new size %lu, final_size %lu, nnz %lu\n", new_size, final_size, result_idx);
  edges.N = final_size;
  edges.loglogN = bsr_word(bsr_word(final_size) + 1);
  edges.logN = 1 << edges.loglogN;
  edges.mask_for_leaf = ~ (edges.logN - 1);
  edges.density_limit = ((double) edges.logN - 1) / edges.logN;
  edges.H = bsr_word(new_size / edges.logN);
  for(unsigned int i = 0; i <= edges.H; i++) {
    lower_density_bound[i] = density_bound(&edges, i).x;
    upper_density_bound[i] = density_bound(&edges, i).y;
  }

  // allocate space for dense
  free((void*)edges.dests);
  free((void*)edges.vals);
  uint32_t* dense_vals = (uint32_t*) malloc(result_idx * sizeof(uint32_t));

  if (dense_vals == NULL) {
    printf("bad alloc for dense_vals\n");
    exit(-1);
  }

  uint32_t* dense_dests = (uint32_t*) malloc(result_idx * sizeof(uint32_t));
  if (dense_dests == NULL) {
    printf("bad alloc for dense_dests\n");
    exit(-1);
  }

  // copy into dense array
  //printf("starting make dense\n");
  parallel_for(uint32_t i = 0; i < num_workers; i++) {
    uint32_t start_output_idx = counts[i];
    uint32_t start_input_idx = output_starts[i];
    uint32_t end_input_idx = start_input_idx + (counts[i+1] - counts[i]);
    for(uint32_t j = start_input_idx; j < end_input_idx; j++) {
      if (new_dests[j] != NULL_VAL) {
        dense_vals[start_output_idx] = new_vals[j];
        dense_dests[start_output_idx] = new_dests[j];
        start_output_idx++;
      }
    }
  }

  free((void*)new_dests);
  free((void*)new_vals);

  uint32_t* final_vals = (uint32_t*) malloc(final_size * sizeof(uint32_t));
  uint32_t* final_dests = (uint32_t*) malloc(final_size * sizeof(uint32_t));

  edges.dests = final_dests;
  edges.vals = final_vals;
  uint64_t num_leaves = final_size >> edges.loglogN;
  uint64_t count_per_leaf = result_idx / num_leaves;
  uint64_t extra = result_idx % num_leaves;
  //printf("num_leaves = %lu, count_per_leaf = %lu, extra = %lu\n", num_leaves, count_per_leaf, extra);
  //printf("starting redistribute\n");
  parallel_for(uint64_t i = 0; i < num_leaves; i++) {
    uint64_t count_for_leaf = count_per_leaf + (i < extra);
    uint64_t out_idx_start = i << edges.loglogN;
    uint64_t in_idx_start = count_per_leaf * i + min(i, extra);
    uint64_t out_idx = out_idx_start;
    //printf("%lu, %lu, %lu, %lu\n", count_for_leaf, out_idx_start, in_idx_start, out_idx);
    for(uint64_t j = in_idx_start; j < in_idx_start + count_for_leaf; j++) {
      if (j > 0) {
        if (dense_vals[j] == dense_vals[j-1] && dense_dests[j] == dense_dests[j-1]) {
          //printf("skipping duplicate edge %u, %u\n", dense_vals[j], dense_dests[j]);
          continue;
        }
      }
      final_vals[out_idx] = dense_vals[j];
      final_dests[out_idx] = dense_dests[j];
      if (final_dests[out_idx] == SENT_VAL) {
        uint64_t node_idx = final_vals[out_idx];
        fix_sentinel(node_idx, out_idx);
      }
      out_idx++;
    }
    for (uint64_t j = out_idx; j < (i+1) << edges.loglogN; j++) {
      final_vals[j] = 0;
    }
    for (uint64_t j = out_idx; j < (i+1) << edges.loglogN; j++) {
      final_dests[j] = NULL_VAL;
    }
  }
  //printf("finished redistribute\n");

  free((void*)dense_dests);
  free((void*)dense_vals);
  //print_array(1);
  //check_sentinals();
  release_all_locks(task_id, true, GENERAL);
  node_lock.unlock_shared(0);
  //printf("after batch insert\n");
  //print_array(0);
  //printf("finished add_edge_batch_update_no_val_parallel");
}

// batch delete exclude merge
uint64_t do_exclude(uint32_t* new_dests, uint32_t* new_vals, uint64_t
output_start, pair_uint *es, uint64_t batch_start, uint64_t batch_end,
uint32_t* dests, uint32_t* vals, uint64_t pma_start, uint64_t pma_end, uint32_t
loglogN, uint64_t num_nodes)
{
  // dests = pma_dests, vals = pma_vals

  uint64_t es_idx = batch_start;
  uint64_t pma_idx = pma_start;
  uint64_t result_idx = output_start;
  // TODO: look around (could be a binary search if its too slow)
  uint32_t current_src_pma = 0;
  while(pma_idx > 0 && dests[pma_idx] != SENT_VAL) {
    pma_idx--;
  }
  if (pma_idx != 0) {
    current_src_pma = vals[pma_idx];
  }

  pma_idx = pma_start;
  if (pma_idx == 0) {
    new_dests[0] = dests[0];
    new_vals[0] = vals[0];
    result_idx++;
    pma_idx++;
  }

  // merge
  while(es_idx < batch_end || pma_idx < pma_end) {
    // merge the rest of the PMA in
    if (es_idx == batch_end) {
      while(pma_idx < pma_end) {
        new_dests[result_idx] = dests[pma_idx];
        new_vals[result_idx] = vals[pma_idx];
        result_idx+=(dests[pma_idx] != NULL_VAL);
        pma_idx++;
      }
      break;
    }

    uint32_t current_src_batch = es[es_idx].x;
    uint32_t current_dest_batch = es[es_idx].y;
    // skip over duplicates
    while(es_idx < batch_end && current_src_batch == es[es_idx+1].x && current_dest_batch == es[es_idx+1].y) {
      es_idx++;
    }
    assert(es_idx < batch_end);
    assert(pma_idx < pma_end);

    if (pma_idx < pma_end && dests[pma_idx] == NULL_VAL) {
      pma_idx = next_leaf(pma_idx, loglogN);
      uint32_t logN = 1 << loglogN;
      if (current_src_pma < num_nodes) {
        while (dests[pma_idx] == NULL_VAL) {
          pma_idx += logN;
        }
      } else {
        while (pma_idx < pma_end && dests[pma_idx] == NULL_VAL) {
          pma_idx += logN;
        }
      }
    }
    // at the last node
    if (pma_idx == pma_end) {
      break;
    }

    uint32_t current_dest_pma = dests[pma_idx];
    assert(current_dest_pma != NULL_VAL);
    // check if pma elt is a sentinel
    if (current_dest_pma == SENT_VAL) {
      // TODO: remove after correctness tests
      /*if (!(current_src_pma == 0 || current_src_pma == vals[pma_idx] - 1 || current_src_pma == vals[pm
        printf("something wrong with pma sentinals giving us a bad pma source, %u, %u\n", current_src_pm
      }*/
      current_src_pma = vals[pma_idx];
    }
    if (current_src_pma < current_src_batch) {
      new_dests[result_idx] = dests[pma_idx];
      new_vals[result_idx] = vals[pma_idx];
      pma_idx++;
      result_idx++;
    } else if (current_src_pma == current_src_batch) {
      if (current_dest_pma < current_dest_batch || current_dest_pma == SENT_VAL) {
        new_dests[result_idx] = dests[pma_idx];
        new_vals[result_idx] = vals[pma_idx];
        pma_idx++;
        result_idx++;
      } else if (current_dest_pma == current_dest_batch) {
        pma_idx++;
        es_idx++;
      } else {
        es_idx++;
      }
    } else { // current_src_pma > current_src_batch
      es_idx++;
    }
  }
  return result_idx - output_start;
}

// do merge or not based on size of batch
void PMA::add_edge_batch_wrapper(pair_uint *es, uint64_t edge_count, int64_t threshold) {
  // default
  if (threshold < 0) {
    if (edges.N > UINT32_MAX) {
      uint64_t top = edges.N >> 32;
      threshold = edges.N / (32+bsr_word(top));
    } else if (edges.N == UINT32_MAX) {
      threshold = edges.N / 32;
    } else {
      threshold = edges.N / bsr_word(edges.N);
    }
  }
  // edge_count < num_cells / log(num_cells)
  //printf("threshold %lu, batch size %lu\n", threshold, edge_count);

  // actual inserts into PMA
  if (true || edge_count < (uint64_t) threshold*2) {
    node_lock.lock_shared(0);
    uint64_t task_id = __sync_fetch_and_add(&next_task_id, 2*edge_count);
    cilk_for (uint32_t i = 0; i < edge_count; i++) {
      uint32_t src = es[i].x;
      uint32_t dest = es[i].y;
      add_edge_update_fast(src,dest,1,task_id+2*i);
    }
    node_lock.unlock_shared(0);
  } else { // large merge
    printf("using large merge\n");
    add_edge_batch_update_no_val_parallel(es, edge_count);
  }
}

// large merge
void PMA::remove_edge_batch_update_no_val_parallel(pair_uint *es, uint64_t edge_count) {
  uint64_t task_id = __sync_fetch_and_add(&next_task_id, 2*edge_count);
  node_lock.lock_shared(0);

  //printf("\nlarge merge in batch\n");
  grab_all_locks(task_id, true, GENERAL);
  //std::sort(es, es + edge_count, sort_helper_x);
  //printf("starting sort\n");
  integerSort_y((pair_els*)es, edge_count, get_n());
  integerSort_x((pair_els*)es, edge_count, get_n());
  //printf("finished sort\n");
  //for (int i = 0; i < edge_count; i++) {
  //  printf("%u, %u\n", es[i].x, es[i].y);
  //}

  uint32_t num_workers = getWorkers()*5;
  uint64_t new_size = edges.N;

  uint32_t* dests = (uint32_t*) edges.dests;
  uint32_t* vals = (uint32_t*) edges.vals;
  uint32_t* new_vals = (uint32_t*) malloc(new_size * sizeof(uint32_t));
  if (new_vals == NULL) {
    printf("bad alloc for new_vals\n");
    exit(-1);
  }
  uint32_t* new_dests = (uint32_t*) malloc(new_size * sizeof(uint32_t));
  if (new_dests == NULL) {
    printf("bad alloc for new_dests\n");
    exit(-1);
  }
  std::vector<uint64_t> counts(num_workers+1);
  counts[0] = 0;
  std::vector<uint64_t> output_starts(num_workers);
  uint32_t increment_batch = edge_count / num_workers;
  //printf("starting merge\n");
  parallel_for(uint32_t i = 0; i < num_workers; i++) {

    uint64_t batch_start = increment_batch * i;
    uint64_t batch_end = std::min(batch_start + increment_batch, (uint64_t) edge_count);
    if (i == num_workers-1) {
      batch_end = edge_count;
    }

    uint64_t pma_start = 0;
    if (i > 0) {
      pma_start = binary_search(&edges, es[batch_start].y, 0, nodes[es[batch_start].x].beginning + 1, nodes[es[batch_start].x].end);
    }
    uint64_t pma_end = edges.N;
    if (i < num_workers - 1) {
      pma_end = binary_search(&edges, es[batch_end].y, 0, nodes[es[batch_end].x].beginning + 1, nodes[es[batch_end].x].end);
    }

    uint64_t output_start = pma_start;
    output_starts[i] = output_start;
    counts[i+1] = do_exclude(new_dests, new_vals, output_start, es, batch_start, batch_end, dests, vals, pma_start, pma_end, edges.loglogN, get_n());
  }
  //printf("finished merge\n");
  // print_array(0);
  // for (auto it : output_starts) {
  //   printf("%lu\n", it);
  // }

  // prefix sum
  for(uint32_t i = 1; i < num_workers + 1; i++) {
    counts[i] += counts[i-1];
  }

  // should be nnz
  uint64_t result_idx = counts[num_workers];

  uint64_t final_size = new_size;
  double cur_density = (double) result_idx / final_size;
  while (cur_density < .25) {
    final_size /= 2;
    cur_density = (double) result_idx / final_size;
  }
  while (cur_density > .75) {
    final_size *= 2;
    cur_density = (double) result_idx / final_size;
  }
  assert(cur_density >= .25 && cur_density <= .75);
  // printf("new size %lu, final_size %lu, nnz %lu\n", new_size, final_size, result_idx);
  edges.N = final_size;
  edges.loglogN = bsr_word(bsr_word(final_size) + 1);
  edges.logN = 1 << edges.loglogN;
  edges.mask_for_leaf = ~ (edges.logN - 1);
  edges.density_limit = ((double) edges.logN - 1) / edges.logN;
  edges.H = bsr_word(new_size / edges.logN);
  for(unsigned int i = 0; i <= edges.H; i++) {
    lower_density_bound[i] = density_bound(&edges, i).x;
    upper_density_bound[i] = density_bound(&edges, i).y;
  }
  // allocate space for dense
  free((void*)edges.dests);
  free((void*)edges.vals);
  uint32_t* dense_vals = (uint32_t*) malloc(result_idx * sizeof(uint32_t));

  if (dense_vals == NULL) {
    printf("bad alloc for dense_vals\n");
    exit(-1);
  }

  uint32_t* dense_dests = (uint32_t*) malloc(result_idx * sizeof(uint32_t));
  if (dense_dests == NULL) {
    printf("bad alloc for dense_dests\n");
    exit(-1);
  }

  // copy into dense array
  //printf("starting make dense\n");
  parallel_for(uint32_t i = 0; i < num_workers; i++) {
    uint32_t start_output_idx = counts[i];
    uint32_t start_input_idx = output_starts[i];
    uint32_t end_input_idx = start_input_idx + (counts[i+1] - counts[i]);
    for(uint32_t j = start_input_idx; j < end_input_idx; j++) {
      if (new_dests[j] != NULL_VAL) {
        dense_vals[start_output_idx] = new_vals[j];
        dense_dests[start_output_idx] = new_dests[j];
        start_output_idx++;
      }
    }
  }

  free((void*)new_dests);
  free((void*)new_vals);

  uint32_t* final_vals = (uint32_t*) malloc(final_size * sizeof(uint32_t));
  uint32_t* final_dests = (uint32_t*) malloc(final_size * sizeof(uint32_t));

  edges.dests = final_dests;
  edges.vals = final_vals;
  uint64_t num_leaves = final_size >> edges.loglogN;
  uint64_t count_per_leaf = result_idx / num_leaves;
  uint64_t extra = result_idx % num_leaves;
  //printf("num_leaves = %lu, count_per_leaf = %lu, extra = %lu\n", num_leaves, count_per_leaf, extra);
  //printf("starting redistribute\n");
  parallel_for(uint64_t i = 0; i < num_leaves; i++) {
    uint64_t count_for_leaf = count_per_leaf + (i < extra);
    uint64_t out_idx_start = i << edges.loglogN;
    uint64_t in_idx_start = count_per_leaf * i + min(i, extra);
    uint64_t out_idx = out_idx_start;
    //printf("%lu, %lu, %lu, %lu\n", count_for_leaf, out_idx_start, in_idx_start, out_idx);
    for(uint64_t j = in_idx_start; j < in_idx_start + count_for_leaf; j++) {
      final_vals[out_idx_start] = dense_vals[j];
      out_idx_start++;
    }
    for (uint64_t j = out_idx_start; j < (i+1) << edges.loglogN; j++) {
      final_vals[j] = 0;
    }
    for(uint64_t j = in_idx_start; j < in_idx_start + count_for_leaf; j++) {
      final_dests[out_idx] = dense_dests[j];
      if (final_dests[out_idx] == SENT_VAL) {
        uint64_t node_idx = final_vals[out_idx];
        fix_sentinel(node_idx, out_idx);
      }
      out_idx++;
    }
    for (uint64_t j = out_idx; j < (i+1) << edges.loglogN; j++) {
      final_dests[j] = NULL_VAL;
    }
  }

  free((void*)dense_dests);
  free((void*)dense_vals);
  //print_array(1);
  //check_sentinals();
  release_all_locks(task_id, true, GENERAL);
  node_lock.unlock_shared(0);
  //printf("after batch insert\n");
  //print_array(0);
}

// do merge or not based on size of batch
void PMA::remove_edge_batch_wrapper(pair_uint *es, uint64_t edge_count, int64_t threshold) {
  // default
  if (threshold < 0) {
    if (edges.N > UINT32_MAX) {
      uint64_t top = edges.N >> 32;
      threshold = edges.N / (32+bsr_word(top));
    } else if (edges.N == UINT32_MAX) {
      threshold = edges.N / 32;
    } else {
      threshold = edges.N / bsr_word(edges.N);
    }
  }

  //printf("threshold%lu, batch size %lu\n", threshold, edge_count);

  // actual inserts into PMA
  if (true || edge_count < (uint64_t) threshold*2) {
    node_lock.lock_shared(0);
    uint64_t task_id = __sync_fetch_and_add(&next_task_id, 2*edge_count);
    cilk_for (uint32_t i = 0; i < edge_count; i++) {
      uint32_t src = es[i].x;
      uint32_t dest = es[i].y;
      remove_edge_fast(src,dest,task_id+2*i);
    }
    node_lock.unlock_shared(0);
  } else { // large merge
    printf("using large merge\n");
    remove_edge_batch_update_no_val_parallel(es, edge_count);
  }
}

bool PMA::remove_edge(uint32_t src, uint32_t dest) {
  //printf(" trying to remove src = %d, dest = %d\n", src,dest);
  uint64_t task_id = __sync_fetch_and_add(&next_task_id, 2);
  #if ENABLE_PMA_LOCK == 1
  assert(check_no_locks_for_me(task_id));

  node_lock.lock_shared(task_id);
  pair_int held_locks = grab_locks_for_leaf_with_resets(task_id, src);
  assert(nodes[src].lock.i_own_lock(task_id));
  #endif
  node_t node = nodes[src];

  
  // edges.list_lock.lock_shared();
  //uint32_t node_begin = node.beginning;
  //uint32_t node_end = node.end;
  // printf("looking in region %u, %u by worker %lu\n", node.beginning + 1, node.end, get_worker_num());
  uint_t loc_to_remove =
      binary_search(&edges, dest, 0, node.beginning + 1, node.end);
  if (edges.dests[loc_to_remove] != dest) {
    //if the edge isn't there
    //printf("not removed\n");
    #if ENABLE_PMA_LOCK == 1
    release_locks_in_range(task_id, held_locks);
    node_lock.unlock_shared(task_id);
    #endif
    return false;
  }
  //printf("removed\n");
  __sync_fetch_and_add(&nodes[src].num_neighbors,-1);
  #if ENABLE_PMA_LOCK == 1
  pair_int needed_locks = which_locks_in_range(find_leaf(&edges, loc_to_remove), edges.logN, src);
  for (uint32_t i = held_locks.x; i < needed_locks.x; i++) {
    nodes[i].lock.unlock(task_id, GENERAL);
    held_locks.x = i+1;
  }
  for (uint32_t i = needed_locks.y+1; i <= held_locks.y; i++) {
    nodes[i].lock.unlock(task_id, GENERAL);
    
  }
  if (needed_locks.y+1 <= held_locks.y) {
    held_locks.y = needed_locks.y;
  }
  // print_array();
  //printf("loc_to_add: %u by worker %lu\n", loc_to_add, get_worker_num());
  // printf("src: %d, dest: %d by worker %lu\n", src, dest, get_worker_num());
  // print_array();

  
  assert(nodes[src].lock.i_own_lock(task_id));
  assert(check_every_lock_in_leaf(task_id, loc_to_remove));
  /*
  if (!check_no_full_leaves(&edges, find_leaf(&edges, loc_to_add), edges.logN)) {
    printf("worker %lu is trying to insert into a full leaf\n", get_worker_num());
    print_array(get_worker_num());
  }
  */
  assert(check_no_full_leaves(&edges, find_leaf(&edges, loc_to_remove), edges.logN));
  #endif

  remove(task_id, loc_to_remove, dest, src
    #if ENABLE_PMA_LOCK == 1
    , held_locks
    #endif
    );
  //printf("worker %lu done\n", get_worker_num());
  // edges.list_lock.unlock_shared();
  //print_array();
  #if ENABLE_PMA_LOCK == 1
  node_lock.unlock_shared(task_id);
  assert(check_no_locks_for_me(task_id));
  #endif
  return true;
}

void PMA::remove_edge_fast(uint32_t src, uint32_t dest, uint64_t task_id) {
  //printf(" trying to remove src = %d, dest = %d\n", src,dest);
  #if ENABLE_PMA_LOCK == 1
  pair_int held_locks = grab_locks_for_leaf_with_resets(task_id, src);
  #endif
  node_t node = nodes[src];

  
  // edges.list_lock.lock_shared();
  //uint32_t node_begin = node.beginning;
  //uint32_t node_end = node.end;
  // printf("looking in region %u, %u by worker %lu\n", node.beginning + 1, node.end, get_worker_num());
  uint_t loc_to_remove =
      binary_search(&edges, dest, 0, node.beginning + 1, node.end);
  if (edges.dests[loc_to_remove] != dest) {
    //if the edge isn't there
    //printf("not removed\n");
    #if ENABLE_PMA_LOCK == 1
    release_locks_in_range(task_id, held_locks);
    #endif
    return;
  }
  //printf("removed\n");
  __sync_fetch_and_add(&nodes[src].num_neighbors,-1);
  #if ENABLE_PMA_LOCK == 1
  pair_int needed_locks = which_locks_in_range(find_leaf(&edges, loc_to_remove), edges.logN, src);
  for (uint32_t i = held_locks.x; i < needed_locks.x; i++) {
    nodes[i].lock.unlock(task_id, GENERAL);
    held_locks.x = i+1;
  }
  for (uint32_t i = needed_locks.y+1; i <= held_locks.y; i++) {
    nodes[i].lock.unlock(task_id, GENERAL);
    
  }
  if (needed_locks.y+1 <= held_locks.y) {
    held_locks.y = needed_locks.y;
  }
  // print_array();
  //printf("loc_to_add: %u by worker %lu\n", loc_to_add, get_worker_num());
  // printf("src: %d, dest: %d by worker %lu\n", src, dest, get_worker_num());
  // print_array();

  
  assert(nodes[src].lock.i_own_lock(task_id));
  /*
  if (!check_no_full_leaves(&edges, find_leaf(&edges, loc_to_add), edges.logN)) {
    printf("worker %lu is trying to insert into a full leaf\n", get_worker_num());
    print_array(get_worker_num());
  }
  */
  #endif
  assert(check_no_full_leaves(&edges, find_leaf(&edges, loc_to_remove), edges.logN));
  remove(task_id, loc_to_remove, dest, src
   #if ENABLE_PMA_LOCK == 1
    , held_locks
    #endif
    );
  //printf("worker %lu done\n", get_worker_num());
  // edges.list_lock.unlock_shared();
  //print_array();
}

void PMA::remove_edge_batch(uint32_t *srcs, uint32_t *dests, uint_t edge_count) {
  uint64_t task_id = __sync_fetch_and_add(&next_task_id, 2*edge_count);
#if ENABLE_PMA_LOCK == 1
  node_lock.lock_shared(0);
  parallel_for (uint_t i = 0; i < edge_count; i++) {
#else
  for (uint_t i = 0; i < edge_count; i++) {
#endif
    uint32_t src = srcs[i];
    uint32_t dest = dests[i];
    remove_edge_fast(src,dest,task_id+2*i);
  }
#if ENABLE_PMA_LOCK == 1
  node_lock.unlock_shared(0);
#endif
}



PMA::PMA(uint32_t init_n) {
  next_task_id = 1;
  //making sure logN is at least 4
  edges.N = max(2UL << bsr_word(init_n*2), 16UL);
  // printf("%d\n", bsf_word(list->N));
  edges.loglogN = bsr_word(bsr_word(edges.N) + 1);
  edges.logN = (1 << edges.loglogN);
  edges.mask_for_leaf = ~(edges.logN - 1);
  assert(edges.logN > 0);
  edges.density_limit = ((double) edges.logN - 1)/edges.logN;
  edges.H = bsr_word(edges.N / edges.logN);
  for (uint32_t i = 0; i <= edges.H; i++) {
    upper_density_bound[i] = density_bound(&edges, i).y;
    lower_density_bound[i] = density_bound(&edges, i).x;
  }

  // printf("N = %d, logN = %d, loglogN = %d, H = %d\n", list->N, list->logN,
  // list->loglogN, list->H);
  
  // edges.list_lock.name("listLock");
#if ENABLE_PMA_LOCK == 1
  node_lock.name("nodeLock");  
#endif
  //edges.vals = (uint32_t *)malloc(edges.N * sizeof(*(edges.vals)));
  //edges.dests = (uint32_t *)malloc(edges.N * sizeof(*(edges.dests)));
  edges.vals = (uint32_t *)aligned_alloc(32, edges.N * sizeof(*(edges.vals)));
  edges.dests = (uint32_t *)aligned_alloc(32, edges.N * sizeof(*(edges.dests)));
  parallel_for (uint_t i = 0; i < edges.N; i++) {
    edges.vals[i] = 0;
    edges.dests[i] = NULL_VAL;
  }
  //TODO might be an issue if we grow it one at a time and let nodes be moved during operation
  nodes.reserve(init_n);
  for (uint32_t i = 0; i < init_n; i++) {
    add_node();
  }
}

PMA::PMA(PMA &other) {
  nodes = other.nodes;
  next_task_id = 1;
  edges.N = other.edges.N;
  edges.loglogN = other.edges.loglogN;
  edges.logN = other.edges.logN;
  edges.mask_for_leaf = other.edges.mask_for_leaf;
  edges.density_limit = other.edges.density_limit;
  edges.H = other.edges.H;
  for (uint32_t i = 0; i <= edges.H; i++) {
    upper_density_bound[i] = other.upper_density_bound[i];
    lower_density_bound[i] = other.lower_density_bound[i];

  }
#if ENABLE_PMA_LOCK == 1
  node_lock.name("nodeLock");
#endif
  //edges.vals = (uint32_t *)malloc(edges.N * sizeof(*(edges.vals)));
  edges.vals = (uint32_t *)aligned_alloc(32, edges.N * sizeof(*(edges.vals)));
  memcpy(__builtin_assume_aligned((void *)edges.vals,32), __builtin_assume_aligned((void *)other.edges.vals, 32), edges.N * sizeof(*(edges.vals)));
  //edges.dests = (uint32_t *)malloc(edges.N * sizeof(*(edges.dests)));
  edges.dests = (uint32_t *)aligned_alloc(32, edges.N * sizeof(*(edges.dests)));
  memcpy(__builtin_assume_aligned((void *)edges.dests, 32), __builtin_assume_aligned((void *)other.edges.dests, 32), edges.N * sizeof(*(edges.dests)));
}

//TODO free lock array
PMA::~PMA() { 
  free((void*)edges.vals);
  free((void*)edges.dests);
}
