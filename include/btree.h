/*
 * ============================================================================
 *
 *       Filename:  btree.h
 *
 *         Author:  Prashant Pandey (), ppandey2@cs.cmu.edu
 *   Organization:  Carnegie Mellon University
 *
 * ============================================================================
 */

#ifndef _BTREE_H_
#define _BTREE_H_

#include <string.h>

#include <utility>
#include <iostream>
#include <stack>

// #define WEIGHTED 0

#define MIN_KEYS 64
#define MAX_KEYS (2*MIN_KEYS-1)
#define MAX_CHILDREN (2*MIN_KEYS)

template<class T, class W>
class BTreeNode {
  public:
    BTreeNode(bool is_leaf);

    const BTreeNode<T, W> *find(T k) const;
#if WEIGHTED
    W get_val(T k) const;
#endif
#if WEIGHTED
    bool insertNonFull(T k, W w);
#else
    bool insertNonFull(T k);
#endif

    bool remove(T k);
    uint32_t findKey(uint32_t k);
    void removeFromLeaf(uint32_t idx);
    void removeFromNonLeaf(uint32_t idx);
#if WEIGHTED
    std::pair<T, W> getPred(uint32_t idx);
    std::pair<T, W> getSucc(uint32_t idx);
#else
    T getPred(uint32_t idx);
    T getSucc(uint32_t idx);
#endif
    void fill(uint32_t idx);
    void borrowFromPrev(uint32_t idx);
    void borrowFromNext(uint32_t idx);
    void merge(uint32_t idx);

    // The Child c must be full when this function is called
    void splitChild(uint32_t i, BTreeNode<T, W> *c);
    void traverse() const;
    uint64_t sum() const;

    uint32_t get_num_nodes() const;

    template <class F> void map_node(F &f) {
      for (uint64_t i = 0; i < num_keys; i++) {
#if WEIGHTED
        f.update(keys[i], weights[i]);
#else
        f.update(keys[i]);
#endif
        }
      }

    template <class F> void map_tree(F &f) {
      map_node(f);
      if (!is_leaf) {
        //parallel_for (uint32_t i = 0; i < num_keys+1; i++) {
        for (uint32_t i = 0; i < num_keys+1; i++) {
          if (children[i] != nullptr)
            children[i]->map_tree(f);
        }
      }
    }

    class NodeIterator {
      public:
        NodeIterator() {};
        NodeIterator(const BTreeNode<T, W> *node);
#if WEIGHTED
        std::pair<T, W> operator*(void) const;
#else
        T operator*(void) const;
#endif
        void operator++(void);
        bool done(void) const;

      private:
        std::stack<std::pair<const BTreeNode<T, W>*, uint32_t>> stack;
        const BTreeNode<T, W> *cur_node;
        uint32_t cur_idx;
    };

    NodeIterator begin(void) const { return NodeIterator(this); }

    template<typename C, typename D>
    friend class BTree;

  private:
    T keys[MAX_KEYS];
#if WEIGHTED
    W weights[MAX_KEYS];
#endif
    BTreeNode<T, W> *children[MAX_CHILDREN];
    uint32_t num_keys;
    bool is_leaf;
};

template <class T, class W>
class BTree {
  public:
    BTree() : root(nullptr) {}

#if WEIGHTED
    W get_val(T k) const {
      return (root == nullptr) ? 0 : root->get_val(k);
    }
#endif
    const BTreeNode<T, W>* find(T k) const {
      return (root == nullptr) ? nullptr : root->find(k);
    }
#if WEIGHTED
    bool insert(T k, W w);
#else
    bool insert(T k);
#endif
    bool remove(T k);

    void traverse() const {
      if (root != nullptr) root->traverse();
    }

    uint64_t sum() const {
      if (root != nullptr)
        return root->sum();
      return 0;
    }

    uint32_t get_num_nodes() const {
      if (root != nullptr) return root->get_num_nodes();
      return 0;
    }

    uint64_t get_size(void) const {
      return get_num_nodes() * sizeof(BTreeNode<T, W>);
    }

    const BTreeNode<T, W>* get_root(void) const {
      return root;
    }

    template<class F> void map(F &f) {
      root->map_tree(f);
    }

    class Iterator {
      public:
        Iterator() {};
        Iterator(const BTree<T, W> *tree) {
          assert(tree->root != nullptr);
          it = tree->root->begin();
        }
#if WEIGHTED
        std::pair<T, W> operator*(void) const { return *it; }
#else
        T operator*(void) const { return *it; }
#endif
        void operator++(void) { ++it; }
        bool done(void) const { return it.done(); }

      private:
        typename BTreeNode<T, W>::NodeIterator it;
    };

    Iterator begin(void) const { return Iterator(this); }

  private:
    BTreeNode<T, W> *root; // Pointer to root node
};


template <class T, class W>
BTreeNode<T, W>::BTreeNode(bool is_leaf) : num_keys(0), is_leaf(is_leaf) {
  memset(keys, 0, MAX_KEYS*sizeof(keys[0]));
#if WEIGHTED
  memset(weights, 0, MAX_KEYS*sizeof(weights[0]));
#endif
  memset(children, 0, MAX_CHILDREN*sizeof(children[0]));
}

#if WEIGHTED
template <class T, class W>
W BTreeNode<T, W>::get_val(T k) const {
  uint32_t i = 0;
  while (i < num_keys && k > keys[i])
    i++;

  if (keys[i] == k)
    return weights[i];

  if (is_leaf)
    return 0;

  return children[i]->get_val(k);
}
#endif

template <class T, class W>
const BTreeNode<T, W> *BTreeNode<T, W>::find(T k) const {
  uint32_t i = 0;
  while (i < num_keys && k > keys[i])
    i++;

  if (i < num_keys && keys[i] == k)
    return this;

  if (is_leaf)
    return nullptr;

  return children[i]->find(k);
}

template <class T, class W>
void BTreeNode<T, W>::traverse() const {
  uint32_t i;
  for (i = 0; i < num_keys; i++)
  {
    // If this is not leaf, then before printing key[i],
    // traverse the subtree rooted with child C[i].
    if (!is_leaf)
      children[i]->traverse();
    std::cout << keys[i] << " ";
  }

  // Print the subtree rooted with last child
  if (!is_leaf)
    children[i]->traverse();
}

template <class T, class W>
uint64_t BTreeNode<T, W>::sum() const {
  uint32_t i;
  uint64_t count{0};
  for (i = 0; i < num_keys; i++)
  {
    // If this is not leaf, then before printing key[i],
    // traverse the subtree rooted with child C[i].
    if (!is_leaf)
      count += children[i]->sum();
    count += keys[i];  
  }

  // Print the subtree rooted with last child
  if (!is_leaf)
    count += children[i]->sum();
  return count;
}

template <class T, class W>
uint32_t BTreeNode<T, W>::get_num_nodes() const {
  uint32_t count{1};
  uint32_t i;
  for (i = 0; i < num_keys; i++) {
    if (!is_leaf)
      count += children[i]->get_num_nodes();
  }
  // Print the subtree rooted with last child
  if (!is_leaf)
    count += children[i]->get_num_nodes();
  return count;
}

template <class T, class W>
#if WEIGHTED
bool BTreeNode<T, W>::insertNonFull(T k, W w) {
#else
bool BTreeNode<T, W>::insertNonFull(T k) {
#endif
  // Initialize index as index of rightmost element
  uint32_t idx;
  for (idx = 0; idx < num_keys; idx++) {
    if (keys[idx] < k)
      continue;
    else if (k == keys[idx]) {
#if WEIGHTED
      weights[idx] = w;
#endif
      return false;
    }
    else
      break;
  }

  if (is_leaf) { // If this is a leaf node
    memmove(&keys[idx+1], &keys[idx], (MAX_KEYS-idx-1)*sizeof(keys[0]));
#if WEIGHTED
    memmove(&weights[idx+1], &weights[idx], (MAX_KEYS-idx-1)*sizeof(weights[0]));
#endif
    // Insert the new key at found location
    keys[idx] = k;
#if WEIGHTED
    weights[idx] = w;
#endif
    num_keys = num_keys+1;
    return true;
  } else { // If this node is not leaf
    if (children[idx]->num_keys == MAX_KEYS) {// See if the found child is full
      // If the child is full, then split it
      splitChild(idx, children[idx]);

      // After split, the middle key of children[idx] goes up and
      // children[idx] is splitted into two. See which of the two
      // is going to have the new key
      if (keys[idx] == k)
        return false;
      else if (keys[idx] < k)
        idx++;
    }
#if WEIGHTED
    return children[idx]->insertNonFull(k, w);
#else
    return children[idx]->insertNonFull(k);
#endif
  }
}

template <class T, class W>
void BTreeNode<T, W>::splitChild(uint32_t i, BTreeNode<T, W> *c) {
  // Create a new node which is going to store (MIN_KEYS-1) keys
  // of c
  BTreeNode<T, W> *z = new BTreeNode<T, W>(c->is_leaf);
  z->num_keys = MIN_KEYS - 1;

  // Copy the last (MIN_KEYS-1) keys of c to z
  memmove(&z->keys[0], &c->keys[MIN_KEYS], (MIN_KEYS-1)*sizeof(keys[0]));
#if WEIGHTED
  memmove(&z->weights[0], &c->weights[MIN_KEYS], (MIN_KEYS-1)*sizeof(weights[0]));
#endif

  // Copy the last t children of y to z
  if (!c->is_leaf)
    memmove(&z->children[0], &c->children[MIN_KEYS],
            (MAX_KEYS)*sizeof(children[0]));

  // Reduce the number of keys in y
  c->num_keys = MIN_KEYS - 1;

  // Since this node is going to have a new child,
  // create space of new child
  memcpy(&children[i+2], &children[i+1],
         (MAX_CHILDREN-i-2)*sizeof(children[0]));

  // Link the new child to this node
  children[i+1] = z;

  // A key of y will move to this node. Find the location of
  // new key and move all greater keys one space ahead
  memcpy(&keys[i+1], &keys[i], (MAX_KEYS-i-1)*sizeof(keys[0]));
#if WEIGHTED
  memcpy(&weights[i+1], &weights[i], (MAX_KEYS-i-1)*sizeof(weights[0]));
#endif

  // Copy the middle key of y to this node
  keys[i] = c->keys[MIN_KEYS-1];
#if WEIGHTED
  weights[i] = c->weights[MIN_KEYS-1];
#endif

  // Increment count of keys in this node
  num_keys = num_keys + 1;
}

template <class T, class W>
BTreeNode<T, W>::NodeIterator::NodeIterator(const BTreeNode<T, W> *node) {
  assert(node != nullptr);
  const BTreeNode<T, W> *cur = node;
  while (cur != nullptr) {
    stack.push(std::make_pair(cur, 0));
    cur = cur->children[0];
  }
  cur_node = stack.top().first;
  cur_idx = 0;
}

template <class T, class W>
#if WEIGHTED
std::pair<T, W> BTreeNode<T, W>::NodeIterator::operator*(void) const {
  return std::make_pair(cur_node->keys[cur_idx], cur_node->weights[cur_idx]);
#else
T BTreeNode<T, W>::NodeIterator::operator*(void) const {
  return cur_node->keys[cur_idx];
#endif
}

template <class T, class W>
void BTreeNode<T, W>::NodeIterator::operator++(void) {
  if (cur_node->is_leaf) {  // if we are traversing a leaf node
    if (cur_idx < cur_node->num_keys-1)
      cur_idx++;
    else {  // leaf node is done. Go to parent
      stack.pop();
      if (stack.size()) {
        cur_node = stack.top().first;
        cur_idx = stack.top().second;
        stack.top().second++;
      }
    }
  } else { // if we are traversing an internal node
    if (cur_idx == cur_node->num_keys-1) { // internal node is done. Go to parent
      stack.pop();
    }
    BTreeNode<T, W> *cur = cur_node->children[cur_idx + 1];
    while (cur != nullptr) {
      stack.push(std::make_pair(cur, 0));
      cur = cur->children[0];
    }
    cur_node = stack.top().first;
    cur_idx = 0;
  }
}

template <class T, class W>
bool BTreeNode<T, W>::NodeIterator::done(void) const {
  return !stack.size();
}

template <class T, class W>
#if WEIGHTED
bool BTree<T, W>::insert(T k, W w) {
#else
bool BTree<T, W>::insert(T k) {
#endif
  // If tree is empty
  if (root == nullptr)  {
    // Allocate memory for root
    root = new BTreeNode<T, W>(true);
    root->keys[0] = k;  // Insert key
#if WEIGHTED
    root->weights[0] = w;  // Insert weight
#endif
    root->num_keys = 1;  // Update number of keys in root
    return true;
  } else {  // If tree is not empty
    // If root is full, then tree grows in height
    if (root->num_keys == MAX_KEYS) {
      // Allocate memory for new root
      BTreeNode<T, W> *s = new BTreeNode<T, W>(false);

      // Make old root as child of new root
      s->children[0] = root;

      // Split the old root and move 1 key to the new root
      s->splitChild(0, root);

      // New root has two children now.  Decide which of the
      // two children is going to have new key
      int i = 0;
      if (s->keys[0] < k)
        i++;

      // Change root
      root = s;
#if WEIGHTED
      return s->children[i]->insertNonFull(k, w);
#else
      return s->children[i]->insertNonFull(k);
#endif
    }
    else  // If root is not full, call insertNonFull for root
#if WEIGHTED
      return root->insertNonFull(k, w);
#else
      return root->insertNonFull(k);
#endif
  }
}

// A utility function that returns the index of the first key that is
// greater than or equal to k
template <class T, class W>
uint32_t BTreeNode<T, W>::findKey(uint32_t k) {
    int idx=0;
    while (idx<num_keys && keys[idx] < k)
        ++idx;
    return idx;
}

template <class T, class W>
void BTreeNode<T, W>::removeFromLeaf(uint32_t idx) {
  // Move all the keys after the idx-th pos one place backward
  memcpy(&keys[idx], &keys[idx+1], (MAX_KEYS-idx-1)*sizeof(keys[0]));
#if WEIGHTED
  memcpy(&weights[idx], &weights[idx+1], (MAX_KEYS-idx-1)*sizeof(weights[0]));
#endif

  //for (uint32_t i=idx+1; i<num_keys; ++i) {
    //keys[i-1] = keys[i];
//#if WEIGHTED
    //weights[i-1] = weights[i];
//#endif
  //}
  // Reduce the count of keys
  num_keys--;
}

// A function to remove the idx-th key from this node - which is a non-leaf node
template <class T, class W>
void BTreeNode<T, W>::removeFromNonLeaf(uint32_t idx) {
    T k = keys[idx];

    // If the child that precedes k (C[idx]) has atleast t keys,
    // find the predecessor 'pred' of k in the subtree rooted at
    // C[idx]. Replace k by pred. Recursively delete pred
    // in C[idx]
    if (children[idx]->num_keys >= MIN_KEYS) {
        auto pred = getPred(idx);
#if WEIGHTED
        keys[idx] = pred.first;
        weights[idx] = pred.second;
        children[idx]->remove(pred.first);
#else
        keys[idx] = pred;
        children[idx]->remove(pred);
#endif
    }

    // If the child C[idx] has less that t keys, examine C[idx+1].
    // If C[idx+1] has atleast t keys, find the successor 'succ' of k in
    // the subtree rooted at C[idx+1]
    // Replace k by succ
    // Recursively delete succ in C[idx+1]
    else if  (children[idx+1]->num_keys >= MIN_KEYS) {
        auto succ = getSucc(idx);
#if WEIGHTED
        keys[idx] = succ.first;
        weights[idx] = succ.second;
        children[idx+1]->remove(succ.first);
#else
        keys[idx] = succ;
        children[idx+1]->remove(succ);
#endif
    }

    // If both C[idx] and C[idx+1] has less that t keys,merge k and all of C[idx+1]
    // into C[idx]
    // Now C[idx] contains 2t-1 keys
    // Free C[idx+1] and recursively delete k from C[idx]
    else {
        merge(idx);
        children[idx]->remove(k);
    }
}

// A function to get predecessor of keys[idx] 
template <class T, class W>
#if WEIGHTED
std::pair<T, W> BTreeNode<T, W>::getPred(uint32_t idx) {
#else
T BTreeNode<T, W>::getPred(uint32_t idx) {
#endif
  // Keep moving to the right most node until we reach a leaf 
  BTreeNode<T, W> *cur = children[idx]; 
  while (!cur->is_leaf) 
    cur = cur->children[cur->num_keys]; 

  // Return the last key of the leaf 
#if WEIGHTED
  return std::make_pair(cur->keys[cur->num_keys-1], cur->weights[cur->num_keys-1]);
#else
  return cur->keys[cur->num_keys-1]; 
#endif
}
  
template <class T, class W>
#if WEIGHTED
std::pair<T, W> BTreeNode<T, W>::getSucc(uint32_t idx) {
#else
T BTreeNode<T, W>::getSucc(uint32_t idx) { 
#endif
  // Keep moving the left most node starting from C[idx+1] until we reach a leaf 
  BTreeNode *cur = children[idx+1]; 
  while (!cur->is_leaf) 
    cur = cur->children[0]; 

  // Return the first key of the leaf 
#if WEIGHTED
  return std::make_pair(cur->keys[0], cur->weights[0]);
#else
  return cur->keys[0]; 
#endif

} 

// A function to fill child C[idx] which has less than t-1 keys 
template <class T, class W>
void BTreeNode<T, W>::fill(uint32_t idx) { 
  // If the previous child(C[idx-1]) has more than t-1 keys, borrow a key 
  // from that child 
  if (idx != 0 && children[idx-1]->num_keys >= MIN_KEYS) 
    borrowFromPrev(idx); 

  // If the next child(C[idx+1]) has more than t-1 keys, borrow a key 
  // from that child 
  else if (idx != num_keys && children[idx+1]->num_keys >= MIN_KEYS) 
    borrowFromNext(idx);

  // Merge C[idx] with its sibling 
  // If C[idx] is the last child, merge it with with its previous sibling 
  // Otherwise merge it with its next sibling 
  else { 
    if (idx != num_keys) 
      merge(idx); 
    else
      merge(idx-1); 
  } 
} 

// A function to borrow a key from C[idx-1] and insert it
// into C[idx]
template <class T, class W>
void BTreeNode<T, W>::borrowFromPrev(uint32_t idx) {
  BTreeNode<T, W> *child = children[idx];
  BTreeNode<T, W> *sibling = children[idx-1];

  // The last key from C[idx-1] goes up to the parent and key[idx-1]
  // from parent is inserted as the first key in C[idx]. Thus, the  loses
  // sibling one key and child gains one key

  // Moving all key in C[idx] one step ahead
  for (int32_t i = child->num_keys-1; i >= 0; --i) {
    child->keys[i+1] = child->keys[i];
#if WEIGHTED
    child->weights[i+1] = child->weights[i];
#endif
  }

  // If C[idx] is not a leaf, move all its child pointers one step ahead
  if (!child->is_leaf) {
    for(int32_t i = child->num_keys; i >= 0; --i)
      child->children[i+1] = child->children[i];
  }

  // Setting child's first key equal to keys[idx-1] from the current node
  child->keys[0] = keys[idx-1];
#if WEIGHTED
  child->weights[0] = weights[idx-1];
#endif

  // Moving sibling's last child as C[idx]'s first child
  if(!child->is_leaf)
    child->children[0] = sibling->children[sibling->num_keys];

  // Moving the key from the sibling to the parent
  // This reduces the number of keys in the sibling
  keys[idx-1] = sibling->keys[sibling->num_keys-1];
#if WEIGHTED
  weights[idx-1] = sibling->weights[sibling->num_keys-1];
#endif

  child->num_keys += 1;
  sibling->num_keys -= 1;
} 

// A function to borrow a key from the C[idx+1] and place
// it in C[idx]
template <class T, class W>
void BTreeNode<T, W>::borrowFromNext(uint32_t idx) {
  BTreeNode<T, W> *child = children[idx];
  BTreeNode<T, W> *sibling = children[idx+1];

  // keys[idx] is inserted as the last key in C[idx]
  child->keys[(child->num_keys)] = keys[idx];
#if WEIGHTED
  child->weights[(child->num_keys)] = weights[idx];
#endif

  // Sibling's first child is inserted as the last child
  // into C[idx]
  if (!(child->is_leaf))
    child->children[(child->num_keys)+1] = sibling->children[0];

  //The first key from sibling is inserted into keys[idx]
  keys[idx] = sibling->keys[0];
#if WEIGHTED
  weights[idx] = sibling->weights[0];
#endif

  // Moving all keys in sibling one step behind
  for (uint32_t i=1; i < sibling->num_keys; ++i) {
    sibling->keys[i-1] = sibling->keys[i];
#if WEIGHTED
    sibling->weights[i-1] = sibling->weights[i];
#endif
  }
  // Moving the child pointers one step behind
  if (!sibling->is_leaf) {
    for(uint32_t i=1; i <= sibling->num_keys; ++i)
      sibling->children[i-1] = sibling->children[i];
  }

  // Increasing and decreasing the key count of C[idx] and C[idx+1]
  // respectively
  child->num_keys += 1;
  sibling->num_keys -= 1;
} 

// A function to merge C[idx] with C[idx+1]
// C[idx+1] is freed after merging
template <class T, class W>
void BTreeNode<T, W>::merge(uint32_t idx) {
  BTreeNode *child = children[idx];
  BTreeNode *sibling = children[idx+1];

  // Pulling a key from the current node and inserting it into (t-1)th
  // position of C[idx]
  child->keys[MIN_KEYS-1] = keys[idx];
#if WEIGHTED
  child->weights[MIN_KEYS-1] = weights[idx];
#endif

  // Copying the keys from C[idx+1] to C[idx] at the end
  for (uint32_t i=0; i < sibling->num_keys; ++i) {
    child->keys[i+MIN_KEYS] = sibling->keys[i];
#if WEIGHTED
    child->weights[i+MIN_KEYS] = sibling->weights[i];
#endif
  }

  // Copying the child pointers from C[idx+1] to C[idx]
  if (!child->is_leaf) {
    for(uint32_t i=0; i <= sibling->num_keys; ++i)
      child->children[i+MIN_KEYS] = sibling->children[i];
  }

  // Moving all keys after idx in the current node one step before -
  // to fill the gap created by moving keys[idx] to C[idx]
  for (uint32_t i=idx+1; i<num_keys; ++i) {
    keys[i-1] = keys[i];
#if WEIGHTED
    weights[i-1] = weights[i];
#endif
  }
  // Moving the child pointers after (idx+1) in the current node one
  // step before
  for (uint32_t i=idx+2; i<=num_keys; ++i)
    children[i-1] = children[i];

  // Updating the key count of child and the current node
  child->num_keys += sibling->num_keys+1;
  num_keys--;

  // Freeing the memory occupied by sibling
  delete(sibling);
}

// A function to remove the key k from the sub-tree rooted with this node
template <class T, class W>
bool BTreeNode<T, W>::remove(T k) {
  uint32_t idx = findKey(k);

  // The key to be removed is present in this node
  if (idx < num_keys && keys[idx] == k) {
    // If the node is a leaf node - removeFromLeaf is called
    // Otherwise, removeFromNonLeaf function is called
    if (is_leaf)
      removeFromLeaf(idx);
    else
      removeFromNonLeaf(idx);
  } else {
    // If this node is a leaf node, then the key is not present in tree
    if (is_leaf) {
      return false;
    }
    // The key to be removed is present in the sub-tree rooted with this node
    // The flag indicates whether the key is present in the sub-tree rooted
    // with the last child of this node
    bool flag = ((idx == num_keys) ? true : false);

    // If the child where the key is supposed to exist has less that t keys,
    // we fill that child
    if (children[idx]->num_keys < MIN_KEYS)
      fill(idx);

    // If the last child has been merged, it must have merged with the previous
    // child and so we recurse on the (idx-1)th child. Else, we recurse on the
    // (idx)th child which now has atleast t keys
    if (flag && idx > num_keys)
      children[idx-1]->remove(k);
    else
      children[idx]->remove(k);
  }

  return true;
}

template <class T, class W>
bool BTree<T, W>::remove(T k) {
    if (!root)
        return false;

    // Call the remove function for root
    bool ret = root->remove(k);

    // If the root node has 0 keys, make its first child as the new root
    //  if it has a child, otherwise set root as NULL
    if (root->num_keys == 0) {
        BTreeNode<T, W> *tmp = root;
        if (root->is_leaf)
            root = NULL;
        else
            root = root->children[0];
        // Free the old root
        delete tmp;
    }
    return ret;
}

#endif
