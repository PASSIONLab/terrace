#include <atomic>
#include <unistd.h>
#include <stdio.h>
#include "util.h"


#define num_tries 3
#define wait_time 1

#define ENABLE_PMA_LOCK 1

//TODO simplify since they don't need to be recursive anymore
//TODO get rid of locks being deleted

typedef enum reasons {GENERAL, REBALANCE, DOUBLE, SAME} REASONS;
class PMA_Lock {
public:
  PMA_Lock();

  REASONS reason;
  uint32_t reason_set_by;
  bool i_own_lock(uint64_t task_id);
  bool lock(uint64_t task_id, REASONS r = GENERAL);
  bool try_lock(uint64_t task_id, REASONS r = GENERAL);
  void unlock(uint64_t task_id, REASONS r = GENERAL);
  bool lock_shared(uint64_t task_id = 0);
  void unlock_shared(uint64_t task_id = 0);
  void name(char const about[10]);
  void number(int i);
  void make_shared(uint64_t task_id);
  void make_exclusive(uint64_t task_id);
  bool check_unlocked();
//private:
#if ENABLE_PMA_LOCK == 1
   uint64_t x;
#ifdef debug
  char info[10];
   int num = 0;
#endif
#ifndef NDEBUG
  uint64_t owner;
#endif
  #endif
};

