#include "PMA_Lock.hpp"

//IMPORTANT!!!!!!!!! locks need to be initialized as all zeros
PMA_Lock::PMA_Lock() {
#if ENABLE_PMA_LOCK == 1
  x = 0;
  reason = GENERAL;
  reason_set_by = 0;
#ifndef NDEBUG
  owner = 0;
#endif
#endif
}

void PMA_Lock::name(char const about[10]){
#if ENABLE_PMA_LOCK == 1
#ifdef debug
  for (int i = 0; i < 10; i++) {
    info[i] = about[i];
  }
#endif
#endif
}


void PMA_Lock::number(int i){
#if ENABLE_PMA_LOCK == 1
#ifdef debug
  num = i;
#endif
#endif
}

// should only be used in sequntial areas
bool PMA_Lock::check_unlocked() {
#if ENABLE_PMA_LOCK == 1
  return x == 0;
#endif
  return true;
}


bool cas( uint64_t *old, uint64_t old_val, uint64_t new_val) {
  return __sync_bool_compare_and_swap(old, old_val, new_val);
} 


// exclusive lock
bool PMA_Lock::lock(uint64_t task_id, REASONS r) {
#if ENABLE_PMA_LOCK == 1
lock_start:
  assert(task_id > 0);
  assert(r != SAME);
  assert(num_tries > 0);
  int tries;
  bool success = false;


  while(!success) {
    tries = 0;
    while (tries < num_tries) {

      success = cas(&x, 0, 1);
      
      
      if(!success) {
        tries++;
      } else {
        break;
      }
    }
    // sleep in ms
    if(!success) {
      //usleep(wait_time);
      //sched_yield();
    }
  }

  if(r < reason) {
    cas(&x, 1, 0);
    dprintf("worker %lu is trying to lock %p, %s,%d, with the wrong reason %d, needed reason %d\n",task_id,this,info,num, r, reason);
    //usleep(wait_time);
    //sched_yield();
    dprintf("###############################################################################\n");
    goto lock_start;
    //return lock(task_id, r);
  }
#ifndef NDEBUG
  assert(x == 1);
  assert(owner == 0);
  owner = task_id;
#endif
#endif
  return true;
}

bool PMA_Lock::try_lock(uint64_t task_id, REASONS r) {
#if ENABLE_PMA_LOCK == 1
  assert(task_id > 0);
  assert(r != SAME);
  assert(num_tries > 0);
  int tries;
  bool success = false;

  tries = 0;
  while (tries < num_tries) {

    if (x == 0) {
      success = cas(&x, 0, 1);
      if(success) { break; }
    }

    tries++;
  }

  if(success && r < reason) {    
    cas(&x, 1, 0);
    dprintf("worker %lu is trying to lock %p, %s,%d, with the wrong reason %d, needed reason %d\n",task_id,this,info,num, r, reason);
    return false;
  }
#ifndef NDEBUG
  if (success) {
    assert(owner == 0);
    owner = task_id;
  }
#endif
  return success;
#endif
  return true;
}


void PMA_Lock::unlock(uint64_t task_id, REASONS r) {
#if ENABLE_PMA_LOCK == 1
  assert(task_id > 0);
#ifdef debug
  assert(num != -1);;
#endif
  assert(num_tries > 0);
  assert(owner == task_id);
  assert(x == 1);

  if (r != SAME && (reason_set_by == task_id || reason == GENERAL || reason_set_by == 0)) {
    reason = r;
    if (r != GENERAL) {
      reason_set_by = task_id;
    } else {
      reason_set_by = 0;
    }
  }
  


#ifdef debug
  assert(num != -1);;
#endif

#ifndef NDEBUG
  owner = 0;
  bool success = cas(&x, 1, 0);
  assert(success);
  
#else
  x = 0;
#endif
#endif
}

bool PMA_Lock::i_own_lock(uint64_t task_id) {
#if ENABLE_PMA_LOCK == 1
#ifdef debug
  assert(num != -1);;
#endif
#ifndef NDEBUG
  return owner == task_id;
#else
  printf("i_own_lock should never be called unless we are in debug mode");
  //exit(0);
  return true;
#endif
#endif
  return true;
}

void PMA_Lock::make_shared(uint64_t task_id) {
#if ENABLE_PMA_LOCK == 1
  assert(task_id > 0);
  assert(task_id == owner);
#ifdef debug
  assert(num != -1);;
#endif
#ifndef NDEBUG
  owner = 0;
#endif
  bool done = cas(&x, 1, 2);
  assert(done);
#endif
}

void PMA_Lock::make_exclusive(uint64_t task_id) {
#if ENABLE_PMA_LOCK == 1
  assert(task_id > 0);
#ifdef debug
  assert(num != -1);;
#endif
  assert(num_tries > 0);
  assert(owner == 0);
  bool success = false;
  int tries = 0;
  while(!success) {
    uint64_t new_x = x;
    dprintf("worker %lu is trying to make exclusive %p,%s,%d with current value %lx\n",task_id,this,info,num, new_x);
    tries = 0;
    if (new_x == 1) {
      while (tries < num_tries) {
        success = cas(&x, 1, 2);
        if (!success) {
          tries++;
          new_x = x;
        } else {
          break;
        }
      }
    }
    //usleep(wait_time);
    //sched_yield();
  }
#ifndef NDEBUG
  owner = task_id;
#endif
#endif
}

bool PMA_Lock::lock_shared(uint64_t task_id) {
#if ENABLE_PMA_LOCK == 1
  int tries;
  bool success = false;
  assert(owner == 0);
  while(!success) {
    tries = 0;
    uint64_t old_x = x;
    dprintf("SHARED worker %lu is trying to grab lock %p,%s,%d with current value %lx\n",task_id,this,info,num, old_x);
    if (x % 2 == 0) {
      while (tries < num_tries) {
        success = cas(&x, old_x, old_x + 2);
        if(!success) {
          tries++;
          old_x = x;
        } else {
          break;
        }
      }
    }
    // sleep in ms
    //usleep(wait_time);
    //sched_yield();
  }
  #ifdef debug
    uint64_t lock_val = x;
  #endif
  dprintf("SHARED worker %lu got lock %p with current value %lx\n",task_id,this, lock_val);
#endif
  return true;
}

void PMA_Lock::unlock_shared(uint64_t task_id) {
#if ENABLE_PMA_LOCK == 1
#ifdef debug
  assert(num != -1);;
#endif
  assert(num_tries > 0);
  #ifdef debug
    uint64_t old_x = x;
  #endif
  dprintf("SHARED worker %lu is trying to unlock %p,%s,%d with current value %lx\n",task_id,this,info,num, old_x);
  assert(x > 1);
  assert(owner == 0);
  assert(x%2==0);
  __sync_fetch_and_add(&x,-2);
#endif
}
