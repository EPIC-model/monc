#ifndef FT_DATA_H_
#define FT_DATA_H_

#include <pthread.h>
#include <stdlib.h>

#include "ft_consts.h"

#if defined(__DARWIN) && !defined(CLOCKID_ALREADY_DECLARED) 
typedef int clockid_t;
#endif

/**
 * A convenient array type
 **/
typedef struct array_tag {
  void **data;
  int size;
  int after;
  pthread_mutex_t mutex;
} array_t;

/**
 * A convenient volatile array type
 **/
typedef struct varray_tag {
  volatile void **data;
  int size;
  int after;
  pthread_mutex_t mutex;
} varray_t;

/**
 * global data structures
 **/

extern int is_initialized;

/**
 * holds the thread IDs
 **/
extern array_t *threads;

/** 
 * holds thread attributes, the index does not
 * necesseraly conincide with the one of threads.
 *
 * This array is just to allow different threads spawn new
 * threads at the same time.
 **/
extern array_t *thread_attrs;

/**
 * holds thread keys for storing thread specific data
 **/
extern array_t *thread_keys;

/*
 * holds the mutex IDs
 **/
extern array_t *mutexes;

/**
 * holds thread mutex attributes, the index does not
 * necesseraly concide with the one of mutexes.
 *
 * This array is just to allow different threads handle
 * their mutexes
 **/
extern array_t *mutex_attrs;

/**
 * an array to hold pthread_once_t structures
 **/
extern array_t *once_ctrls;

/**
 * an array to hold thread condition variable (pthread_cond_t) 
 * structures
 **/
extern array_t *conds;

/**
 * holds thread condition variable attributes, the index does not
 * necesseraly concide with the one of conds.
 *
 * This array is just to allow different threads handle
 * their condition variables
 **/
extern array_t *cond_attrs ;

/**
 * an array to hold thread barrier variable (pthread_barrier_t) 
 * structures
 **/
extern array_t *barriers;

/**
 * holds thread barrier variable attributes, the index does not
 * necesseraly concide with the one of conds.
 *
 * This array is just to allow different threads handle
 * their condition variables
 **/
extern array_t *barrier_attrs;

/**
 * an array to hold spinlock variable (pthread_spinlock_t) 
 * structures
 **/
extern varray_t *spinlocks;

/**
 * an array to hold thread read-write lock variable (pthread_rwlock_t) 
 * structures
 **/
extern array_t *rwlocks;

/**
 * holds thread rwlock variable attributes, the index does not
 * necesseraly concide with the one of conds.
 *
 * This array is just to allow different threads handle
 * their condition variables
 **/
extern array_t *rwlock_attrs;

void thread_init_internal(int *);
void thread_destroy_internal(int *);

///////////////////////////////////////////////////////////////

void array_init(array_t **array,int size);

void varray_init(varray_t **array,int size);

void array_resize(array_t **array,int size);

void varray_resize(varray_t **array,int size);

void array_delete(array_t *array);

void varray_delete(varray_t *array);

// This only works for pointer arrays!!
int is_valid(array_t *arr, int id);

// varray version
int vis_valid(varray_t *arr, int id);


#endif //FT_DATA_H_
