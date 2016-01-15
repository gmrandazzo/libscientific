#include "memwrapper.h"
#include <stdio.h>
#include <stdlib.h>

#ifdef WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <unistd.h>
#endif


inline void *xmalloc(size_t size)
{
  void *ptr = NULL;
  ptr = malloc(size);
  if(ptr == NULL){
    fprintf(stderr, "Memory Exhausted!\n");
    abort();
  }
  return ptr;
}

void *xrealloc(void *ptr, size_t size)
{
  register void *value = realloc (ptr, size);
  if (value == 0){
    fprintf(stderr, "Memory Exhausted!\n");
    abort();
  }
  return value;
}

void xfree(void *ptr)
{
  free(ptr);
}

void GetNProcessor(size_t *nprocs_online, size_t *nprocs_max)
{
  if(nprocs_online != NULL)
    (*nprocs_online) = -1;
  
  if(nprocs_max != NULL)
    (*nprocs_max) = -1;
  #ifdef WIN32
  #ifndef _SC_NPROCESSORS_ONLN
  SYSTEM_INFO info;
  GetSystemInfo(&info);
  #define sysconf(a) info.dwNumberOfProcessors
  #define _SC_NPROCESSORS_ONLN
  #endif
  #endif
  
  #ifdef _SC_NPROCESSORS_ONLN
  if(nprocs_online != NULL){
    (*nprocs_online) = sysconf(_SC_NPROCESSORS_ONLN);
    if ((*nprocs_online) < 1){
      (*nprocs_online) = 1;
    }
  }
  
  if(nprocs_max != NULL){
    (*nprocs_max) = sysconf(_SC_NPROCESSORS_CONF);
    if((*nprocs_max) < 1){
      (*nprocs_max) = 1;
    }
  }
  #else
  if(nprocs_online != NULL)
    (*nprocs_online) = 1;
  
  if(nprocs_max != NULL)
    (*nprocs_max) = 1
  #endif
}
