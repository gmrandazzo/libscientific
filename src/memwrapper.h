#ifndef MEMWRAPPER_H
#define MEMWRAPPER_H

#include <stdio.h>
#include <stdlib.h>
#include "scientificinfo.h"

void *xmalloc(size_t size);
void *xrealloc(void *ptr, size_t size);
void xfree(void *ptr);

void GetNProcessor(size_t *nprocs_online, size_t *nprocs_max);

#endif
