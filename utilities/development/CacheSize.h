
#ifndef CACHESIZE_H_INCLUDED
#define CACHESIZE_H_INCLUDED

// Returns the cache line size (in bytes) of the processor, or 0 on failure

#include <stddef.h>
size_t CacheLineSize();
size_t CacheL0Size();
size_t CacheL1Size();
size_t CacheL2Size();
size_t CacheL3Size();

#endif
