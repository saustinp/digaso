
#include "CacheSize.h"

#if defined(__APPLE__)

#include <sys/sysctl.h>
size_t CacheLineSize() {
	size_t lineSize = 0;
	size_t sizeOfLineSize = sizeof(lineSize);
	sysctlbyname("hw.cachelinesize", &lineSize, &sizeOfLineSize, 0, 0);
	return lineSize;
}

size_t CacheL0Size() {
    size_t L0size = 0;
    return L0size;
}

size_t CacheL1Size() {
    size_t L1size = 0;
    return L1size;
}

size_t CacheL2Size() {
    size_t L2size = 0;
    return L2size;
}

size_t CacheL2Size() {
    size_t L2size = 0;
    return L2size;
}


#elif defined(_WIN32)

#include <stdlib.h>
#include <windows.h>
size_t CacheLineSize() {
	size_t lineSize = 0;
	DWORD bufferSize = 0;
	DWORD i = 0;
	SYSTEM_LOGICAL_PROCESSOR_INFORMATION * buffer = 0;

	GetLogicalProcessorInformation(0, &bufferSize);
	buffer = (SYSTEM_LOGICAL_PROCESSOR_INFORMATION *) malloc(bufferSize);
	GetLogicalProcessorInformation(&buffer[0], &bufferSize);

	for (i = 0; i != bufferSize / sizeof(SYSTEM_LOGICAL_PROCESSOR_INFORMATION); ++i) {
		if (buffer[i].Relationship == RelationCache && buffer[i].Cache.Level == 1) {
			lineSize = buffer[i].Cache.LineSize;
			break;
		}
	}

	free(buffer);
	return lineSize;
}

size_t CacheL0Size() {
    size_t L0size = 0;
    return L0size;
}

size_t CacheL1Size() {
    size_t L1size = 0;
    return L1size;
}

size_t CacheL2Size() {
    size_t L2size = 0;
    return L2size;
}

size_t CacheL2Size() {
    size_t L2size = 0;
    return L2size;
}


#elif defined(__linux__)

#include <stdio.h>
size_t CacheLineSize() {
	FILE * p = 0;
	p = fopen("/sys/devices/system/cpu/cpu0/cache/index0/coherency_line_size", "r");
	unsigned int lineSize = 0;
	if (p) {
		fscanf(p, "%d", &lineSize);
		fclose(p);
	}
	return lineSize;
}

size_t CacheL0Size() {
	FILE * p = 0;
	p = fopen("/sys/devices/system/cpu/cpu0/cache/index0/size", "r");
    char units [1];
	unsigned int L0size = 0;
	if (p) {
		fscanf(p, "%d%s", &L0size, units);
		fclose(p);
	}
    if (units[0] == 'K')
        L0size *= 1024;
    else if (units[0] == 'M')
        L0size *= 1024*1024;
    else if (units[0] == 'G')
        L0size *= 1024*1024*1024;
    else
        L0size = 0;
	return L0size;
}

size_t CacheL1Size() {
	FILE * p = 0;
	p = fopen("/sys/devices/system/cpu/cpu0/cache/index1/size", "r");
    char units [1];
	unsigned int L1size = 0;
	if (p) {
		fscanf(p, "%d%s", &L1size, units);
		fclose(p);
	}
    if (units[0] == 'K')
        L1size *= 1024;
    else if (units[0] == 'M')
        L1size *= 1024*1024;
    else if (units[0] == 'G')
        L1size *= 1024*1024*1024;
    else
        L1size = 0;
	return L1size;
}

size_t CacheL2Size() {
	FILE * p = 0;
	p = fopen("/sys/devices/system/cpu/cpu0/cache/index2/size", "r");
    char units [1];
	unsigned int L2size = 0;
	if (p) {
		fscanf(p, "%d%s", &L2size, units);
		fclose(p);
	}
    if (units[0] == 'K')
        L2size *= 1024;
    else if (units[0] == 'M')
        L2size *= 1024*1024;
    else if (units[0] == 'G')
        L2size *= 1024*1024*1024;
    else
        L2size = 0;
	return L2size;
}

size_t CacheL3Size() {
	FILE * p = 0;
	p = fopen("/sys/devices/system/cpu/cpu0/cache/index3/size", "r");
    char units [1];
	unsigned int L3size = 0;
	if (p) {
		fscanf(p, "%d%s", &L3size, units);
		fclose(p);
	}
    if (units[0] == 'K')
        L3size *= 1024;
    else if (units[0] == 'M')
        L3size *= 1024*1024;
    else if (units[0] == 'G')
        L3size *= 1024*1024*1024;
    else
        L3size = 0;
	return L3size;
}

#else
#error Unrecognized platform
// size_t CacheLineSize() {
// 	size_t lineSize = 0;
// 	return lineSize;
// }
// 
// size_t CacheL0Size() {
//     size_t L0size = 0;
//     return L0size;
// }
// 
// size_t CacheL1Size() {
//     size_t L1size = 0;
//     return L1size;
// }
// 
// size_t CacheL2Size() {
//     size_t L2size = 0;
//     return L2size;
// }
// 
// size_t CacheL2Size() {
//     size_t L2size = 0;
//     return L2size;
// }

#endif
