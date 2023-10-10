/************************************************************
 * Script for benchmarking different matvec() implementations
 ***********************************************************/

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <mpi.h>

#ifndef HAVE_MPI
#define HAVE_MPI
#endif

#define TIME(stmts) ({long startTime = getMicrotime(); stmts; getMicrotime() - startTime;})

#include "wrapper.h"
// #include "MKLwrapper.h"
// #include "../blas/mkl/include/mkl.h"
#include "matvec1.cpp"
#include "matvec2.cpp"
#include "matvec3.cpp"
#include "matvec4.cpp"
#include "matvec4_prefetch.cpp"

void prematvec() {   
// #ifdef  HAVE_MPI
//     MPI_Barrier( MPI_COMM_WORLD );
// #endif
    
    usleep(100);  
}


/*************************************
 * Returns the current time in microseconds.
 ************************************/
long getMicrotime(){
	struct timeval currentTime;
	gettimeofday(&currentTime, NULL);
	return currentTime.tv_sec * (int)1e6 + currentTime.tv_usec;
}


/************************************************************************
 * Generate a sparse block matrix with random values but a predictable
 * structure, and a random vector. `matrixHeight` and `matrixEffectiveWidth` are the
 * height and width of the desired matrix (in terms of individual rows and
 * columns, not whole blocks).  The matrix will be stored in condensed form in
 * `*Hg`. Random values will also be generated and stored in `*x`, which will
 * be the vector that Hg should be multiplied into. If `generateFullMatrix` is
 * true, then a full version of the matrix (with all of the zeros) will be
 * stored in `*fullHg` (this is useful if you want to generate the canonically
 * correct result of `v = Hg * x` using a single call to `SGEMV()`, but is not
 * useful/practical for large enough matrices).
 ***********************************************************************/
sysstruct generateSystem(
	int matrixHeight, int matrixEffectiveWidth, float **fullHg, float **Hg, float **x){
	int blockDimension = 5;
	int blockSize2 = blockDimension * blockDimension;
	
	int numBlocksPerRow = matrixEffectiveWidth / blockDimension;
	int numBlockRows = matrixHeight / blockDimension;
	int numBlocks = numBlockRows * numBlocksPerRow;

	*Hg = (float *)calloc(numBlocks * blockSize2, sizeof(float));
	*x = (float *)calloc(matrixHeight, sizeof(float));
	Int *row_blk = (Int *)calloc(numBlockRows + 1, sizeof(Int));
	Int *col_ind = (Int *)calloc(numBlocks, sizeof(Int));
    
	int numElementsTotal = 0;
	int numBlocksTotal = 0;
	for(int blockRowNum = 0; blockRowNum < numBlockRows; blockRowNum++){
		row_blk[blockRowNum] = numBlocksTotal;
        
		for(int numBlockInRow = 0; numBlockInRow < numBlocksPerRow; numBlockInRow++){
			int blockIV = rand() % 10;
			for(int elementNum = 0; elementNum < blockSize2; elementNum++){
				(*Hg)[numElementsTotal++] = (float) (blockIV + elementNum);
			}
            
            col_ind[numBlocksTotal++] = rand() % numBlockRows;//numBlockInRow * blockSpacingInterval / blockDimension;
		}
	}
    
	for(int ind = 0; ind < matrixHeight; ind++){
		(*x)[ind] = (float) (rand() % 10);
	}

	row_blk[numBlockRows] = numBlocksTotal;
    
	sysstruct sys;
	sys.numRows = numBlockRows;
	sys.blkSize = 5;
	sys.row_blk = row_blk;
	sys.col_ind = col_ind;
	sys.xDenseRow = (float *)calloc(numBlocksPerRow*blockDimension, sizeof(float));
	return sys;
}


/************************************************************************
 * Test the first `length` elements of `actual` and `expected` for equality
 * (within a small float tolerance). If they're equal, return `true`; otherwise,
 * print out a message and return `false`.
 ***********************************************************************/
bool testVecEquality(float *expected, float *actual, int length){
	for(int ind = 0; ind < length; ind++){
		float valDiff = actual[ind] - expected[ind];
		if(valDiff < 0){
			valDiff = -valDiff;
		}
		if(valDiff > 1.0e-6){
			printf("Error at index %d: %f (actual) != %f (expected)\n", ind, actual[ind], expected[ind]);
			return false;
		}
	}

	return true;
}



/***********************************
 * Benchmarks matvec implementations
 ***********************************/
void test(){
    
    int world_size, world_rank;
    
#ifdef  HAVE_MPI
    // Get the number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
 
    // Get the rank of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    if (world_size <= 1) {
        printf("Number of processors must be more than 1.");
        exit(-1);
    }
#else
    world_size = 1;
    world_rank = 0;
#endif
    
	if (world_rank == 0) {
        printf("\nStarting benchmark...\n\n");
    }
    
	srand(0);
    int blockSize = 5;
	int matrixHeight = 5000 * blockSize;
    int matrixEffectiveWidth = 100 * blockSize;
	int numTrials = 500;
    
	float *Hg, *x;
	sysstruct sys = generateSystem( matrixHeight, matrixEffectiveWidth, NULL, &Hg, &x);
	float *v1 = (float *)calloc(matrixHeight, sizeof(float));
// 	float *v2 = (float *)calloc(matrixHeight, sizeof(float));
// 	float *v3 = (float *)calloc(matrixHeight, sizeof(float));
// 	float *v4 = (float *)calloc(matrixHeight, sizeof(float));
    
	double dataLoaded = ((double) matrixHeight / (double) blockSize) * ((double) matrixEffectiveWidth / (double) blockSize) * ((double) (blockSize * blockSize + blockSize)) * (double) sizeof(float);	// In Bytes
	double flops = 2.0 * (double) matrixHeight * (double) matrixEffectiveWidth;
    
// 	matvec1(v1, Hg, x, sys);
// 	matvec2(v2, Hg, x, sys);
// 	matvec3(v3, Hg, x, sys);
// 	matvec4(v4, Hg, x, sys);

	// Check whether the new implementation generated correct values
// 	if(testVecEquality(v2, v3, matrixHeight)){
// 		puts("Vector values match");
// 	}
// 	else {
// 		puts("ERROR: vector values mismatch");
// 		return;
// 	}
    
	long long totalTime[12] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    double avgTime[12] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    int numPreftch = 4;
    int preftchDist[4] = { 4*1 , 4*5 , 4*20 , 4*100 };
	for(int trial = 0; trial < numTrials; trial++){
		// TODO: Remove everything from Cache before calling each of the matvec routines
        if (world_rank == 0)
            fprintf(stderr, "\rTrial %d / %d", trial+1, numTrials);
        
        prematvec();
		totalTime[0] += TIME(matvec1(v1, Hg, x, sys));
        
        prematvec();
        totalTime[1] += TIME(matvec2(v1, Hg, x, sys, 0));
        
        prematvec();
		totalTime[2] += TIME(matvec3(v1, Hg, x, sys));
        
        prematvec();
		totalTime[3] += TIME(matvec4(v1, Hg, x, sys));
        
        for (int preftchIndex = 0; preftchIndex < numPreftch; preftchIndex++) { 
            prematvec();
            totalTime[4+preftchIndex] += TIME(matvec4_prefetch(v1, Hg, x, sys, preftchDist[preftchIndex]));
        }
	}
	puts("");
    
	avgTime[0] = (double) totalTime[0] / (double) numTrials;
    avgTime[1] = (double) totalTime[1] / (double) numTrials;
	avgTime[2] = (double) totalTime[2] / (double) numTrials;
    avgTime[3] = (double) totalTime[3] / (double) numTrials;
    for (int preftchIndex = 0; preftchIndex < numPreftch; preftchIndex++)
        avgTime[4+preftchIndex] =(double)  totalTime[4+preftchIndex] / (double) numTrials;
    
    if (world_rank == 0) {
        printf("\n***************************\n***************************\n");
        printf("CASE DESCRIPTION:\n");
        printf("Matrix size: %d x %d (%d x %d blocks).\n", matrixHeight, matrixHeight, matrixHeight / blockSize, matrixHeight / blockSize);
        printf("Number of non-zero blocks per row: %d.\n", matrixEffectiveWidth / blockSize);
        printf("Block size: %d.\n", blockSize);
        printf("Number of trials to generate average timings: %d.\n\n", numTrials);
        printf("AVERAGE TIMINGS ON PROCESSOR NO. %d:\n", world_rank);
        printf("Method #1 (BLAS per block, no prefetching): %g us\n", avgTime[0]);
        printf("Method #2 (BLAS per row, no prefetching): %g us\n", avgTime[1]);
        printf("Method #3 (in-house, long code, no prefetching): %g us\n", avgTime[2]);
        printf("Method #4 (in-house, short code, no prefetching): %g us\n", avgTime[3]);
        for (int preftchIndex = 0; preftchIndex < numPreftch; preftchIndex++)
            printf("Method #%d (in-house, short code, prefetching %d blocks in advance): %g us\n", 5+preftchIndex, preftchDist[preftchIndex], avgTime[4+preftchIndex]);
        printf("\n");
        printf("\n***************************\n***************************\n");
        printf("Data loaded to cache per matrix-vector product: %g MBytes.\n", dataLoaded / (1024.0*1024.0));
        printf("Performance with Method No. %d: %g GBytes/s | %g GFLOPS\n", 4, 
                (dataLoaded / (1024.0*1024.0*1024.0)) / (avgTime[3] / 1.0e6), 
                (flops / (1024.0*1024.0*1024.0)) / (avgTime[3] / 1.0e6));
        printf("\n");
    }
    
#ifdef  HAVE_MPI
    if (world_rank == 0)
        printf("\nAVERAGE PERFORMANCE OF METHOD %d ON EACH PROCESSOR:\n", 4);
    MPI_Barrier( MPI_COMM_WORLD );
    printf("Processor No. %d / %d: %g GBytes/s\n", world_rank, world_size,
            (dataLoaded / (1024.0*1024.0*1024.0)) / (avgTime[3] / 1.0e6));
#endif
    printf("\n");
}


int main() {
    
#ifdef  HAVE_MPI
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);
#endif
    
    // Run benchmark:
	test();
    
#ifdef  HAVE_MPI
    // Finalize the MPI environment:
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();
#endif
    
	return 0;
}
