#ifndef __ERRORMSG
#define __ERRORMSG

static void error(const char* errstr)
{
#ifdef  HAVE_MPI            // TODO: Print only with master processor
    cout << "Error: ";
    cout << errstr;
    cout << endl;

//    free(sys.requests); free(sys.statuses);       // TODO: Free variables allocated dynamically
//    free(sys.requestsEnt2entWeight); free(sys.statusesEnt2entWeight);
//     MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();
#else
    cout << "Error: ";
    cout << errstr;
    cout << endl;
#endif

    exit(1);
}

static void error(string errstr)
{
  error(errstr.c_str());
}

static void splitArray1D(Int arrayLen, Int this_thread, Int noThreads, Int* chunkLen, Int* startPosition)
{
    *startPosition = ((arrayLen * this_thread) / noThreads);
    Int endPosition = ((arrayLen * (this_thread+1)) / noThreads);
    
    if (this_thread == (noThreads-1)) {   // This is reduntant unless we go over the range in which double represent integers exactly
        if (endPosition != arrayLen) {
            printf("Warning 4TG67N in splitArray1D.\n");
            endPosition = arrayLen;
        }
    }
    
    *chunkLen = endPosition - *startPosition;
}

#endif
