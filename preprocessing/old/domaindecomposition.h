#ifndef __DOMAINDECOMPOSITION_H
#define __DOMAINDECOMPOSITION_H

#include "typedefint.h"
        
struct dmdstruct {    
    Int my_rank; // index of the processor
    vector<Int> nbsd;    // neighboring cpus
    vector<Int> intelem; // nonoverlapping global elements
    vector<Int> intent;  // nonoverlapping global entities
    vector<Int> elempart; // overlapping global elements
    vector<Int> elempartpts; // classifiers of overlapping global elements: (interior, interface, exterior)
    vector<Int> entpart;  // overlapping global entities
    vector<Int> entpartpts; // classifiers of overlapping global entities: (interior, interface, exterior)
    vector<Int> elemrecv;  // local elements received from neighboring cpus
    vector<Int> elemrecvpts; // classifiers of local elements received from neighboring cpus: (# elements from ncpu1, ...) 
    vector<Int> elemsend;  // local elements sent to neighboring cpus
    vector<Int> elemsendpts; // classifiers of local elements sent to neighboring cpus: (# elements to ncpu1, ...)    
    vector<Int> entrecv;    // local entities received from neighboring cpus
    vector<Int> entrecvpts; // classifiers of local entities received from neighboring cpus: (# entities from ncpu1, ...)          
    vector<Int> entsend;  // local entities sent to neighboring cpus
    vector<Int> entsendpts; // classifiers of local entities sent to neighboring cpus: (# elements to ncpu1, ...)    
    vector<Int> vecrecv;  // local vectors received from neighboring cpus to perform matrix-vector product
    vector<Int> vecrecvpts; // classifiers of local vectors received from neighboring cpus: (# vectors from ncpu1, ...)          
    vector<Int> vecsend;    // local vectors sent to neighboring cpus to perform matrix-vector product
    vector<Int> vecsendpts; // classifiers of local vectors sent to neighboring cpus: (# vectors to ncpu1, ...)             
    vector<Int> matrecv;    // local matrices received from neighboring cpus to construct the preconditioner
    vector<Int> matrecvpts; // classifiers of local matrices received from neighboring cpus: (# matrices from ncpu1, ...)                
    vector<Int> matsend;   // local matrices sent to neighboring cpus to construct the preconditioner
    vector<Int> matsendpts;  // classifiers of local matrices sent to neighboring cpus: (# matrices to ncpu1, ...)                                         
    vector<Int> rowent2elem; // global entity-to-element connectivities
    vector<Int> colent2elem; // global entity-to-element connectivities
    vector<Int> rowent2ent;  // global entity-to-entity connectivities
    vector<Int> colent2ent;  // global entity-to-entity connectivities
    vector<Int> bcrs_rowent2elem; // local entity-to-element connectivities
    vector<Int> bcrs_colent2elem; // local entity-to-element connectivities
    vector<Int> bcrs_rowent2ent;  // local entity-to-entity connectivities
    vector<Int> bcrs_colent2ent;  // local entity-to-entity connectivities    
    vector<Int> ent2ind;  // global-to-local entity mapping  
    vector<Int> elcon;    // local element-to-entity connectivities
    vector<Int> t2f;      // local elemeent-to-face connectivities
    vector<Int> t2t;      // local elemeent-to-element  connectivities
    Int  maxBlocksPerRow; // maximum number of entities per row
    Int  minBlocksPerRow; // minimum number of entities per row    
};

//int domaindecomposition(dmdstruct &dmd, const vector<Int> &elcon, const vector<Int> &t2t, const vector<Int> &t2f, Int nproc, Int preconditioner, const vector<Int> &elconhdg);

//typedef struct dmdstruct dmdstruct;
//int domaindecomposition(dmdstruct dmd, const vector<Int> elcon, const vector<Int> t2t, const vector<Int> t2f, Int nproc, Int preconditioner, const vector<Int> elconhdg)
// dmd = domaindecomposition(elcon,t2t,t2f,nproc,preconditioner,elconhdg)

#endif