#ifndef __READBINARYFILES
#define __READBINARYFILES

void readappstruct(string filename, appstruct &app)
{
    // Open file to read
    ifstream in(filename.c_str(), ios::in | ios::binary);

    if (!in) {
        error("Unable to open file " + filename);
    }

    if (in) {
        Int N = 30;
        
        /* Read dimensions to ndims array */
        readiarrayfromdouble(in, app.ndims, N);
        
        /* Read data to app structure */        
        readiarrayfromdouble(in, app.meshsize, app.ndims[0]);
        readiarrayfromdouble(in, app.solsize, app.ndims[1]);
        readiarrayfromdouble(in, app.porder, app.ndims[2]);
        readiarrayfromdouble(in, app.elemtype, app.ndims[3]);
        readiarrayfromdouble(in, app.nodetype, app.ndims[4]);
        readiarrayfromdouble(in, app.pgauss, app.ndims[5]);
        readiarrayfromdouble(in, app.pgaussR, app.ndims[6]);
        readiarrayfromdouble(in, app.quadtype, app.ndims[7]);
        readiarrayfromdouble(in, app.bcm, app.ndims[8]);
        readdarray(in, app.bcs, app.ndims[9]);
        readiarrayfromdouble(in, app.bcd, app.ndims[10]);
        readdarray(in, app.bcv, app.ndims[11]);
        readdarray(in, app.dt, app.ndims[12]);        
        readiarrayfromdouble(in, app.flag, app.ndims[13]);
        readdarray(in, app.factor, app.ndims[14]);       
        readiarrayfromdouble(in, app.problem, app.ndims[15]);
        readdarray(in, app.physicsparam, app.ndims[16]);       
        readdarray(in, app.solversparam, app.ndims[17]);        
        
        app.nproc = app.ndims[18];
        app.nfile = app.ndims[19];
    }

    // Close file:
    in.close();
}

void writeappstruct(string filename, appstruct &app)
{
    // Open file to write
    ofstream out(filename.c_str(), ios::out | ios::binary);

    if (!out) {
        error("Unable to open file " + filename);
    }

    if (out) {
        Int N = 30;
        
        /* Read dimensions to ndims array */
        writeiarraytodouble(out, app.ndims, N);
        
        /* Read data to app structure */        
        writeiarraytodouble(out, app.meshsize, app.ndims[0]);
        writeiarraytodouble(out, app.solsize, app.ndims[1]);
        writeiarraytodouble(out, app.porder, app.ndims[2]);
        writeiarraytodouble(out, app.elemtype, app.ndims[3]);
        writeiarraytodouble(out, app.nodetype, app.ndims[4]);
        writeiarraytodouble(out, app.pgauss, app.ndims[5]);
        writeiarraytodouble(out, app.pgaussR, app.ndims[6]);
        writeiarraytodouble(out, app.quadtype, app.ndims[7]);
        writeiarraytodouble(out, app.bcm, app.ndims[8]);
        writedarray(out, app.bcs, app.ndims[9]);
        writeiarraytodouble(out, app.bcd, app.ndims[10]);
        writedarray(out, app.bcv, app.ndims[11]);
        writedarray(out, app.dt, app.ndims[12]);        
        writeiarraytodouble(out, app.flag, app.ndims[13]);
        writedarray(out, app.factor, app.ndims[14]);       
        writeiarraytodouble(out, app.problem, app.ndims[15]);
        writedarray(out, app.physicsparam, app.ndims[16]);       
        writedarray(out, app.solversparam, app.ndims[17]);               
    }

    // Close file:
    out.close();
}

void readmasterstruct(string filename, vector< masterstruct > &master, appstruct &app)
{
    // Open file to read
    ifstream in(filename.c_str(), ios::in | ios::binary);

    if (!in) 
        error("Unable to open file " + filename);

    // begin reading data from the binary file
    Int i, j, elem, N = 20;        
    Int nelem = app.elemtype.size();
    Int dim, elemtype, nodetype, npe, nge, ngeR, nfe, nve, nle, nplmax;
    Int npf, ngf, ngfR, nvf;
            
    for (i=0; i<nelem; i++) {        
        elem = app.elemtype[i];        
        master[elem].porder = app.porder;
        master[elem].pgauss = app.pgauss;
        master[elem].pgaussR = app.pgaussR;
        
        /* Read dimensions to ndims array */
        readiarrayfromdouble(in, master[elem].ndims, N); 
        dim = master[elem].ndims[0];
        elemtype = master[elem].ndims[1];
        nodetype = master[elem].ndims[2];
        npe = master[elem].ndims[3];
        nge = master[elem].ndims[4];
        ngeR = master[elem].ndims[5];
        nfe = master[elem].ndims[6];
        nve = master[elem].ndims[7];
        nle = master[elem].ndims[8];
        nplmax = master[elem].ndims[9];
        
        master[elem].dim = master[elem].ndims[0];
        master[elem].elemtype = master[elem].ndims[1];
        master[elem].nodetype = master[elem].ndims[2];
        master[elem].npe = master[elem].ndims[3];
        master[elem].nge = master[elem].ndims[4];
        master[elem].ngeR = master[elem].ndims[5];
        master[elem].nfe = master[elem].ndims[6];
        master[elem].nve = master[elem].ndims[7];
        master[elem].nle = master[elem].ndims[8];
        master[elem].nplmax = master[elem].ndims[9];
        
        /* Read data to mesh structure */        
        readiarrayfromdouble(in, master[elem].npf, nfe);
        readiarrayfromdouble(in, master[elem].ngf, nfe);
        readiarrayfromdouble(in, master[elem].ngfR, nfe);
        readiarrayfromdouble(in, master[elem].nvf, nfe);
        readiarrayfromdouble(in, master[elem].npl, nle);        

        readdarray(in, master[elem].plocvl, npe*dim);
        readdarray(in, master[elem].gpvl, nge*dim);
        readdarray(in, master[elem].gwvl, nge);
        readdarray(in, master[elem].shapvt, npe*nge*(dim+1));
        readdarray(in, master[elem].shapvg, npe*nge*(dim+1));
        readdarray(in, master[elem].shapvgdotshapvl, npe*npe*nge*(dim+1));
        readdarray(in, master[elem].gpvlR, ngeR*dim);
        readdarray(in, master[elem].gwvlR, ngeR);
        readdarray(in, master[elem].shapvtR, npe*ngeR*(dim+1));
        readdarray(in, master[elem].shapvgR, npe*ngeR*(dim+1));        
        readiarrayfromdouble(in, master[elem].permnode, nve);
        readiarrayfromdouble(in, master[elem].permedge, nle*nplmax);                                    
        
        master[elem].plocfc.resize(nfe);
        master[elem].gpfc.resize(nfe);
        master[elem].gwfc.resize(nfe);
        master[elem].shapft.resize(nfe);
        master[elem].shapfg.resize(nfe);
        master[elem].shapfgdotshapfc.resize(nfe);
        master[elem].gpfcR.resize(nfe);
        master[elem].gwfcR.resize(nfe);
        master[elem].shapftR.resize(nfe);
        master[elem].shapfgR.resize(nfe);
        master[elem].perm.resize(nfe);
        master[elem].face.resize(nfe);
        //readdarray(in, master[elem].philocvl, npe*nve);
        for (j=0; j<nfe; j++) {
            npf = master[elem].npf[j];
            ngf = master[elem].ngf[j];            
            ngfR = master[elem].ngfR[j];
            nvf = master[elem].nvf[j];  
            readdarray(in, master[elem].plocfc[j], npf*(dim-1));
            readdarray(in, master[elem].gpfc[j], ngf*(dim-1));
            readdarray(in, master[elem].gwfc[j], ngf);
            readdarray(in, master[elem].shapft[j], npf*ngf*dim);
            readdarray(in, master[elem].shapfg[j], npf*ngf*dim);
            readdarray(in, master[elem].shapfgdotshapfc[j], npf*npf*ngf*dim);
            readdarray(in, master[elem].gpfcR[j], ngfR*(dim-1));
            readdarray(in, master[elem].gwfcR[j], ngfR);
            readdarray(in, master[elem].shapftR[j], npf*ngfR*dim);
            readdarray(in, master[elem].shapfgR[j], npf*ngfR*dim);                    
            readiarrayfromdouble(in, master[elem].perm[j], npf);
            readiarrayfromdouble(in, master[elem].face[j], nvf);    
            //readiarrayfromdouble(in, master[elem].philocfc[j], npf*nvf);            
        }                
    }    
    
    // Close file:
    in.close();
}


void writemasterstruct(string filename, vector< masterstruct > &master, appstruct &app)
{
    // Open file to write
    ofstream out(filename.c_str(), ios::out | ios::binary);

    if (!out) {
        error("Unable to open file " + filename);
    }

    // begin writing data to the binary file
    Int i, j, elem, N = 20;        
    Int nelem = app.elemtype.size();
    Int dim, elemtype, nodetype, npe, nge, ngeR, nfe, nve, nle, nplmax;
    Int npf, ngf, ngfR, nvf;
            
    for (i=0; i<nelem; i++) {        
        elem = app.elemtype[i];        
        
        /* write dimensions to ndims array */
        writeiarraytodouble(out, master[elem].ndims, N); 
        dim = master[elem].ndims[0];
        elemtype = master[elem].ndims[1];
        nodetype = master[elem].ndims[2];
        npe = master[elem].ndims[3];
        nge = master[elem].ndims[4];
        ngeR = master[elem].ndims[5];
        nfe = master[elem].ndims[6];
        nve = master[elem].ndims[7];
        nle = master[elem].ndims[8];
        nplmax = master[elem].ndims[9];
                
        /* write data to mesh structure */        
        writeiarraytodouble(out, master[elem].npf, nfe);
        writeiarraytodouble(out, master[elem].ngf, nfe);
        writeiarraytodouble(out, master[elem].ngfR, nfe);
        writeiarraytodouble(out, master[elem].nvf, nfe);
        writeiarraytodouble(out, master[elem].npl, nle);        

        writedarray(out, master[elem].plocvl, npe*dim);
        writedarray(out, master[elem].gpvl, nge*dim);
        writedarray(out, master[elem].gwvl, nge);
        writedarray(out, master[elem].shapvt, npe*nge*(dim+1));
        writedarray(out, master[elem].shapvg, npe*nge*(dim+1));
        writedarray(out, master[elem].shapvgdotshapvl, npe*npe*nge*(dim+1));
        writedarray(out, master[elem].gpvlR, ngeR*dim);
        writedarray(out, master[elem].gwvlR, ngeR);
        writedarray(out, master[elem].shapvtR, npe*ngeR*(dim+1));
        writedarray(out, master[elem].shapvgR, npe*ngeR*(dim+1));        
        writeiarraytodouble(out, master[elem].permnode, nve);
        writeiarraytodouble(out, master[elem].permedge, nle*nplmax);                                    
        
        //writedarray(out, master[elem].philocvl, npe*nve);
        for (j=0; j<nfe; j++) {
            npf = master[elem].npf[j];
            ngf = master[elem].ngf[j];            
            ngfR = master[elem].ngfR[j];
            nvf = master[elem].nvf[j];  
            writedarray(out, master[elem].plocfc[j], npf*(dim-1));
            writedarray(out, master[elem].gpfc[j], ngf*(dim-1));
            writedarray(out, master[elem].gwfc[j], ngf);
            writedarray(out, master[elem].shapft[j], npf*ngf*dim);
            writedarray(out, master[elem].shapfg[j], npf*ngf*dim);
            writedarray(out, master[elem].shapfgdotshapfc[j], npf*npf*ngf*dim);
            writedarray(out, master[elem].gpfcR[j], ngfR*(dim-1));
            writedarray(out, master[elem].gwfcR[j], ngfR);
            writedarray(out, master[elem].shapftR[j], npf*ngfR*dim);
            writedarray(out, master[elem].shapfgR[j], npf*ngfR*dim);                    
            writeiarraytodouble(out, master[elem].perm[j], npf);
            writeiarraytodouble(out, master[elem].face[j], nvf);    
            //writeiarraytodouble(out, master[elem].philocfc[j], npf*nvf);            
        }                
    }    
    
    // Close file:
    out.close();
}

void readmeshstruct(string filename, meshstruct &mesh, appstruct &app)
{
    // Open file to read
    ifstream in(filename.c_str(), ios::in | ios::binary);

    if (!in) 
        error("Unable to open file " + filename);
    
    if (app.nproc>1) { // parallel reading from the binary file       
        Int my_rank = app.my_rank;        // cpu rank
        Int nprocperfile = app.nproc/app.nfile;    
        Int filenumber = my_rank/nprocperfile+1; // file number        
        Int file_rank = my_rank%nprocperfile;    // file rank
        Int pos = (nprocperfile+1)*(filenumber-1);        
        Int *meshsize = &app.meshsize[pos];
                
        // jump to the right position in the binary file
        in.seekg(meshsize[file_rank]*sizeof(double), ios::beg);
        
        // number of entries to be read
        Int datasize = meshsize[file_rank+1]-meshsize[file_rank];
        
        // begin reading data from the binary file
        Int N = 30;        
        /* Read dimensions to ndims array */
        readiarrayfromdouble(in, mesh.ndims, N); 
        Int dim = mesh.ndims[0];
        Int ncd = mesh.ndims[1];
        Int nfe = mesh.ndims[2];
        Int nve = mesh.ndims[3];
        Int nvf = mesh.ndims[4];
        Int ne = mesh.ndims[5];
        Int nf = mesh.ndims[6];
        Int nv = mesh.ndims[7];
        Int ndh = mesh.ndims[8];
        Int npe = mesh.ndims[9];
        Int npf = mesh.ndims[10];
        Int nc = mesh.ndims[11];
        Int ncu = mesh.ndims[12];
        Int ncq = mesh.ndims[13];
        Int ncp = mesh.ndims[14];
        Int nch = mesh.ndims[15];
        Int nn = mesh.ndims[16];   // number of entities
        Int nie = mesh.ndims[17];  // number of interior elements
        Int nif = mesh.ndims[18];  // number of interior faces                    
        Int nin = mesh.ndims[19];  // number of interior entities                     
        Int nnbsd = mesh.ndims[20];
        Int nxe = mesh.ndims[21];  // number of exterior elements
        Int nxf = mesh.ndims[22];  // number of exterior faces     
        Int nxn = mesh.ndims[23];  // number of exterior entities            
        Int hybrid = mesh.ndims[N-2];
        Int nproc = mesh.ndims[N-1];        
                
        if (app.nproc != nproc) {
            printiarray(mesh.ndims);
            cout<<app.nproc<<"   "<<nproc<<endl;
            error("MPI world size must be equal to nproc");
        }
                        
        Int n = N+(1+nve+npe*ncd+npf*nfe+nfe+nfe+nfe+1+1)*ne+(nvf+3)*nf+3*nf+nnbsd;        
        if (hybrid != 0)
            n = n + 2*ndh;
        if (n != datasize) {
            cout<<my_rank<<"   "<<n<<"  "<<datasize<<endl;            
            printiarray(app.meshsize);
            printiarray(mesh.ndims);            
            error("Something wrong in the binary mesh file when reading.");
        }        
        
        /* Read data to mesh structure */        
        readiarrayfromdouble(in, mesh.elementtype, ne);
        readiarrayfromdouble(in, mesh.t, nve*ne);
        readdarray(in, mesh.dgnodes, npe*ncd*ne);
        readiarrayfromdouble(in, mesh.elcon, npf*nfe*ne);
        readiarrayfromdouble(in, mesh.bf, nfe*ne);
        readiarrayfromdouble(in, mesh.t2f, nfe*ne);
        readiarrayfromdouble(in, mesh.t2t, nfe*ne);
        readiarrayfromdouble(in, mesh.f, (nvf+3)*nf);         
        readiarrayfromdouble(in, mesh.isEDGface, nf);
        readiarrayfromdouble(in, mesh.nbsd, nnbsd);   
        readiarrayfromdouble(in, mesh.extintelem, ne);                           
        readiarrayfromdouble(in, mesh.elem2cpu, ne);
        readiarrayfromdouble(in, mesh.extintface, nf);
        readiarrayfromdouble(in, mesh.face2cpu, nf);         
        if (hybrid != 0) {
            readiarrayfromdouble(in, mesh.extintent, ndh);
            readiarrayfromdouble(in, mesh.ent2cpu, ndh);         
        }
    }
    else {
        // begin reading data from the binary file
        Int N = 30;        
        /* Read dimensions to ndims array */
        readiarrayfromdouble(in, mesh.ndims, N); 
        Int dim = mesh.ndims[0];
        Int ncd = mesh.ndims[1];
        Int nfe = mesh.ndims[2];
        Int nve = mesh.ndims[3];
        Int nvf = mesh.ndims[4];
        Int ne = mesh.ndims[5];
        Int nf = mesh.ndims[6];
        Int nv = mesh.ndims[7];
        Int ndh = mesh.ndims[8];
        Int npe = mesh.ndims[9];
        Int npf = mesh.ndims[10];
        Int nc = mesh.ndims[11];
        Int ncu = mesh.ndims[12];
        Int ncq = mesh.ndims[13];
        Int ncp = mesh.ndims[14];
        Int nch = mesh.ndims[15];
        Int nie = mesh.ndims[17];
        Int nif = mesh.ndims[18];
        Int hybrid = mesh.ndims[N-2];
        Int nproc = mesh.ndims[N-1];
                        
        /* Read data to mesh structure */        
        readiarrayfromdouble(in, mesh.elementtype, ne);
        readiarrayfromdouble(in, mesh.t, nve*ne);
        readdarray(in, mesh.dgnodes, npe*ncd*ne);
        readiarrayfromdouble(in, mesh.elcon, npf*nfe*ne);
        readiarrayfromdouble(in, mesh.bf, nfe*ne);
        readiarrayfromdouble(in, mesh.t2f, nfe*ne);
        readiarrayfromdouble(in, mesh.t2t, nfe*ne);                
        readiarrayfromdouble(in, mesh.f, (nvf+3)*nf);
        readiarrayfromdouble(in, mesh.isEDGface, nf);
    }
    
    // Close file:
    in.close();
}

void writemeshstruct(string filename, meshstruct &mesh, appstruct &app)
{    
    if (app.nproc>1) { // parallel writing to the binary file   
        
        //Int my_rank = app.my_rank;
        //Int *meshsize = &app.meshsize[0];
        Int my_rank = app.my_rank;        // cpu rank
        Int nprocperfile = app.nproc/app.nfile;    
        Int filenumber = my_rank/nprocperfile+1; // file number        
        Int file_rank = my_rank%nprocperfile;    // file rank
        Int pos = (nprocperfile+1)*(filenumber-1);        
        Int *meshsize = &app.meshsize[pos];            
        Int N = 30;                                
        
        //cout<<app.my_rank<<"   "<<app.nproc<<"   "<<app.nfile<<"   "<<nprocperfile<<"   "<<filenumber<<"   "<<file_rank<<endl;    
        //printiarray(mesh.ndims);
                
        // begin writeing data from the binary file        
        Int dim = mesh.ndims[0];
        Int ncd = mesh.ndims[1];
        Int nfe = mesh.ndims[2];
        Int nve = mesh.ndims[3];
        Int nvf = mesh.ndims[4];
        Int ne = mesh.ndims[5];
        Int nf = mesh.ndims[6];
        Int nv = mesh.ndims[7];
        Int ndh = mesh.ndims[8];
        Int npe = mesh.ndims[9];
        Int npf = mesh.ndims[10];
        Int nc = mesh.ndims[11];
        Int ncu = mesh.ndims[12];
        Int ncq = mesh.ndims[13];
        Int ncp = mesh.ndims[14];
        Int nch = mesh.ndims[15];
        Int nie = mesh.ndims[17];
        Int nif = mesh.ndims[18];
        Int nin = mesh.ndims[19];
        Int nnbsd = mesh.ndims[20];        
        Int hybrid = mesh.ndims[N-2];
        Int nproc = mesh.ndims[N-1];                        
        Int n = N+(1+nve+npe*ncd+npf*nfe+nfe+nfe+nfe+1+1)*ne+(nvf+3)*nf+3*nf+nnbsd;      
        Int datasize = meshsize[file_rank+1]-meshsize[file_rank];
                
        if (hybrid != 0)
            n = n + 2*ndh;
        if (n != datasize) {
            cout<<my_rank<<"   "<<n<<"  "<<datasize<<endl;
            printiarray(app.meshsize);
            printiarray(mesh.ndims);            
            error("Something wrong in the binary mesh file when writing");
        }        
                
//         char *fileout;
//         fileout = new char[filename.size() + 1];
//         memcpy(fileout, filename.c_str(), filename.size() + 1);   
//         MPI_File file;
//         MPI_Status status;          
//         MPI_File_open(MPI_COMM_WORLD, fileout, MPI_MODE_CREATE|MPI_MODE_WRONLY,
//                           MPI_INFO_NULL, &file);                
//         //MPI_Offset offset = sizeof(double)*meshsize[file_rank];
//         MPI_Offset offset = sizeof(double)*file_rank*N;
//         MPI_File_seek(file, offset, MPI_SEEK_SET);
//         vector<double> temp(N,0);
//         for (Int i=0; i<N; i++)
//             temp[i] = mesh.ndims[i];
//         //print1darray(&temp[0],N);
//         MPI_File_write(file, &temp[0], N, MPI_DOUBLE, &status);
//         MPI_File_close(&file);                                           
//         delete [] fileout;

        MPI_Barrier(MPI_COMM_WORLD);
        
        // Open file to write
        ofstream out(filename.c_str(), ios::out | ios::binary);

        if (!out) {
            error("Unable to open file " + filename);
        }
        
        // jump to the right position in the binary file
        out.seekp(meshsize[file_rank]*sizeof(double), ios::beg);
        
        MPI_Barrier(MPI_COMM_WORLD);
        
        /* write data to mesh structure */      
        writeiarraytodouble(out, mesh.ndims, N);         
        writeiarraytodouble(out, mesh.elementtype, ne);
        writeiarraytodouble(out, mesh.t, nve*ne);
        writedarray(out, mesh.dgnodes, npe*ncd*ne);
        writeiarraytodouble(out, mesh.elcon, npf*nfe*ne);
        writeiarraytodouble(out, mesh.bf, nfe*ne);
        writeiarraytodouble(out, mesh.t2f, nfe*ne);
        writeiarraytodouble(out, mesh.t2t, nfe*ne);
        writeiarraytodouble(out, mesh.f, (nvf+3)*nf);
        writeiarraytodouble(out, mesh.isEDGface, nf);
        writeiarraytodouble(out, mesh.nbsd, nnbsd);   
        writeiarraytodouble(out, mesh.extintelem, ne);      
        writeiarraytodouble(out, mesh.elem2cpu, ne);
        writeiarraytodouble(out, mesh.extintface, nf);              
        writeiarraytodouble(out, mesh.face2cpu, nf);                             
        if (hybrid != 0) {
            writeiarraytodouble(out, mesh.extintent, ndh);
            writeiarraytodouble(out, mesh.ent2cpu, ndh);         
        }                
        out.close();                
    }
    else {
        // Open file to write
        ofstream out(filename.c_str(), ios::out | ios::binary);

        if (!out) {
            error("Unable to open file " + filename);
        }
        
        // begin writeing data from the binary file
        Int N = 30;                
        Int dim = mesh.ndims[0];
        Int ncd = mesh.ndims[1];
        Int nfe = mesh.ndims[2];
        Int nve = mesh.ndims[3];
        Int nvf = mesh.ndims[4];
        Int ne = mesh.ndims[5];
        Int nf = mesh.ndims[6];
        Int nv = mesh.ndims[7];
        Int ndh = mesh.ndims[8];
        Int npe = mesh.ndims[9];
        Int npf = mesh.ndims[10];
        Int nc = mesh.ndims[11];
        Int ncu = mesh.ndims[12];
        Int ncq = mesh.ndims[13];
        Int ncp = mesh.ndims[14];
        Int nch = mesh.ndims[15];
        Int nie = mesh.ndims[17];
        Int nif = mesh.ndims[18];
        Int nproc = mesh.ndims[N-1];
                        
        /* write data to mesh structure */        
        writeiarraytodouble(out, mesh.ndims, N); 
        writeiarraytodouble(out, mesh.elementtype, ne);
        writeiarraytodouble(out, mesh.t, nve*ne);
        writedarray(out, mesh.dgnodes, npe*ncd*ne);
        writeiarraytodouble(out, mesh.elcon, npf*nfe*ne);
        writeiarraytodouble(out, mesh.bf, nfe*ne);
        writeiarraytodouble(out, mesh.t2f, nfe*ne);
        writeiarraytodouble(out, mesh.t2t, nfe*ne);                
        writeiarraytodouble(out, mesh.f, (nvf+3)*nf);
        writeiarraytodouble(out, mesh.isEDGface, nf);
        
        // Close file:
        out.close();        
    }        
}

void readsolstruct(string filename, solstruct &sol, appstruct &app)
{
    // Open file to read
    ifstream in(filename.c_str(), ios::in | ios::binary);

    if (!in) 
        error("Unable to open file " + filename);
    
    if (app.nproc>1) { // parallel reading from the binary file       
        Int my_rank = app.my_rank;        // cpu rank
        Int nprocperfile = app.nproc/app.nfile;    
        Int filenumber = my_rank/nprocperfile+1; // file number        
        Int file_rank = my_rank%nprocperfile;    // file rank
        Int pos = (nprocperfile+1)*(filenumber-1);        
        Int *solsize = &app.solsize[pos];
                
        // jump to the right position in the binary file
        in.seekg(solsize[file_rank]*sizeof(double), ios::beg);
        
        // number of entries to be read
        Int datasize = solsize[file_rank+1]-solsize[file_rank];
        
        // begin reading data from the binary file
        Int N = 10;        
        /* Read dimensions to ndims array */
        readiarrayfromdouble(in, sol.ndims, N); 
        Int nudg = sol.ndims[0];
        Int nuh = sol.ndims[1];
        Int npdg = sol.ndims[2];        
        Int n = N+nudg+nuh+npdg;        
        
        if (n != datasize) {
            cout<<my_rank<<"   "<<n<<"  "<<datasize<<endl;
            printiarray(app.solsize);
            printiarray(sol.ndims);            
            error("Something wrong in the binary file");
        }        
        
        /* Read data to mesh structure */                
        readdarray(in, sol.UDG, nudg);
        readdarray(in, sol.UH, nuh);
        readdarray(in, sol.PDG, npdg);
    }
    else {
        // begin reading data from the binary file
        Int N = 10;        
        /* Read dimensions to ndims array */
        readiarrayfromdouble(in, sol.ndims, N); 
        Int nudg = sol.ndims[0];
        Int nuh = sol.ndims[1];
        Int npdg = sol.ndims[2];        
                        
        /* Read data to mesh structure */        
        readdarray(in, sol.UDG, nudg);
        readdarray(in, sol.UH, nuh);
        readdarray(in, sol.PDG, npdg);
    }
    
    // Close file:
    in.close();
}


void writesolstruct(string filename, solstruct &sol, appstruct &app)
{   
    if (app.nproc>1) { // parallel writing from the binary file       
        Int my_rank = app.my_rank;        // cpu rank
        Int nprocperfile = app.nproc/app.nfile;    
        Int filenumber = my_rank/nprocperfile+1; // file number        
        Int file_rank = my_rank%nprocperfile;    // file rank
        Int pos = (nprocperfile+1)*(filenumber-1);        
        Int *solsize = &app.solsize[pos];
                        
        // begin writeing data from the binary file
        Int N = 10;        
        Int nudg = sol.ndims[0];
        Int nuh = sol.ndims[1];
        Int npdg = sol.ndims[2];        
        Int n = N+nudg+nuh+npdg;                        
        Int datasize = solsize[file_rank+1]-solsize[file_rank];
        
        if (n != datasize) {
            cout<<my_rank<<"   "<<n<<"  "<<datasize<<endl;
            printiarray(app.solsize);
            printiarray(sol.ndims);            
            error("Something wrong in the binary file");
        }        

        MPI_Barrier(MPI_COMM_WORLD);
        
        // Open file to write
        ofstream out(filename.c_str(), ios::out | ios::binary);

        if (!out) {
            error("Unable to open file " + filename);
        }
                
	// jump to the right position in the binary file
        out.seekp(solsize[file_rank]*sizeof(double), ios::beg);
        
        MPI_Barrier(MPI_COMM_WORLD);
                
        /* write data to mesh structure */                
        writeiarraytodouble(out, sol.ndims, N);         
        writedarray(out, sol.UDG, nudg);
        writedarray(out, sol.UH, nuh);
        writedarray(out, sol.PDG, npdg);
        out.close();    
    }
    else {
        // Open file to write
        ofstream out(filename.c_str(), ios::out | ios::binary);

        if (!out) {
            error("Unable to open file " + filename);
        }

        // begin writeing data from the binary file
        Int N = 10;        
        Int nudg = sol.ndims[0];
        Int nuh = sol.ndims[1];
        Int npdg = sol.ndims[2];        
                                                
        /* write data to mesh structure */        
        writeiarraytodouble(out, sol.ndims, N);         
        writedarray(out, sol.UDG, nudg);
        writedarray(out, sol.UH, nuh);
        writedarray(out, sol.PDG, npdg);       
        out.close();    
    }    
}

void writedmdstruct(string filename, dmdstruct &dmd)
{
    // Open file to write
    ofstream out(filename.c_str(), ios::out | ios::binary);

    if (!out) {
        error("Unable to open file " + filename);
    }

    if (out) {                                
        Int N = 50; 
        vector<Int> ndims(N,0);                    
        ndims[0] = N;
        ndims[1] = dmd.my_rank;
        ndims[2] = (Int) dmd.nbsd.size();
        ndims[3] = (Int) dmd.intelem.size();
        ndims[4] = (Int) dmd.intent.size();
        ndims[5] = (Int) dmd.elempart.size();
        ndims[6] = (Int) dmd.elempartpts.size();
        ndims[7] = (Int) dmd.entpart.size();
        ndims[8] = (Int) dmd.entpartpts.size();
        ndims[9] = (Int) dmd.elemrecv.size();
        ndims[10] = (Int) dmd.elemrecvpts.size();
        ndims[11] = (Int) dmd.elemsend.size();
        ndims[12] = (Int) dmd.elemsendpts.size();        
        ndims[13] = (Int) dmd.entrecv.size();
        ndims[14] = (Int) dmd.entrecvpts.size();
        ndims[15] = (Int) dmd.entsend.size();
        ndims[16] = (Int) dmd.entsendpts.size();
        ndims[17] = (Int) dmd.vecrecv.size();
        ndims[18] = (Int) dmd.vecrecvpts.size();
        ndims[19] = (Int) dmd.vecsend.size();
        ndims[20] = (Int) dmd.vecsendpts.size();
        ndims[21] = (Int) dmd.matrecv.size();
        ndims[22] = (Int) dmd.matrecvpts.size();
        ndims[23] = (Int) dmd.matsend.size();
        ndims[24] = (Int) dmd.matsendpts.size();
        ndims[25] = (Int) dmd.rowent2elem.size();
        ndims[26] = (Int) dmd.colent2elem.size();
        ndims[27] = (Int) dmd.rowent2ent.size();
        ndims[28] = (Int) dmd.colent2ent.size();
        ndims[29] = (Int) dmd.bcrs_rowent2elem.size();
        ndims[30] = (Int) dmd.bcrs_colent2elem.size();
        ndims[31] = (Int) dmd.bcrs_rowent2ent.size();
        ndims[32] = (Int) dmd.bcrs_colent2ent.size();        
        ndims[33] = (Int) dmd.ent2ind.size();
        ndims[34] = (Int) dmd.elcon.size();
        ndims[35] = (Int) dmd.t2f.size();
        
        /* Write dmd structure to a file */    
        writeiarray(out, ndims);
        writeiarray(out, dmd.nbsd);
        writeiarray(out, dmd.intelem);
        writeiarray(out, dmd.intent);
        writeiarray(out, dmd.elempart);
        writeiarray(out, dmd.elempartpts);
        writeiarray(out, dmd.entpart);
        writeiarray(out, dmd.entpartpts);
        writeiarray(out, dmd.elemrecv);
        writeiarray(out, dmd.elemrecvpts);
        writeiarray(out, dmd.elemsend);
        writeiarray(out, dmd.elemsendpts);
        writeiarray(out, dmd.entrecv);
        writeiarray(out, dmd.entrecvpts);
        writeiarray(out, dmd.entsend);
        writeiarray(out, dmd.entsendpts);
        writeiarray(out, dmd.vecrecv);
        writeiarray(out, dmd.vecrecvpts);
        writeiarray(out, dmd.vecsend);
        writeiarray(out, dmd.vecsendpts);
        writeiarray(out, dmd.matrecv);
        writeiarray(out, dmd.matrecvpts);
        writeiarray(out, dmd.matsend);
        writeiarray(out, dmd.matsendpts);
        writeiarray(out, dmd.rowent2elem);
        writeiarray(out, dmd.colent2elem);
        writeiarray(out, dmd.rowent2ent);
        writeiarray(out, dmd.colent2ent);
        writeiarray(out, dmd.bcrs_rowent2elem);
        writeiarray(out, dmd.bcrs_colent2elem);
        writeiarray(out, dmd.bcrs_rowent2ent);
        writeiarray(out, dmd.bcrs_colent2ent);
        writeiarray(out, dmd.ent2ind);
        writeiarray(out, dmd.elcon);
        writeiarray(out, dmd.t2f);
    }

    // Close file:
    out.close();
}

#endif