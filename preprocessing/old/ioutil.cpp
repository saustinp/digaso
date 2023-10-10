#ifndef __IOUTIL
#define __IOUTIL

template <typename T> string NumberToString ( T Number )
{
    ostringstream ss;
    ss << Number;
    return ss.str();
}

void printiarray(vector<Int> a)
{    
    Int m = (Int) a.size();
    for (Int i=0; i<m; i++)
        cout << a[i] << "   ";
    cout << endl;
}

void print1iarray(Int* a, Int m)
{    
    for (Int i=0; i<m; i++)
        cout << a[i] << "   ";
    cout << endl;
}

void print2iarray(Int* a, Int m, Int n)
{
    for (Int i=0; i<m; i++) {
        for (Int j=0; j<n; j++)
            cout << a[j*m+i] << "   ";
        cout << endl;
    }
    cout << endl;
}

void print3iarray(Int* a, Int m, Int n, Int p)
{
    for (Int k=0; k<p; k++) {
        for (Int i=0; i<m; i++) {
            for (Int j=0; j<n; j++)
                cout << a[k*n*m+j*m+i] << "   ";
            cout << endl;
        }
        cout << endl;
    }
    cout << endl;
}

void print1darray(double* a, Int m)
{
    cout.precision(4);
    for (Int i=0; i<m; i++)
        cout << scientific << a[i] << "   ";
    cout << endl;
}

void print2darray(double* a, Int m, Int n)
{
    cout.precision(4);
    for (Int i=0; i<m; i++) {
        for (Int j=0; j<n; j++)
            cout << scientific << a[j*m+i] << "   ";
        cout << endl;
    }
    cout << endl;
}

void print3darray(double* a, Int m, Int n, Int p)
{
    cout.precision(4);
    for (Int k=0; k<p; k++) {
        for (Int i=0; i<m; i++) {
            for (Int j=0; j<n; j++)
                cout << scientific << a[k*n*m+j*m+i] << "   ";
            cout << endl;
        }
        cout << endl;
    }
    cout << endl;
}

void writeiarray(ofstream &out, vector<Int> &a, Int N)
{
    out.write( reinterpret_cast<char*>( &a[0] ), sizeof(Int) * N );
}

void writedarray(ofstream &out, vector<double> &a, Int N)
{
    out.write( reinterpret_cast<char*>( &a[0] ), sizeof(double) * N );
}

void writeiarray(ofstream &out, vector<Int> &a)
{
    Int N = (Int) a.size();
    if (N > 0)
        out.write( reinterpret_cast<char*>( &a[0] ), sizeof(Int) * N );
}

void writedarray(ofstream &out, vector<double> &a)
{
    Int N = (Int) a.size();
    if (N > 0)
        out.write( reinterpret_cast<char*>( &a[0] ), sizeof(double) * N );
}

void writeiarraytodouble(ofstream &out, vector<Int> &a, Int N)
{
    double b;
    for (unsigned i = 0; i < N; i++) {
        b = (double) a[i];
        out.write( reinterpret_cast<char*>( &b ), sizeof(double) );
    }
}

void writedarray2file(string filename, double *a, Int N)
{
    // Open file to read
    ofstream out(filename.c_str(), ios::out | ios::binary);

    if (!out) {
        error("Unable to open file " + filename);
    }

    out.write( reinterpret_cast<char*>( &a[0] ), sizeof(double) * N );

    out.close();
}

void writeiarray2file(string filename, Int *a, Int N)
{
    // Open file to read
    ofstream out(filename.c_str(), ios::out | ios::binary);
            
    if (!out) {
        error("Unable to open file " + filename);
    }

    out.write( reinterpret_cast<char*>( &a[0] ), sizeof(Int) * N );

    out.close();
}

void readcarray(ifstream &in, vector<char> &a, Int N)
{
    a.resize(N);

    in.read( reinterpret_cast<char*>( &a[0] ), sizeof(char)*N );
}

void readdarray(ifstream &in, vector<double> &a, Int N)
{
    a.resize(N);

    in.read( reinterpret_cast<char*>( &a[0] ), sizeof(double)*N );
}

void readiarray(ifstream &in, vector<Int> &a, Int N)
{
    a.resize(N);
    in.read( reinterpret_cast<char*>( &a[0] ), sizeof(Int) * N );
}

void readiarrayfromdouble(ifstream &in, vector<Int> &a, Int N)
{
    a.resize(N);

    double read;
    for (unsigned i = 0; i < N; i++) {
        in.read( reinterpret_cast<char*>( &read ), sizeof read );
        a[i] = (Int) round(read);
    }
}


//readInput(filein, elcon, t2t, t2f, elemall, entall, param);
void readInput(string filename, vector<Int> &param, vector<Int> &elcon, vector<Int> &t2t, 
        vector<Int> &t2f, vector<Int> &elemall, vector<Int> &entall)
{
    // Open file to read
    ifstream in(filename.c_str(), ios::in | ios::binary);

    if (!in) {
        error("Unable to open file " + filename);
    }

    if (in) {
        Int N = 20;
        
        /* Read dimensions to ndims array */
        readiarrayfromdouble(in, param, N);
        
        /* Get dimensions from param array */
        Int ne = param[0];
        Int nf = param[1];
        Int ndof = param[2];
        Int nfe = param[3];
        Int nee = param[4];
        Int npe = param[5];
        Int npf = param[6];                
        Int preconditioner = param[11];
        Int overlappinglevel = param[12]; 
        Int hdg = param[13];
        Int nproc = param[14];

        /* Read mesh data to mesh structure */
        N = npf*nfe*ne;
        readiarrayfromdouble(in, elcon, N);
        N = ne*nee;
        readiarrayfromdouble(in, t2t, N);        
        N = ne*nfe;
        readiarrayfromdouble(in, t2f, N);        
        N = ne;
        readiarrayfromdouble(in, elemall, N);        
        N = ndof;
        readiarrayfromdouble(in, entall, N);                
    }

    // Close file:
    in.close();
}

void writeOutput(string filename, dmdstruct &dmd)
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