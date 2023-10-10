#ifndef __IOUTILITIES
#define __IOUTILITIES

// Written by: C. Nguyen & P. Fernandez

template <typename T> string NumberToString ( T Number )
{
    ostringstream ss;
    ss << Number;
    return ss.str();
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

void solstruct::writeSol2File(meshstruct &mesh, string filename, Int ne, Int lenUH, Int* ndims, Int writeQflag)
{
    // Open file to write
    ofstream out(filename.c_str(), ios::out | ios::binary);

    if (!out)
        cout <<"Unable to open file" << filename << endl;

    if (out) {
        // Write udg
        if (writeQflag == 0) {
            Int UDG2WriteStart = 0, UDGstart = 0, inc = 1;
            for (int i=0; i<ne; i++) {
                Int elementype = mesh.elementype[ie];
                Int npv = mesh.elements[elementype].npv;
                Int nc  = ndims[19];
                Int ncu = ndims[20];
                Int len = npv*ncu;
                
                DCOPY(&len, &UDG[UDGstart], &inc, &UDG2Write[UDG2WriteStart], &inc);
                UDGstart += npv*nc;
                UDG2WriteStart += npv*ncu;
            }
            out.write( reinterpret_cast<char*>( &UDG2Write[0] ), sizeof(double) * UDG2WriteStart );
        }
        else if (writeQflag == 1) {
            Int lenUDG = 0;
            for (int i=0; i<ne; i++) {
                Int elementype = mesh.elementype[ie];
                Int npv = mesh.elements[elementype].npv;
                Int nc  = ndims[19];
                lenUDG += npv*nc;
            }
            out.write( reinterpret_cast<char*>( &UDG[0] ), sizeof(double) * lenUDG );
        }
        else {
            printf("writeQflag has invalid value in writeSol2File.\n");
            exit(-1);
        }
        
        // Write UH:
        out.write( reinterpret_cast<char*>( &UH[0] ), sizeof(double) * lenUH );
    }

    // Close file:
    out.close();
}

void solstruct::writeAvgSol2File(meshstruct &mesh, string filename_avg, Int ne, Int lenUH_avg, Int* ndims, Int writeQflag)
{
    // Open file to write:
    ofstream out(filename_avg.c_str(), ios::out | ios::binary);

    if (!out)
        cout <<"Unable to open file" << filename_avg << endl;

    // If this is a valid file:
    if (out) {
        // Write udg_avg:
        if (writeQflag == 0) {
            Int UDG2WriteStart = 0, UDG_avgStart = 0, inc = 1;
            for (int i=0; i<ne; i++) {
                Int elementype = mesh.elementype[ie];
                Int npv = mesh.elements[elementype].npv;
                Int nc  = ndims[19];
                Int ncu = ndims[20];
                Int len = npv*ncu;
                
                DCOPY(&len, &UDG_avg[UDG_avgStart], &inc, &UDG2Write[UDG2WriteStart], &inc);
                UDG_avgStart += npv*nc;
                UDG2WriteStart += npv*ncu;
            }
            out.write( reinterpret_cast<char*>( &UDG2Write[0] ), sizeof(double) * UDG2WriteStart );
        }
        else if (writeQflag == 1) {
            Int lenUDG_avg = 0;
            for (int i=0; i<ne; i++) {
                Int elementype = mesh.elementype[ie];
                Int npv = mesh.elements[elementype].npv;
                Int nc  = ndims[19];
                lenUDG_avg += npv*nc;
            }
            out.write( reinterpret_cast<char*>( &UDG_avg[0] ), sizeof(double) * lenUDG_avg );
        }
        else {
            printf("writeQflag has invalid value in writeSol2File.\n");
            exit(-1);
        }

        // Write UH_avg:
        out.write( reinterpret_cast<char*>( &UH_avg[0] ), sizeof(double) * lenUH_avg );
    }

    // Close file
    out.close();
}

void writeTimeStepSize2File(meshstruct &mesh, string filename, Int timeStep, Int DIRKstage, double time, double dt)
{
    ofstream out(filename.c_str(), ios::out | ios::app);
    if (!out)
        cout <<"Unable to open file" << filename << endl;
    else
        out << NumberToString(timeStep) + "\t" + NumberToString(DIRKstage) + "\t" + NumberToString(time) + "\t" + NumberToString(dt) + "\n";
    out.close();
}

void writeScalarField2File(meshstruct &mesh, string filename, double* field, Int ne, Int* ndims)
{
    ofstream out(filename.c_str(), ios::out | ios::binary);
    if (!out)
        cout <<"Unable to open file" << filename << endl;
    else if (out) {
        Int lenfield = 0;
        for (int i=0; i<ne; i++) {
            Int elementype = mesh.elementype[ie];
            Int npv = mesh.elements[elementype].npv;
            lenfield  += npv;
        }
        out.write( reinterpret_cast<char*>( &field[0] ), sizeof(double) * lenfield );
    }
    out.close();
}

#endif
