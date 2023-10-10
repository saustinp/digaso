#ifndef __IOUTILITIES
#define __IOUTILITIES

template <typename T> string NumberToString ( T Number )
{
    ostringstream ss;
    ss << Number;
    return ss.str();
}

void printiarray(vector<Int> &a)
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

void printdarray(vector<double> &a, int prec)
{
    Int m = (Int) a.size();
    cout.precision(prec);
    for (Int i=0; i<m; i++)
        cout << scientific << a[i] << "   ";
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


void readcarray(ifstream &in, vector<char> &a, Int N)
{
    if (N>0) {
        a.resize(N);
        in.read( reinterpret_cast<char*>( &a[0] ), sizeof(char)*N );
    }
}

void readdarray(ifstream &in, vector<double> &a, Int N)
{
    if (N>0) {
        a.resize(N);
        in.read( reinterpret_cast<char*>( &a[0] ), sizeof(double)*N );
    }
}

void readiarray(ifstream &in, vector<Int> &a, Int N)
{
    if (N>0) {
        a.resize(N);
        in.read( reinterpret_cast<char*>( &a[0] ), sizeof(Int) * N );
    }
}

void readiarrayfromdouble(ifstream &in, vector<Int> &a, Int N)
{
    if (N>0) {
        a.resize(N);
        double read;
        for (unsigned i = 0; i < N; i++) {
            in.read( reinterpret_cast<char*>( &read ), sizeof read );
            a[i] = (Int) round(read);
        }
    }
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

#endif