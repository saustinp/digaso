#ifndef __MATHUTILITIES
#define __MATHUTILITIES

#include <vector>
#include <algorithm>
#include <numeric>

// Written by: C. Nguyen & P. Fernandez

vector<Int> uniqueiarray(vector<Int> & b)
{	
    // input: integer array b
    // output: sorted integer array a without duplicate elements
                     
    // make a new copy of b
    vector<Int> a = b;
            
	// First Sort the given range to bring duplicate
	// elements at consecutive positions
	sort(a.begin(), a.end());
  
	vector<Int>::iterator newEnd;
 
	// Override duplicate elements
	newEnd = unique(a.begin(), a.end());
 
    // remove duplicate elements
	a.erase(newEnd, a.end());
    
    return a;
}

void uniqueiarray(vector<Int> & a, vector<Int> & b)
{	
    // input: integer array b
    // output: sorted integer array a without duplicate elements
    
    // make a new copy of b
    a = b;
            
	// First Sort the given range to bring duplicate
	// elements at consecutive positions
	sort(a.begin(), a.end());
  
	vector<Int>::iterator newEnd;
 
	// Override duplicate elements
	newEnd = unique(a.begin(), a.end());
 
    // remove duplicate elements
	a.erase(newEnd, a.end());
}

struct IdxCompare
{
    const vector<Int>& target;

    IdxCompare(const vector<Int>& target): target(target) {}

    bool operator()(Int a, Int b) const { return target[a] < target[b]; }
};

vector<Int> sort_indexes(const vector<Int> &v) 
{
  // initialize original index locations
  vector<Int> idx(v.size());
  for (Int i=0; i<v.size();i++)
      idx[i] = i;  

  // sort indexes based on comparing values in v
  //sort(idx.begin(), idx.end(), [&v](Int i1, Int i2){return v[i1] < v[i2];});
  sort(idx.begin(), idx.end(), IdxCompare(v));
  
  return idx;
}

// sorting a particular column of a 2D array
void sort_columnk(vector<Int> &v, Int m, Int n, Int k) 
{
    Int i, j;
    
    // copy k-th column of v to a
    vector<Int> a(m);    
    for(i=0; i<m; i++)
        a[i] = v[k*m+i];
        
    // sort a and return indexes
    vector<Int> ind;
    ind = sort_indexes(a);
        
    // make new sorted 2D array 
    vector<Int> w(m*n);
    for (j=0; j<n; j++)
        for (i=0; i<m; i++)                    
           w[j*m+i] = v[j*m+ind[i]];
    
    // assign w back to v
    v = w;
}

void sort_rows(vector<Int> &v, Int m, Int n) 
{
    vector< vector<Int> > w(m, vector<Int>(n));
    
    Int i, j;
    for (i=0; i<m; i++) 
        for (j=0; j<n; j++)
            w[i][j] = v[m*j+i];
        
    sort(w.begin(), w.end());
    
    for (i=0; i<m; i++) 
        for (j=0; j<n; j++)
            v[m*j+i] = w[i][j];    
}

vector<Int> find(const vector<Int> &x, Int a)
{
    Int i, j = 0;
    Int n = x.size();    
    
    vector<Int> idx(n,0);        
    for (i=0; i<n; i++) 
        if (x[i] == a) {
            idx[j] = i;
            j += 1;
        }    
    idx.erase(idx.begin()+j,idx.end());    
    
    return idx;
}

vector<Int> setintersection(vector<Int> v1, vector<Int> v2)
{

    vector<Int> v3;

    sort(v1.begin(), v1.end());
    sort(v2.begin(), v2.end());

    set_intersection(v1.begin(),v1.end(),v2.begin(),v2.end(),back_inserter(v3));

    return v3;
}

vector<Int> setdifference(vector<Int> v1, vector<Int> v2)
{

    vector<Int> v3;

    sort(v1.begin(), v1.end());
    sort(v2.begin(), v2.end());

    set_difference(v1.begin(),v1.end(),v2.begin(),v2.end(),back_inserter(v3));

    return v3;
}

vector<Int> setunion(vector<Int> v1, vector<Int> v2)
{

    vector<Int> v3;

    sort(v1.begin(), v1.end());
    sort(v2.begin(), v2.end());

    set_union(v1.begin(),v1.end(),v2.begin(),v2.end(),back_inserter(v3));

    return v3;
}

bool IsSubset(vector<Int> A, vector<Int> B)
{
    vector<Int> C = setdifference(A, B);
    
    Int AinB = 0;
    if (C.empty()==1)
        AinB = 1;
    
    return AinB;
}


#endif
