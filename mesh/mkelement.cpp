#ifndef __MKELEMENT
#define __MKELEMENT

elementinfo mkeleminfo(Int dim, Int elemtype)
{    
    Int nfe, nle, nve, nqf, ntf;
    vector< Int > facetype;
    vector< vector< Int > > face;
    
    switch (dim) {
        case 1: // line element
            nfe = 2;
            nle = 0;
            nve = 2;
            nqf = 0;
            ntf = 0;
            facetype.resize(nfe,0);
            face.resize(nfe);
            face[0].resize(1);
            face[1].resize(1);
            face[0][0] = 0;
            face[1][0] = 1;
            break;
        case 2:            
            if (elemtype==0) { // triangular
                nfe = dim+1;  
                face.resize(nfe);
                
                // [[1,2];[2 0]; [0,1]]
                face[0].resize(2);
                face[0][0] = 1;
                face[0][1] = 2;
                
                face[1].resize(2);
                face[1][0] = 2;
                face[1][1] = 0;
                
                face[2].resize(2);
                face[2][0] = 0;
                face[2][1] = 1;     
            }
            else if (elemtype==1) { // quadrilateral
                nfe = 2*dim;   
                face.resize(nfe);
                
                // [[1,2];[2,3];[3,4];[4,1]] - 1;
                face[0].resize(2);
                face[0][0] = 0;
                face[0][1] = 1;
                
                face[1].resize(2);
                face[1][0] = 1;
                face[1][1] = 2;
                
                face[2].resize(2);
                face[2][0] = 2;
                face[2][1] = 3;                            
                
                face[3].resize(2);
                face[3][0] = 3;
                face[3][1] = 0;    
            }
            nle = nfe;
            nve = nfe;
            nqf = 0;
            ntf = 0;
            facetype.resize(nfe,0);
            break;
        case 3:
            if (elemtype==0) { // tet
                nfe = dim+1;
                nve = 4;
                nle = 6;
                ntf = 4;
                nqf = 0;
                
                facetype.resize(nfe,0);                
                face.resize(nfe);
                        
                // [[2,3,4];[1,4,3];[1,2,4];[1,3,2]] - 1
                face[0].resize(3);
                face[0][0] = 1;
                face[0][1] = 2;
                face[0][2] = 3;
                
                face[1].resize(3);
                face[1][0] = 0;
                face[1][1] = 3;
                face[1][2] = 2;
                
                face[2].resize(3);
                face[2][0] = 0;
                face[2][1] = 1;
                face[2][2] = 3;
                
                face[3].resize(3);
                face[3][0] = 0;
                face[3][1] = 2;
                face[3][2] = 1;                
            }
            else if (elemtype==1) { // hex
                nfe = 2*dim;       
                nve = 8;
                nle = 12;
                nqf = nfe;
                ntf = 0;
                
                facetype.resize(nfe,1);                
                face.resize(nfe);
                
                //face=[[1,4,3,2];[5,6,7,8];[1,2,6,5];[3,4,8,7];[2,3,7,6];[4,1,5,8]] - 1;                
                face[0].resize(4);
                face[0][0] = 0;
                face[0][1] = 3;
                face[0][2] = 2;
                face[0][3] = 1;
                
                face[1].resize(4);
                face[1][0] = 4;
                face[1][1] = 5;
                face[1][2] = 6;
                face[1][3] = 7;
                
                face[2].resize(4);
                face[2][0] = 0;
                face[2][1] = 1;
                face[2][2] = 5;
                face[2][3] = 4;
                
                face[3].resize(4);
                face[3][0] = 2;
                face[3][1] = 3;
                face[3][2] = 7;                
                face[3][3] = 6;                
                
                face[4].resize(4);
                face[4][0] = 1;
                face[4][1] = 2;
                face[4][2] = 6;                
                face[4][3] = 5;         
                
                face[5].resize(4);
                face[5][0] = 3;
                face[5][1] = 0;
                face[5][2] = 4;                
                face[5][3] = 7;         
                //face=[[1,4,3,2];[5,6,7,8];[1,2,6,5];[3,4,8,7];[2,3,7,6];[4,1,5,8]] - 1;                
            }
            else if (elemtype==2) { // prism
                nfe = 5;                          
                nve = 6;
                nle = 9;
                nqf = 3;
                ntf = 2;
                
                facetype.resize(nfe,0);                
                facetype[0] = 0;
                facetype[1] = 0;
                facetype[2] = 1;
                facetype[3] = 1;
                facetype[4] = 1;
                
                face.resize(nfe);                
                //face=[[0,2,1];[3,4,5];[1,2,5,4];[2,0,3,5];[0,1,4,3]];                
                face[0].resize(3);
                face[0][0] = 0;
                face[0][1] = 2;
                face[0][2] = 1;                
                
                face[1].resize(3);
                face[1][0] = 3;
                face[1][1] = 4;
                face[1][2] = 5;                
                
                face[3].resize(4);
                face[3][0] = 1;
                face[3][1] = 2;
                face[3][2] = 5;                
                face[3][3] = 4;                
                
                face[2].resize(4);
                face[2][0] = 2;
                face[2][1] = 0;
                face[2][2] = 3;
                face[2][3] = 5;                                
                
                face[4].resize(4);
                face[4][0] = 0;
                face[4][1] = 1;
                face[4][2] = 4;                
                face[4][3] = 3;                         
                //face=[[0,2,1];[3,4,5];[1,2,5,4];[2,0,3,5];[0,1,4,3]];                                
            }
            else if (elemtype==3) { // pyramid
                nfe = 5;
                nve = 5;
                nle = 8;
                nqf = 1;
                ntf = 4;
                
                facetype.resize(nfe,0);                
                facetype[0] = 1;
                facetype[1] = 0;
                facetype[2] = 0;
                facetype[3] = 0;
                facetype[4] = 0;
                
                face.resize(nfe);
                //face=[[0,3,2,1];[0,1,4];[1,2,4];[2,3,4];[3,0,4]];                
                face[0].resize(4);
                face[0][0] = 0;
                face[0][1] = 3;
                face[0][2] = 2;
                face[0][3] = 1;
                
                face[1].resize(3);
                face[1][0] = 0;
                face[1][1] = 1;
                face[1][2] = 4;
                
                face[2].resize(3);
                face[2][0] = 1;
                face[2][1] = 2;
                face[2][2] = 4;
                
                face[3].resize(3);
                face[3][0] = 2;
                face[3][1] = 3;
                face[3][2] = 4;                
                
                face[4].resize(3);
                face[4][0] = 3;
                face[4][1] = 0;
                face[4][2] = 4;                
            }
            break;
        default:
            error("Only can handle dim=1, dim=2 or dim=3\n");
    }
    
    elementinfo eleminfo;
    eleminfo.nfe = nfe;
    eleminfo.nle = nle;
    eleminfo.nve = nve;
    eleminfo.nqf = nqf;
    eleminfo.ntf = ntf;
    eleminfo.facetype = facetype;
    eleminfo.face = face;
    return eleminfo;
}

elementstruct mkelementstruct(Int dim, Int elemtype, vector<Int> &porder, Int nodetype)
{        
    elementstruct element;
    
    elementinfo eleminfo = mkeleminfo(dim,elemtype);    
    element.dim = eleminfo.dim;     
    element.elemtype = eleminfo.elemtype; // type         
    element.nfe = eleminfo.nfe; // number of faces    
    element.nle = eleminfo.nle; // number of edges 
    element.nve = eleminfo.nve; // number of vertices     
    element.facetype = eleminfo.facetype; // 0 triangular, 1 quadrilateral 
    element.face = eleminfo.face; // face connectivity     
    
    element.nodetype = nodetype;
    element.porder = porder;    
    
    masternodes(element.plocvl, element.plocfc, element.perm, element.npe, element.npf, porder, dim, elemtype, nodetype);
        
    return element;
}

vector< elementstruct > mkelementstructs(Int dim, vector<Int> &porder, Int nodetype)
{
    Int nelemtype;
    switch (dim) {
        case 1: // line element            
            nelemtype = 1;
            break;            
        case 2:
            nelemtype = 2;
            break;
        case 3:
            nelemtype = 4;
            break;
        default:
            error("Only can handle dim=1, dim=2 or dim=3\n");
    }
    
    vector< elementstruct > elements(nelemtype, elementstruct());
    
    //mkelemstruct(Int dim, Int elemtype, Int porder, Int nodetype)
    for (Int i=0; i<nelemtype; i++)
        elements[i] = mkelementstruct(dim, i, porder, nodetype);
    
    return elements;    
}

elementstruct mkelementstruct(Int dim, Int elemtype, Int porder, Int nodetype)
{   
    vector<Int> p(dim,0);
    for (Int i=0; i<dim; i++)
        p[i] = porder;    
    elementstruct element = mkelementstruct(dim, elemtype, p, nodetype);
    
    return element;
}

vector< elementstruct > mkelementstructs(Int dim, Int porder, Int nodetype)
{ 
    vector<Int> p(dim,0);
    for (Int i=0; i<dim; i++)
        p[i] = porder;    
    vector< elementstruct > elements = mkelementstructs(dim, p, nodetype);
    
    return elements;
}

#endif
