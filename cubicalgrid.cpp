#include "cubicalgrid.h"


CubicalGrid::CubicalGrid(vector<float> val, unsigned int x, unsigned int y, unsigned int z, unsigned int fSize)
{
    fieldValues = vector<vector<float> >(val.size(), vector<float>(fSize,0));


    for(int i=0; i<val.size(); i++){
        cout << val[i] << " ";
        fieldValues[i][0]=val[i];
        cout << fieldValues[i][0] << endl;
    }

    xRes=x;
    yRes=y;
    zRes=z;
}


CubicalGrid::CubicalGrid(vector<char*> fileNames, unsigned int x, unsigned int y, unsigned int z)
{

    xRes=x;
    yRes=y;
    zRes=z;

    fieldValues = vector<vector<float> >(x*y*z, vector<float>(fileNames.size(),0));

    for(int f=0; f<fileNames.size(); f++){
        string name(fileNames[f]);
        if(name.find(".raw") !=std::string::npos)   readRAW(name,f);
        else if(name.find(".ascii") !=std::string::npos)   readASCII(name,f);
        else if(name.find(".bin") !=std::string::npos)   readBIN(name,f);
        else if(name.find(".dat") !=std::string::npos)   readDAT(name,f);
        else    cout << "Uknown input format" << endl;
    }

}

CubicalGrid::CubicalGrid(){

}

CubicalGrid::~CubicalGrid()
{

}

void CubicalGrid::readASCII(string fName, int nField){

    FILE* filew = fopen(fName.c_str(),"r");

    for(int i=0;i<nVerts();i++){
        float val;
        fscanf(filew,"%f",&val);
        fieldValues[i][nField] = val;
    }

    fclose(filew);
}

void CubicalGrid::readRAW(string fName, int nField){

    char *fileBuf=(char*) malloc(nVerts()*sizeof(char));	// Pointer to our buffered data
    FILE* filew = fopen(fName.c_str(),"rb");

    fread((char*)fileBuf, nVerts()*sizeof(char), 1, filew);

    for(int i=0;i<nVerts();i++){

        int val = fileBuf[i];
        if(val < 0) val = 128-val;
        fieldValues[i][nField] = val;
        cout << fieldValues[i][nField] << endl;
    }

    fclose(filew);
}

void CubicalGrid::readBIN(string fName, int nField){

    float *fileBuf=(float*) malloc(nVerts()*sizeof(float));	// Pointer to our buffered data
    FILE* filew = fopen(fName.c_str(),"rb");

    fread((float*)fileBuf, nVerts()*sizeof(float), 1, filew);

    for(int i=0;i<nVerts();i++){
        float val = fileBuf[i];

        //update this please (conversion from BigEndian to LittleEndian)
        char* a = (char*)&val;
        char temp1,temp2;
        temp1=a[0];
        temp2=a[1];
        a[0]=a[3];
        a[1]=a[2];
        a[3]=temp1;
        a[2]=temp2;

        fieldValues[i][nField] = val;
    }

    fclose(filew);
}

void CubicalGrid::readDAT(string fName, int nField){

    float *fileBuf=(float*) malloc(nVerts()*sizeof(float));	// Pointer to our buffered data
    FILE* filew = fopen(fName.c_str(),"rb");

    fread((float*)fileBuf, nVerts()*sizeof(float), 1, filew);

    for(int i=0;i<nVerts();i++){
        float val = fileBuf[i];
        fieldValues[i][nField] = val;
    }

    fclose(filew);
}



vector<unsigned int> CubicalGrid::vertexToCoords(unsigned int ind){
    assert(ind < nVerts());
    unsigned int z = ind/(xRes*yRes);
    unsigned int y = ind%(xRes*yRes)/xRes;
    unsigned int x = (ind%(xRes*yRes))%xRes;
    vector<unsigned int> vec = {x,y,z};
    return vec;
}

unsigned int CubicalGrid::coordsToVertex(vector<unsigned int> coords){
    return coords[0]+coords[1]*xRes+coords[2]*(xRes*yRes);
}

unsigned int CubicalGrid::cubeIndexToVindex(unsigned int cube){
    return cube-7*nVerts();
}

GCell CubicalGrid::indexToCell(unsigned int ind){

    if(ind < nVerts()){
        //vertex
        GCell v = {ind};
        return v;
    }
    else if(ind >= nVerts() && ind < 4*nVerts()){
        //edge
        unsigned int v1_ind = (ind - nVerts())/3;
        unsigned int internal_ind = (ind - nVerts())%3;

        GCell edge;
        switch(internal_ind){
            case 0:
                edge = {v1_ind,v1_ind+1};
                break;
            case 1:
                edge = {v1_ind,v1_ind+xRes};
                break;
            case 2:
                edge = {v1_ind,v1_ind+(xRes*yRes)};
                break;
            default:
                cout << "Error! wrong internal index" << endl;
                break;
        }
        return edge;

    }
    else if(ind >= 4*nVerts() && ind < 7*nVerts()){
        //face
        unsigned int v1_ind = (ind - 4*nVerts())/3;
        unsigned int internal_ind = (ind - 4*nVerts())%3;

        GCell face;
        switch(internal_ind){
            case 0:
                face = {v1_ind,v1_ind+1,v1_ind+(xRes*yRes),v1_ind+(xRes*yRes)+1};
                break;
            case 1:
                face = {v1_ind,v1_ind+1,v1_ind+xRes,v1_ind+xRes+1};
                break;
            case 2:
                face = {v1_ind,v1_ind+xRes,v1_ind+(xRes*yRes),v1_ind+(xRes*yRes)+xRes};
                break;
            default:
                cout << "Error! wrong internal index" << endl;
                break;
        }
        return face;
    }
    else{
        //voxel
        assert(ind < 8*nVerts());
        GCell cube;
        unsigned int v1_ind = (ind - 7*nVerts());
        cube = {v1_ind, v1_ind+1, v1_ind+xRes, v1_ind+xRes+1, v1_ind+(xRes*yRes), v1_ind+(xRes*yRes)+1, v1_ind+(xRes*yRes)+xRes, v1_ind+(xRes*yRes)+xRes+1};
        return cube;
    }
}

unsigned int CubicalGrid::cellToIndex(GCell cell){

    unsigned int internal_i=0;

    switch(cell.size()){
        case 1:
            return cell[0];

        case 2:
            if(cell[1]-cell[0] == 1){
                internal_i=0;
            }
            else if(cell[1]-cell[0] == xRes){
                internal_i=1;
            }
            else if(cell[1]-cell[0] == (xRes*yRes)){
                internal_i=2;
            }
            else{
                cout << "Erro! wrong vertex index inside edge" << endl;
            }

            return nVerts()+3*cell[0]+internal_i;

        case 4:
            if(cell[3]-cell[0] == (xRes*yRes)+1){
                internal_i=0;
            }
            else if(cell[3]-cell[0] == xRes+1){
                internal_i=1;
            }
            else if(cell[3]-cell[0] == (xRes*yRes)+xRes){
                internal_i=2;
            }
            else{
                cout << "Erro! wrong vertex index inside face" << endl;
            }
            return 4*nVerts()+3*cell[0]+internal_i;

        case 8:
            return 7*nVerts()+cell[0];

        default:
            cout << "Erro! wrong cell" << endl;
            return 0;
    }
}

unsigned short CubicalGrid::internalIndex(GCell cellUp, GCell cellDown){

    switch(cellUp.size()){
    case 2:
        assert(cellDown.size()==1);
        if(cellDown[0]==cellUp[0]) return 0;
        if(cellDown[0]==cellUp[1]) return 1;
        return 6;
    case 4:
        assert(cellDown.size()==2);
        if(cellDown[1]==cellUp[1]) return 0;
        if(cellDown[1]==cellUp[2]) return 1;
        if(cellDown[0]==cellUp[1]) return 2;
        if(cellDown[0]==cellUp[2]) return 3;
        return 7;
    case 8:
        assert(cellDown.size()==4);
        if(cellDown[2]==cellUp[2]) return 0;
        if(cellDown[3]==cellUp[5]) return 1;
        if(cellDown[3]==cellUp[6]) return 2;
        if(cellDown[0]==cellUp[1]) return 3;
        if(cellDown[0]==cellUp[2]) return 4;
        if(cellDown[0]==cellUp[4]) return 5;
        return 8;
    default:
        cout << "Wrong cell pair" << endl;
        return 9;
    }
}

GCell CubicalGrid::internalIndexToCell(GCell cellUp, unsigned short index){
    switch(cellUp.size()){
    case 2:
    {
        GCell vertex = {cellUp[index]};
        return vertex;
    }
    case 4:
    {
        GCell edge;
        switch(index){
            case 0: edge = {cellUp[0],cellUp[1]}; break;
            case 1: edge = {cellUp[0],cellUp[2]}; break;
            case 2: edge = {cellUp[1],cellUp[3]}; break;
            case 3: edge = {cellUp[2],cellUp[3]}; break;
            default: cout << "Wrong internal index edge!" << endl; break;
        }
        return edge;
    }
    case 8:
    {
        GCell face;
        switch(index){
            case 0: face = {cellUp[0],cellUp[1],cellUp[2],cellUp[3]}; break;
            case 1: face = {cellUp[0],cellUp[1],cellUp[4],cellUp[5]}; break;
            case 2: face = {cellUp[0],cellUp[2],cellUp[4],cellUp[6]}; break;
            case 3: face = {cellUp[1],cellUp[3],cellUp[5],cellUp[7]}; break;
            case 4: face = {cellUp[2],cellUp[3],cellUp[6],cellUp[7]}; break;
            case 5: face = {cellUp[4],cellUp[5],cellUp[6],cellUp[7]}; break;
            default: cout << "Wrong internal index face!" << endl; break;
        }
        return face;
    }
    default:
        cout << "Wrong internal index!" << endl;
        return GCell();
    }
}

void CubicalGrid::cubeFaces(GCell cube, list<GCell> &faces){

    GCell face;
    face = {cube[0],cube[1],cube[2],cube[3]}; faces.push_back(face);
    face = {cube[0],cube[1],cube[4],cube[5]}; faces.push_back(face);
    face = {cube[0],cube[2],cube[4],cube[6]}; faces.push_back(face);
    face = {cube[1],cube[3],cube[5],cube[7]}; faces.push_back(face);
    face = {cube[2],cube[3],cube[6],cube[7]}; faces.push_back(face);
    face = {cube[4],cube[5],cube[6],cube[7]}; faces.push_back(face);
}

void CubicalGrid::faceEdges(GCell face, list<GCell> &edges){

    GCell edge;
    edge = {face[0],face[1]}; edges.push_back(edge);
    edge = {face[0],face[2]}; edges.push_back(edge);
    edge = {face[1],face[3]}; edges.push_back(edge);
    edge = {face[2],face[3]}; edges.push_back(edge);
}

void CubicalGrid::edgeVertices(GCell edge, list<GCell> &vertices){

    GCell vert;
    vert = {edge[0]}; vertices.push_back(vert);
    vert = {edge[1]}; vertices.push_back(vert);
}


void CubicalGrid::faceCubes(GCell face,list<GCell>& cubes){

    GCell cube;
    vector<unsigned int> coords1 = vertexToCoords(face[1]);
    vector<unsigned int> coords2 = vertexToCoords(face[2]);

    if(abs((int)coords1[0]-(int)coords2[0])==0){
        //chage x
        if(coords1[0]>0){
            cube = {face[0]-1, face[0], face[1]-1, face[1], face[2]-1, face[2], face[3]-1, face[3]}; cubes.push_back(cube);
        }
        if(coords1[0]<xRes-1){
            cube = {face[0], face[0]+1, face[1], face[1]+1, face[2], face[2]+1, face[3], face[3]+1}; cubes.push_back(cube);
        }
    }
    else if(abs((int)coords1[1]-(int)coords2[1])==0){
        //chage y
        if(coords1[1]>0){
            cube = {face[0]-xRes, face[1]-xRes, face[0], face[1], face[2]-xRes, face[3]-xRes, face[2], face[3]}; cubes.push_back(cube);
        }
        if(coords1[1]<yRes-1){
            cube = {face[0], face[1], face[0]+xRes, face[1]+xRes, face[2], face[3], face[2]+xRes, face[3]+xRes}; cubes.push_back(cube);
        }
    }
    else{
        assert(abs((int)coords1[2]-(int)coords2[2])==0);
        //chage z
        if(coords1[2]>0){
            cube = {face[0]-(xRes*yRes), face[1]-(xRes*yRes), face[2]-(xRes*yRes), face[3]-(xRes*yRes), face[0], face[1], face[2], face[3]}; cubes.push_back(cube);
        }
        if(coords1[2]<zRes-1){
            cube = {face[0], face[1], face[2], face[3], face[0]+(xRes*yRes), face[1]+(xRes*yRes), face[2]+(xRes*yRes), face[3]+(xRes*yRes)}; cubes.push_back(cube);
        }
    }
}


void CubicalGrid::edgeFaces(GCell edge,list<GCell>& faces){

    GCell face;
    vector<unsigned int> coords1 = vertexToCoords(edge[0]);
    vector<unsigned int> coords2 = vertexToCoords(edge[1]);

    if(abs((int)coords1[0]-(int)coords2[0])==1){
        //x is changing
        if(coords1[1]>0){
            face = {edge[0]-xRes,edge[1]-xRes,edge[0],edge[1]}; faces.push_back(face);
        }
        if(coords1[1]<yRes-1){
            face = {edge[0],edge[1],edge[0]+xRes,edge[1]+xRes}; faces.push_back(face);
        }
        if(coords1[2]>0){
            face = {edge[0]-(xRes*yRes),edge[1]-(xRes*yRes),edge[0],edge[1]}; faces.push_back(face);
        }
        if(coords1[2]<zRes-1){
            face = {edge[0],edge[1],edge[0]+(xRes*yRes),edge[1]+(xRes*yRes)}; faces.push_back(face);
        }
    }
    else if(abs((int)coords1[1]-(int)coords2[1])==1){
        //y is changing
        if(coords1[0]>0){
            face = {edge[0]-1,edge[0],edge[1]-1,edge[1]}; faces.push_back(face);
        }
        if(coords1[0]<xRes-1){
            face = {edge[0],edge[0]+1,edge[1],edge[1]+1}; faces.push_back(face);
        }
        if(coords1[2]>0){
            face = {edge[0]-(xRes*yRes),edge[1]-(xRes*yRes),edge[0],edge[1]}; faces.push_back(face);
        }
        if(coords1[2]<zRes-1){
            face = {edge[0],edge[1],edge[0]+(xRes*yRes),edge[1]+(xRes*yRes)}; faces.push_back(face);
        }
    }
    else{
        assert(abs((int)coords1[2]-(int)coords2[2])==1);
        //z is changing
        if(coords1[0]>0){
            face = {edge[0]-1,edge[0],edge[1]-1,edge[1]}; faces.push_back(face);
        }
        if(coords1[0]<xRes-1){
            face = {edge[0],edge[0]+1,edge[1],edge[1]+1}; faces.push_back(face);
        }
        if(coords1[1]>0){
            face = {edge[0]-xRes,edge[0],edge[1]-xRes,edge[1]}; faces.push_back(face);
        }
        if(coords1[1]<yRes-1){
            face = {edge[0],edge[0]+xRes,edge[1],edge[1]+xRes}; faces.push_back(face);
        }
    }
}

void CubicalGrid::vertexEdges(GCell vertex, list<GCell> &edges){

    GCell edge;
    vector<unsigned int> coords1 = vertexToCoords(vertex[0]);
    if(coords1[0]>0){
        edge={vertex[0]-1,vertex[0]}; edges.push_back(edge);
    }
    if(coords1[0]<xRes-1){
        edge={vertex[0],vertex[0]+1}; edges.push_back(edge);
    }
    if(coords1[1]>0){
        edge={vertex[0]-xRes,vertex[0]}; edges.push_back(edge);
    }
    if(coords1[1]<yRes-1){
        edge={vertex[0],vertex[0]+xRes}; edges.push_back(edge);
    }
    if(coords1[2]>0){
        edge={vertex[0]-(xRes*yRes),vertex[0]}; edges.push_back(edge);
    }
    if(coords1[2]<zRes-1){
        edge={vertex[0],vertex[0]+(xRes*yRes)}; edges.push_back(edge);
    }
}

void CubicalGrid::vertexVertices(GCell vertex, list<GCell> &vertices){

    GCell v;
    vector<unsigned int> coords1 = vertexToCoords(vertex[0]);
    if(coords1[0]>0){
        v={vertex[0]-1}; vertices.push_back(v);
    }
    if(coords1[0]<xRes-1){
        v={vertex[0]+1}; vertices.push_back(v);
    }
    if(coords1[1]>0){
        v={vertex[0]-xRes}; vertices.push_back(v);
    }
    if(coords1[1]<yRes-1){
        v={vertex[0]+xRes}; vertices.push_back(v);
    }
    if(coords1[2]>0){
        v={vertex[0]-(xRes*yRes)}; vertices.push_back(v);
    }
    if(coords1[2]<zRes-1){
        v={vertex[0]+(xRes*yRes)}; vertices.push_back(v);
    }
}


void CubicalGrid::vertexStar(GCell v, list<GCell> &cells){

    //rewrite hardcoded for optimization

    set<GCell> newFaces;
    vertexEdges(v, cells);

    list<GCell> newWithRep;
    for(auto e : cells){
        edgeFaces(e,newWithRep);
        newFaces.insert(newWithRep.begin(),newWithRep.end());
        newWithRep.clear();
    }
    cells.insert(cells.end(),newFaces.begin(),newFaces.end());

    set<GCell> newCubes;
    for(auto f : newFaces){
        faceCubes(f,newWithRep);
        newCubes.insert(newWithRep.begin(),newWithRep.end());
        newWithRep.clear();
    }
    cells.insert(cells.end(),newCubes.begin(),newCubes.end());
}

void CubicalGrid::immediateBoundary(GCell cell, list<GCell> &cells){

    if(cell.size() == 2)
        edgeVertices(cell,cells);
    else if(cell.size() == 4)
        faceEdges(cell,cells);
    else if(cell.size() == 8)
        cubeFaces(cell,cells);
    else{
        cout << "No boundary cells on a vertex" << endl;
    }
}

void CubicalGrid::immediateCoboundary(GCell cell, list<GCell> &cells){

    if(cell.size() == 1)
        vertexEdges(cell,cells);
    else if(cell.size() == 2)
        edgeFaces(cell,cells);
    else if(cell.size() == 4)
        faceCubes(cell,cells);
    else{
        cout << "No coboundary cells on a cube" << endl;
    }
}
