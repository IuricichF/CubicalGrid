#ifndef CUBICALGRID_H
#define CUBICALGRID_H

#include <vector>
#include <list>
#include <cmath>
#include <set>
#include <map>
#include <queue>

#include <boost/function.hpp>
#include <boost/bind.hpp>

#include <assert.h>
#include <iostream>
#include <fstream>

#include "Timer.h"
#include "Usage.h"

using namespace std;

typedef vector<unsigned int> GCell;
typedef set<GCell, boost::function<bool(const GCell &, const GCell &)> > GCellSet;

inline bool operator== (const GCell& lhs, const GCell& rhs){
    if(lhs.size() != rhs.size())
        return false;
    else{
        for(int i=0; i<rhs.size(); i++)
            if(rhs[i] != lhs[i])
                return false;
    }
    return true;
}

class CubicalGrid
{
private:
    vector<vector<float> > fieldValues;

    unsigned int xRes;
    unsigned int yRes;
    unsigned int zRes;

public:
    CubicalGrid(vector<float>,unsigned int, unsigned int, unsigned int, unsigned int);
    CubicalGrid(vector<char*> fileNames, unsigned int, unsigned int, unsigned int);
    CubicalGrid();
    ~CubicalGrid();


    //Functions for cubical grid input/output
    void readRAW(string fName, int nField);
    void readBIN(string fName, int nField);
    void readDAT(string fName, int nField);
    void readASCII(string fName, int nField);

    void writeVTK(char*);

    //Basic functions for geometry
    void computeCentroid(unsigned iCell, vector<float> &coords);
    void resizeGrid();

    inline float getFieldValue(unsigned int i_vertex, unsigned int i_function){return fieldValues[i_vertex][i_function];}
    inline unsigned int nFields(){return fieldValues[0].size();}

    //Basic functions for topological relations
    inline unsigned int nVerts(){return xRes*yRes*zRes;}
    inline unsigned int nCubes(){return nVerts()-(xRes*yRes)-((xRes-1)*(zRes-1))-((yRes-1)*(zRes-1))-(zRes-1);}

    GCell indexToCell(unsigned int ind);
    unsigned int cellToIndex(GCell);
    unsigned int cubeIndexToVindex(unsigned int cube);

    vector<unsigned int> vertexToCoords(unsigned int v);
    unsigned int coordsToVertex(vector<unsigned int> coords);

    void faceCubes(GCell face,list<GCell>& cubes);
    void edgeFaces(GCell edge,list<GCell>& faces);
    void vertexEdges(GCell vertex,list<GCell>& edges);
    void vertexVertices(GCell vertex, list<GCell> &vertices);

    void cubeFaces(GCell cube,list<GCell>& faces);
    void faceEdges(GCell face,list<GCell>& edges);
    void edgeVertices(GCell edge,list<GCell>& vertices);

    void vertexStar(GCell v, list<GCell>& cells);
    void immediateBoundary(GCell cell, list<GCell>& cells);
    void immediateCoboundary(GCell cell, list<GCell>& cells);

    unsigned short internalIndex(GCell cellUp, GCell cellDown);
    GCell internalIndexToCell(GCell cellUp, unsigned short index);

};

#endif // CUBICALGRID_H
