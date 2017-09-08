#include "cubicalgrid.h"

void CubicalGrid::writeVTK(char * filename){

    FILE* file = fopen(filename,"w");

    fprintf(file, "# vtk DataFile Version 2.0\n\n");
    fprintf(file, "ASCII\n");

    fprintf(file, "DATASET STRUCTURED_GRID\n");
    fprintf(file, "DIMENSIONS %d %d %d\n",xRes,yRes,zRes);
    fprintf(file, "POINTS %d float\n",nVerts());


    for(int i=0; i<nVerts(); i++){
        vector<unsigned int> coords = vertexToCoords(i);
        fprintf(file, "%d %d %d\n", coords[0],coords[1],coords[2]);
    }

    fprintf(file,"POINT_DATA %d\n", nVerts());
    fprintf(file,"FIELD FieldData %d\n", fieldValues[0].size());

    string fName("originalField");

    for(unsigned int i=0; i<fieldValues[0].size(); i++){
        fprintf(file,"%s 1 %d float\n", (fName + std::to_string(i)).c_str(), fieldValues.size());
        for(unsigned int j=0; j<fieldValues.size(); j++){
            //fprintf(file,"%d ", vFiltration[i][j]); //filtration values are better for visualization
            fprintf(file,"%f ", fieldValues[j][i]); //field values are better for data analysis
        }
        fprintf(file,"\n");
    }

    fclose(file);
}
