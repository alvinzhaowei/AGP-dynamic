#ifndef DYNAGP_GRAPH_H
#define DYNAGP_GRAPH_H

#include "Vertex.h"
#include "../MyLib/MyVector.h"
#include "../DS/DS.h"
#include <cmath>
#include <random>

using namespace std;

class Graph {
protected:

    // The list of all the vertices.
    MyVector<dynagp::Vertex *> vList;
    int vertex_number;
    MyVector<int> T;


public:
    /*
     *  The constructor of Graph.
     */
    Graph(MyVector<dynagp::Vertex *> &_vList);
    /*
 *  Insert an edge to the graph.
 */
    int insertEdge(int _vID1, int _vID2);

    /*
     *  Delete an edge from the graph.
     */
    int removeEdge(int _vID1, int _vID2);

    void naive_initialize(double* pi, double* x);

    void construct_DS(float a, float b);

    void construct_DS_static(float a, float b);

    double* query(double x[], float a, float b, int L, double* w, double eps, char type);

    int update_DS(int _vID1, int _vID2, float a, float b);

    int naive_insert(int _vID1, int _vID2, float a, float b);

    int naive_delete(int _vID1, int _vID2, float a, float b);

    void DAGP_initialize();

    int DAGP_insert(int _vID1, int _vID2, float a, float b);

    int DAGP_delete(int _vID1, int _vID2, float a, float b);


protected:
    /**
     *
     * @param _id
     * @return create new vertex with vertex factory and return its pointer.
     */
    dynagp::Vertex* createVertex(const int &_id){
        dynagp::Vertex *newVertex = new dynagp::Vertex(_id);
        vList[_id - 1] = newVertex;
        return newVertex;
    }

};


#endif //DYNAGP_GRAPH_H
