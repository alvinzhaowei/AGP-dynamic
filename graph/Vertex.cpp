#include "Vertex.h"

using namespace dynagp;

static const int adj_hash_initialization_size = 11;

Vertex::Vertex(const int &_id) {
    id = _id;
    adjacentList.reserve(adj_hash_initialization_size);
    neighborIDAdjacentIndexMap.reserve(adj_hash_initialization_size);
}

void Vertex::insertNeighbor(const int &_neighborID) {
    adjacentList.push_back(_neighborID);
    neighborIDAdjacentIndexMap[_neighborID] = adjacentList.size() - 1;
}


void Vertex::deleteNeighbor(const int _neighborID) {
    int _index = neighborIDAdjacentIndexMap.at(_neighborID);
    neighborIDAdjacentIndexMap.erase(_neighborID);
    int length = adjacentList.size();
    if (length == _index + 1) {
        adjacentList.pop_back();
    } else {

        adjacentList[_index] = adjacentList[length - 1];
        adjacentList.pop_back();
        neighborIDAdjacentIndexMap.at(Vertex::adjacentList[_index]) = _index;
    }
}


void Vertex::assignDS(DSBucket* ptr){
    DsBucketPtr = ptr;
}

DSBucket *Vertex::getDS() {
    return DsBucketPtr;
}
