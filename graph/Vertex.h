#ifndef DYNAGP_VERTEX_H
#define DYNAGP_VERTEX_H

#include "../MyLib/MyVector.h"
#include "../Tessil_robin_map/robin_map.h"
#include "../DS/DS.h"

#define hash_map tsl::robin_map

namespace dynagp {

    class Vertex {

    public:
        int id;
    protected:

        MyVector<int> adjacentList;

        hash_map<int, int> neighborIDAdjacentIndexMap;

        // A pointer to DS
        DSBucket* DsBucketPtr;

        MyVector<int>* NonEmptyIndex;

    public:

        hash_map<int, pair<int, int>> neighborElementIndexMap;

        inline int getNeighborID(const int &_index) const {
            return adjacentList[_index];
        }


        /**
         *
         * @return return current vertex's degree
         */
        inline int getDegree() const {
            return adjacentList.size();
        };


        Vertex(const int &_id);

        /**
         * Insert a new neighbor into adjacent list
         * @param _neighborID
         */

        void insertNeighbor(const int &_neighborID);


        /**
         * Delete a neighbor at give index
         * @param _index
         */
        void deleteNeighbor(const int _neighborID);

        inline int *getAdjacentList() {
            return Vertex::adjacentList.get_list();
        };

        /**
         * @param _neighborID
         * @return the index of the given neighbor ID in the adjacent list of the current vertex.
         */
        inline const int getAdjacentIndex(const int &_neighborID) const {
            auto it = neighborIDAdjacentIndexMap.find(_neighborID);
            return it == neighborIDAdjacentIndexMap.end() ? -1 : it->second;
        }

        /**
         * @param _neighbor
         * @return
         */
        inline bool has_neighbor(const int &_neighbor) const {
            return neighborIDAdjacentIndexMap.find(_neighbor) != neighborIDAdjacentIndexMap.end()
                   || this->id == _neighbor;
        }

        void assignDS(DSBucket* ptr);

        void assignNoEmptyIndex(MyVector<int>* ptr) {
            NonEmptyIndex = ptr;
        };

        DSBucket* getDS();

        MyVector<int>* getNoEmtInx() {
            return NonEmptyIndex;
        };

        pair<int, int> getElementIndex(const int &_neighborID) const {
            auto it = neighborElementIndexMap.find(_neighborID);
            if (it != neighborElementIndexMap.end()) {
                return it->second;
            } else {
                return {-1, -1};
            }
        }
    };

}
#endif //DYNAGP_VERTEX_H
