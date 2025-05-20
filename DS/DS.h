#ifndef DAGP_DS_H
#define DAGP_DS_H

#include "../MyLib/MyVector.h"

using namespace std;

/*
 *  This is the element class of the buckets for DS.
 */
class DSBucketElement {
    friend class DSBucket;

protected:
    // records the IDs of current vertex's neighbor.
    int neighborID;
    // store the element index
    unsigned int element_index;

public:
    DSBucketElement(int id) {
        neighborID = id;
//        p_ = 0;
        element_index = 0;
    }

//    inline const int &get_p() const {
//        return p_;
//    }

    inline const int &get_neighbor_id() {
        return neighborID;
    };

//    inline void update_p(float p) {
//        this->p_ = p;
//    }

    inline const unsigned int &get_element_index() const {
        return element_index;
    }

};


/*
 *  The class of buckets for DS.
 */
class DSBucket {
public:

    MyVector<MyVector<DSBucketElement *>> buckList;


    //constructor
    DSBucket() {

    }

    static MyVector<DSBucketElement *> EmptyBucket;

    int listSize();

    int sizeByIndex(int i);

    void extendListSize(int size);

    void shrinkToFit();

    int InsertNewELement(int i, DSBucketElement *e);

    bool CheckEmptyByIndex(int index);

    // delete element by bucket index and element index
    void DeleteElement(int bucketIndex, int elementIndex);

    DSBucketElement *getElement(int i, int j) const;

    int getCnt(int i) const;

    ~DSBucket();

    void ReleaseSpace();

};

#endif //DAGP_DS_H
