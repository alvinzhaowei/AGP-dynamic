#include "Graph.h"
#include <random>
#include "../MyLib/MyTimer.h"


Graph::Graph(MyVector<dynagp::Vertex *> &_vList) {
    vList.swap(_vList);
    vertex_number = (int) vList.size();
}

int Graph::insertEdge(int _vID1, int _vID2) {

    auto *v1 = (dynagp::Vertex *) vList[_vID1 - 1];
    auto *v2 = (dynagp::Vertex *) vList[_vID2 - 1];
    // check if these vertices already exist. If not, create new vertex/vertices
    // with vertex factory and insert it/them into unordered map.
    if (v1 == NULL) {
        v1 = (dynagp::Vertex *) createVertex(_vID1);
    }
    if (v2 == NULL) {
        v2 = (dynagp::Vertex *) createVertex(_vID2);
    }
    v1->insertNeighbor(_vID2);
    v2->insertNeighbor(_vID1);

    return 0;
}

int Graph::removeEdge(int _vID1, int _vID2) {
    auto *v1 = (dynagp::Vertex *) vList[_vID1 - 1];
    auto *v2 = (dynagp::Vertex *) vList[_vID2 - 1];

    if (v1 == NULL || v2 == NULL)
        return 1;
    v1->deleteNeighbor(_vID2);
    v2->deleteNeighbor(_vID1);


    return 0;
}

void Graph::naive_initialize(double *r, double *x) {
    // Naive method
    for (int i = 0; i < vertex_number; i++) {
        r[i] = x[i];
    }
}

void Graph::construct_DS(float a, float b) {
    //construct subset sampling data structure
    int *u_adjacentList;
    int d_u;
    for (int i = 0; i < vertex_number; i++) {
        auto *u = (dynagp::Vertex *) vList[i];
        d_u = u->getDegree();
        u_adjacentList = u->getAdjacentList();
        auto *buckets = new DSBucket();
        u->assignDS(buckets);
        for (int j = 0; j < d_u; j++) {
            const int &neighborID = u_adjacentList[j];
            auto *neighbor_v = (dynagp::Vertex *) vList[neighborID - 1];
            int d_v = neighbor_v->getDegree();
            int ind = floor(log2(d_v));
            auto *dse = new DSBucketElement(neighborID);
            buckets->InsertNewELement(ind, dse);
            u->neighborElementIndexMap[neighborID] = {ind, dse->get_element_index()};
        }
    }
}

void Graph::construct_DS_static(float a, float b) {
    //construct subset sampling data structure
    int *u_adjacentList;
    int d_u;
    for (int i = 0; i < vertex_number; i++) {
        auto *u = (dynagp::Vertex *) vList[i];
        d_u = u->getDegree();
        u_adjacentList = u->getAdjacentList();
        auto *buckets = new DSBucket();
        u->assignDS(buckets);
        for (int j = 0; j < d_u; j++) {
            const int &neighborID = u_adjacentList[j];
            auto *neighbor_v = (dynagp::Vertex *) vList[neighborID - 1];
            int d_v = neighbor_v->getDegree();
            int ind = floor(log2(d_v));
            auto *dse = new DSBucketElement(neighborID);
            buckets->InsertNewELement(ind, dse);
            u->neighborElementIndexMap[neighborID] = {ind, dse->get_element_index()};
        }
        auto *intIndex = new MyVector<int>;
        for (int j = 0; j < buckets->listSize(); j++) {
            if (buckets->sizeByIndex(j) != 0) {
                intIndex->push_back(j);
            }
        }
        u->assignNoEmptyIndex(intIndex);
    }
}

int generateRandomInt(int min, int max) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dis(min, max);

    return dis(gen);
}

double *Graph::query(double *x, float a, float b, int L, double *w, double eps, char type) {
    // type = 'N': naive method;
    auto *pi = new double[vertex_number];
    auto *r = new double[vertex_number];
    auto *r_ = new double[vertex_number];
    auto *q = new double[vertex_number];
    auto *y = new double[L + 1];

    //calculate Y
    y[0] = 1.0;
    for (int i = 1; i < L + 1; i++) {
        y[i] = y[i - 1] - w[i - 1];
    }


    for (int i = 0; i < vertex_number; i++) {
        r_[i] = 0;
        q[i] = 0;
        pi[i] = 0;
    }
    double start;
    int *u_adjacentList;
    int d_u;
    double value;
    if (type == 'N') {
        naive_initialize(r, x);
        start = getCurrentTime();
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < vertex_number; j++) {
                if (r[j] == 0) {
                    continue;
                } else {
                    auto *u = (dynagp::Vertex *) vList[j];
                    d_u = u->getDegree();
                    value = y[i + 1] * r[j] / y[i] / pow(d_u, b);
                    u_adjacentList = u->getAdjacentList();
                    for (int s = 0; s < d_u; s++) {
                        const int &neighborID = u_adjacentList[s];
                        auto *neighbor_v = (dynagp::Vertex *) vList[neighborID - 1];
                        r_[neighborID] = r_[neighborID] + value / pow(neighbor_v->getDegree(), a);
                    }
                    q[j] = q[j] + w[i] * r[j] / y[i];
                }
            }
            for (int j = 0; j < vertex_number; j++) {
                pi[j] = pi[j] + q[j];
            }
            for (int j = 0; j < vertex_number; j++) {
                r[j] = r_[j];
                q[j] = 0;
                r_[j] = 0;
            }
        }
        for (int j = 0; j < vertex_number; j++) {
            q[j] = w[L] * r[j] / y[L];
        }

        for (int j = 0; j < vertex_number; j++) {
            pi[j] = pi[j] + q[j];
        }

    } else {
        construct_DS_static(a, b);
//        printf("constructed\n");
        naive_initialize(r, x);
        double su;
        double p;
        double p_a;
        DSBucket *dsb;
        MyVector<int> *intIndex;
//        printf("initialized\n");
        start = getCurrentTime();
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < vertex_number; j++) {
                if (r[j] == 0) {
                    continue;
                } else {
                    auto *u = (dynagp::Vertex *) vList[j];
                    d_u = u->getDegree();
                    if (d_u < 500) {
                        value = y[i + 1] * r[j] / y[i] / pow(d_u, b);
                        u_adjacentList = u->getAdjacentList();
                        for (int s = 0; s < d_u; s++) {
                            const int &neighborID = u_adjacentList[s];
                            auto *neighbor_v = (dynagp::Vertex *) vList[neighborID - 1];
                            r_[neighborID] = r_[neighborID] + value / pow(neighbor_v->getDegree(), a);
                        }
                    } else {
                        su = y[i + 1] * r[j] / eps / y[i] / pow(d_u, b);
                        dsb = u->getDS();
                        intIndex = u->getNoEmtInx();
                        int bn = intIndex->size();
                        for (int o = 0; o < bn; o++) {
                            int s = intIndex->get_list()[o];
                            if (s == 0) {
                                p = su;
                            } else {
                                //leave the budget for dynamic
                                p = su / pow(pow(2, s - 1), b);
                                // static
                                // p = su / pow(pow(2, s), b);
                            }
                            int bs = dsb->sizeByIndex(s);
                            if (p >= 1) {
                                for (int k = 0; k < bs; k++) {
                                    DSBucketElement *e = dsb->getElement(s, k);
                                    const int &neighborID = e->get_neighbor_id();
                                    auto *neighbor_v = (dynagp::Vertex *) vList[neighborID - 1];
                                    int d_v = neighbor_v->getDegree();
                                    p_a = su / pow(d_v, a);
                                    if (p_a >= 1) {
                                        r_[neighborID] = r_[neighborID] + p_a * eps;
                                    } else {
                                        double randomValue = generateRandomInt(0, 10000) / 10000.0;
                                        if (randomValue < p_a) {
                                            r_[neighborID] = r_[neighborID] + eps;
                                        }
                                    }
                                }
                            } else {
                                int k = 0;
                                while (true) {
                                    double randomValue = generateRandomInt(1, 10000) / 10000.0;
                                    int z = ceil(log2(randomValue) / log2(1 - p));
                                    k = k + z;
                                    // k < 0 is when z is too large to be held and 'become' negative
                                    if (k >= bs || k < 0) {
                                        break;
                                    }
                                    DSBucketElement *e = dsb->getElement(s, k);
                                    const int &neighborID = e->get_neighbor_id();
                                    auto *neighbor_v = (dynagp::Vertex *) vList[neighborID - 1];
                                    int d_v = neighbor_v->getDegree();
                                    p_a = su / pow(d_v, a);
                                    randomValue = generateRandomInt(0, 10000) / 10000.0;
                                    if (randomValue < p_a) {
                                        r_[neighborID] = r_[neighborID] + eps;
                                    }
                                }
                            }
                        }
                    }
                }
                q[j] = q[j] + w[i] * r[j] / y[i];
            }
            for (int j = 0; j < vertex_number; j++) {
                pi[j] = pi[j] + q[j];
            }
            for (int j = 0; j < vertex_number; j++) {
                r[j] = r_[j];
                q[j] = 0;
                r_[j] = 0;
            }

            // for debug: check the value changes
            // double sum = 0;
            // for (int ii = 0; ii < vertex_number; ii++) {
            //     sum += pi[ii];
            // }
            // printf("%.9lf \n", sum);
        }

        for (int j = 0; j < vertex_number; j++) {
            q[j] = w[L] * r[j] / y[L];
        }

        for (int j = 0; j < vertex_number; j++) {
            pi[j] = pi[j] + q[j];
        }

    }

    delete[] r;
    delete[] r_;
    delete[] q;
    delete[] y;

    // for debug: check the value changes
    // double sum = 0;
    // for (int ii = 0; ii < vertex_number; ii++) {
    //     sum += pi[ii];
    // }
    // printf("%.9lf \n", sum);

    double end = getCurrentTime();
    printf("query time: %.9lf\n", end - start);


    return pi;
}

void Graph::DAGP_initialize() {
    T.reserve(vertex_number);
    for (int i = 0; i < vertex_number; i++) {
        T.push_back(vList[i]->getDegree());
    }
}

int Graph::update_DS(int _vID1, int _vID2, float a, float b) {
    auto *v1 = (dynagp::Vertex *) vList[_vID1 - 1];
    auto *v2 = (dynagp::Vertex *) vList[_vID2 - 1];

    int d_v1 = v1->getDegree();
    int d_v2 = v2->getDegree();
    auto *buckets1 = v1->getDS();
    auto *buckets2 = v2->getDS();

    int ori_ind1 = v1->getElementIndex(_vID2).first;
    int ori_eleInd1 = v1->getElementIndex(_vID2).second;

    int ori_ind2 = v2->getElementIndex(_vID1).first;
    int ori_eleInd2 = v2->getElementIndex(_vID1).second;

    int l1 = buckets1->listSize();
    int l2 = buckets2->listSize();
    int ind1 = floor(log2(d_v1));
    int ind2 = floor(log2(d_v2));
    if (ind1 != ori_ind1) {
        auto *dse1 = new DSBucketElement(_vID2);

        //remove from the previous bucket
        int s1 = buckets1->sizeByIndex(ori_ind1);
        if (ori_eleInd1 == s1 - 1) {
            buckets1->DeleteElement(ori_ind1, ori_eleInd1);
            v1->neighborElementIndexMap.erase(_vID2);
        } else {
            int swappedNeighId = buckets1->getElement(ori_ind1, s1 - 1)->get_neighbor_id();
            buckets1->DeleteElement(ori_ind1, ori_eleInd1);
            v1->neighborElementIndexMap.erase(_vID2);
            v1->neighborElementIndexMap[swappedNeighId] = {ori_ind1, ori_eleInd1};
        }

        //insert to the new bucket
        buckets1->InsertNewELement(ind1, dse1);
        v1->neighborElementIndexMap[_vID2] = {ind1, dse1->get_element_index()};


    }
    if (ind2 != ori_ind2) {
        //remove from the previous bucket
        int s2 = buckets2->sizeByIndex(ori_ind2);
        if (ori_eleInd2 == s2 - 1) {
            buckets2->DeleteElement(ori_ind2, ori_eleInd2);
            v2->neighborElementIndexMap.erase(_vID1);
        } else {
            int swappedNeighId = buckets2->getElement(ori_ind2, s2 - 1)->get_neighbor_id();
            buckets2->DeleteElement(ori_ind2, ori_eleInd2);
            v2->neighborElementIndexMap.erase(_vID1);
            v2->neighborElementIndexMap[swappedNeighId] = {ori_ind2, ori_eleInd2};
        }

        //insert to the new bucket
        auto *dse2 = new DSBucketElement(_vID1);
        buckets2->InsertNewELement(ind2, dse2);
        v2->neighborElementIndexMap[_vID1] = {ind2, dse2->get_element_index()};

    }
    return 0;
}

int Graph::naive_insert(int _vID1, int _vID2, float a, float b) {
    auto *v1 = (dynagp::Vertex *) vList[_vID1 - 1];
    auto *v2 = (dynagp::Vertex *) vList[_vID2 - 1];
    // check if these vertices already exist. If not, create new vertex/vertices
    // with vertex factory and insert it/them into unordered map.
    if (v1 == NULL) {
        v1 = (dynagp::Vertex *) createVertex(_vID1);
        auto *buckets = new DSBucket();
        v1->assignDS(buckets);
    }
    if (v2 == NULL) {
        v2 = (dynagp::Vertex *) createVertex(_vID2);
        auto *buckets = new DSBucket();
        v2->assignDS(buckets);
    }
    v1->insertNeighbor(_vID2);
    v2->insertNeighbor(_vID1);

    int d_v1 = v1->getDegree();
    int d_v2 = v2->getDegree();
    auto *buckets1 = v1->getDS();
    auto *buckets2 = v2->getDS();

    int ind1 = floor(log2(d_v1));
    auto *dse1 = new DSBucketElement(_vID2);
    buckets1->InsertNewELement(ind1, dse1);
    v1->neighborElementIndexMap[_vID2] = {ind1, dse1->get_element_index()};

//    if (ind1 != floor(log2(d_v1 - 1))) {
//        for (int i = 0; i < d_v1; i++) {
//            int neighID = v1->getAdjacentList()[i];
//            if (neighID == _vID2) {
//                continue;
//            } else {
//                update_DS(_vID1, neighID, a, b);
//            }
//        }
//    }
//
    for (int i = 0; i < d_v1; i++) {
        int neighID = v1->getAdjacentList()[i];
        if (neighID == _vID2) {
            continue;
        } else {
            update_DS(_vID1, neighID, a, b);
        }
    }
    for (int i = 0; i < floor(log2(d_v1 - 1)); i++) {
        if (d_v1 == 1);
    }

    int ind2 = floor(log2(d_v2));
    auto *dse2 = new DSBucketElement(_vID1);

    buckets2->InsertNewELement(ind2, dse2);
    v2->neighborElementIndexMap[_vID1] = {ind2, dse2->get_element_index()};

//    if (ind2 != floor(log2(d_v2 - 1))) {
//
//        for (int i = 0; i < d_v2; i++) {
//            int neighID = v2->getAdjacentList()[i];
//            if (neighID == _vID1) {
//                continue;
//            } else {
//                update_DS(_vID2, neighID, a, b);
//            }
//        }
//    }

    for (int i = 0; i < d_v2; i++) {
        int neighID = v2->getAdjacentList()[i];
        if (neighID == _vID1) {
            continue;
        } else {
            update_DS(_vID2, neighID, a, b);
        }
    }

    for (int i = 0; i < floor(log2(d_v2 - 1)); i++) {
        if (d_v2 == 1);
    }
    return 0;
}

int Graph::naive_delete(int _vID1, int _vID2, float a, float b) {
    auto *v1 = (dynagp::Vertex *) vList[_vID1 - 1];
    auto *v2 = (dynagp::Vertex *) vList[_vID2 - 1];

    if (v1 == NULL || v2 == NULL)
        return 1;
    v1->deleteNeighbor(_vID2);
    v2->deleteNeighbor(_vID1);
    auto *buckets1 = v1->getDS();
    auto *buckets2 = v2->getDS();
    int ind1 = v1->getElementIndex(_vID2).first;
    int ind2 = v2->getElementIndex(_vID1).first;
    int eleInd1 = v1->getElementIndex(_vID2).second;
    int eleInd2 = v2->getElementIndex(_vID1).second;
    int l1 = buckets1->sizeByIndex(ind1);
    int l2 = buckets2->sizeByIndex(ind2);

    if (eleInd1 == l1 - 1) {
        buckets1->DeleteElement(ind1, eleInd1);
        v1->neighborElementIndexMap.erase(_vID2);
    } else {
        int swappedNeighId = buckets1->getElement(ind1, l1 - 1)->get_neighbor_id();
        buckets1->DeleteElement(ind1, eleInd1);
        v1->neighborElementIndexMap.erase(_vID2);
        v1->neighborElementIndexMap[swappedNeighId] = {ind1, eleInd1};
    }

    if (eleInd2 == l2 - 1) {
        buckets2->DeleteElement(ind2, eleInd2);
        v2->neighborElementIndexMap.erase(_vID1);
    } else {
        int swappedNeighId = buckets2->getElement(ind2, l2 - 1)->get_neighbor_id();
        buckets2->DeleteElement(ind2, eleInd2);
        v2->neighborElementIndexMap.erase(_vID1);
        v2->neighborElementIndexMap[swappedNeighId] = {ind2, eleInd2};
    }
    int d_v1 = v1->getDegree();
    int d_v2 = v2->getDegree();

//    if (ind1 != floor(log2(d_v1))) {
//        for (int i = 0; i < d_v1; i++) {
//            int neighID = v1->getAdjacentList()[i];
//            if (neighID == _vID2) {
//                continue;
//            } else {
//                update_DS(_vID1, neighID, a, b);
//            }
//        }
//    }

    for (int i = 0; i < d_v1; i++) {
        int neighID = v1->getAdjacentList()[i];
        if (neighID == _vID2) {
            continue;
        } else {
            update_DS(_vID1, neighID, a, b);
        }
    }

//    if (ind2 != floor(log2(d_v2))) {
//        for (int i = 0; i < d_v2; i++) {
//            int neighID = v2->getAdjacentList()[i];
//            if (neighID == _vID1) {
//                continue;
//            } else {
//                update_DS(_vID2, neighID, a, b);
//            }
//        }
//    }

    for (int i = 0; i < d_v2; i++) {
        int neighID = v2->getAdjacentList()[i];
        if (neighID == _vID1) {
            continue;
        } else {
            update_DS(_vID2, neighID, a, b);
        }
    }
    return 0;
}

int Graph::DAGP_insert(int _vID1, int _vID2, float a, float b) {
    auto *v1 = (dynagp::Vertex *) vList[_vID1 - 1];
    auto *v2 = (dynagp::Vertex *) vList[_vID2 - 1];
    // check if these vertices already exist. If not, create new vertex/vertices
    // with vertex factory and insert it/them into unordered map.
    if (v1 == NULL) {
        v1 = (dynagp::Vertex *) createVertex(_vID1);
        auto *buckets = new DSBucket();
        buckets->extendListSize(1);
        v1->assignDS(buckets);
        T.push_back(1);
    }
    if (v2 == NULL) {
        v2 = (dynagp::Vertex *) createVertex(_vID2);
        auto *buckets = new DSBucket();
        buckets->extendListSize(1);
        v2->assignDS(buckets);
        T.push_back(1);
    }

    v1->insertNeighbor(_vID2);
    v2->insertNeighbor(_vID1);

    int d_v1 = v1->getDegree();
    int d_v2 = v2->getDegree();
    auto *buckets1 = v1->getDS();
    auto *buckets2 = v2->getDS();

    int l1 = buckets1->listSize();
    int l2 = buckets2->listSize();
    int ind1 = floor(log2(d_v1));
    auto *dse1 = new DSBucketElement(_vID2);
    buckets1->InsertNewELement(ind1, dse1);
    v1->neighborElementIndexMap[_vID2] = {ind1, dse1->get_element_index()};


    int ind2 = floor(log2(d_v2));
    auto *dse2 = new DSBucketElement(_vID1);

    buckets2->InsertNewELement(ind2, dse2);
    v2->neighborElementIndexMap[_vID1] = {ind2, dse2->get_element_index()};

    if (d_v1 > pow(2, a + b) * T[_vID1 - 1]) {
        for (int i = 0; i < d_v1; i++) {
            int neighID = v1->getAdjacentList()[i];
            if (neighID == _vID2) {
                continue;
            } else {
                update_DS(_vID1, neighID, a, b);
            }
        }
        T[_vID1 - 1] = d_v1;
    }

    if (d_v2 > pow(2, a + b) * T[_vID2 - 1]) {
        for (int i = 0; i < d_v2; i++) {
            int neighID = v2->getAdjacentList()[i];
            if (neighID == _vID1) {
                continue;
            } else {
                update_DS(_vID2, neighID, a, b);
            }
        }
        T[_vID2 - 1] = d_v2;
    }
    return 0;
}

int Graph::DAGP_delete(int _vID1, int _vID2, float a, float b) {
    auto *v1 = (dynagp::Vertex *) vList[_vID1 - 1];
    auto *v2 = (dynagp::Vertex *) vList[_vID2 - 1];

    if (v1 == NULL || v2 == NULL)
        return 1;
    v1->deleteNeighbor(_vID2);
    v2->deleteNeighbor(_vID1);
    auto *buckets1 = v1->getDS();
    auto *buckets2 = v2->getDS();
    int ind1 = v1->getElementIndex(_vID2).first;
    int ind2 = v2->getElementIndex(_vID1).first;
    int eleInd1 = v1->getElementIndex(_vID2).second;
    int eleInd2 = v2->getElementIndex(_vID1).second;
    int l1 = buckets1->sizeByIndex(ind1);
    int l2 = buckets2->sizeByIndex(ind2);

    if (eleInd1 == l1 - 1) {
        buckets1->DeleteElement(ind1, eleInd1);
        v1->neighborElementIndexMap.erase(_vID2);
    } else {
        int swappedNeighId = buckets1->getElement(ind1, l1 - 1)->get_neighbor_id();
        buckets1->DeleteElement(ind1, eleInd1);
        v1->neighborElementIndexMap.erase(_vID2);
        v1->neighborElementIndexMap[swappedNeighId] = {ind1, eleInd1};
    }

    if (eleInd2 == l2 - 1) {
        buckets2->DeleteElement(ind2, eleInd2);
        v2->neighborElementIndexMap.erase(_vID1);
    } else {
        int swappedNeighId = buckets2->getElement(ind2, l2 - 1)->get_neighbor_id();
        buckets2->DeleteElement(ind2, eleInd2);
        v2->neighborElementIndexMap.erase(_vID1);
        v2->neighborElementIndexMap[swappedNeighId] = {ind2, eleInd2};
    }
    int d_v1 = v1->getDegree();
    int d_v2 = v2->getDegree();


    if (d_v1 > pow(2, a + b) * T[_vID1 - 1]) {
        for (int i = 0; i < d_v1; i++) {
            int neighID = v1->getAdjacentList()[i];
            if (neighID == _vID2) {
                continue;
            } else {
                update_DS(_vID1, neighID, a, b);
            }
        }
        T[_vID1 - 1] = d_v1;
    }
    if (d_v2 > pow(2, a + b) * T[_vID2 - 1]) {
        for (int i = 0; i < d_v2; i++) {
            int neighID = v2->getAdjacentList()[i];
            if (neighID == _vID1) {
                continue;
            } else {
                update_DS(_vID2, neighID, a, b);
            }
        }
        T[_vID2 - 1] = d_v2;
    }

    return 0;
}

