#include "DS.h"

MyVector<DSBucketElement*> DSBucket::EmptyBucket;

void DSBucket::DeleteElement(int bucketIndex, int elementIndex) {
    MyVector<DSBucketElement*>& bList = buckList[bucketIndex];
    if (bList.size() == 1) {
        bList.pop_back();
        //this->shrinkToFit();
        return;
    }
    int old_size = bList.size();
    DSBucketElement* swap_e = bList[old_size - 1];
    bList[elementIndex] = swap_e;
    swap_e->element_index = elementIndex;
    bList.pop_back();
    if (bList.size() < bList.capacity() * 0.5) {
        bList.shrink_to_fit();
    }
}

DSBucket::~DSBucket() {
    for (int i = 0; i < buckList.size(); ++i) {
        buckList[i].release_space();
    }
    buckList.release_space();

}

int DSBucket::listSize() {
    return buckList.size();
}

void DSBucket::extendListSize(int index) {
    while (buckList.size() < index) {
        buckList.push_back(EmptyBucket);
    }
}

int DSBucket::InsertNewELement(int i, DSBucketElement* e) {
    this->extendListSize(i + 1);
    buckList[i].push_back(e);
    e->element_index = buckList[i].size() - 1;
    return buckList[i].size();
}

bool DSBucket::CheckEmptyByIndex(int index) {
    if (buckList.size() < index) {
        return true;
    } else if (buckList[index].size() == 0) {
        return true;
    }
    return false;

}

int DSBucket::sizeByIndex(int i) {
    if (buckList.size() <= i) {
        return 0;
    }
    return buckList[i].size();
}

void DSBucket::shrinkToFit() {
    int i = buckList.size() - 1;
    while (buckList[i].size() == 0) {
        buckList.pop_back();
        i--;
    }
}

void DSBucket::ReleaseSpace() {
    for (int i = 0; i < buckList.size(); ++i) {
        buckList[i].release_space();
    }
    buckList.release_space();
}

DSBucketElement* DSBucket::getElement(int i, int j) const {
    return buckList[i][j];
}




