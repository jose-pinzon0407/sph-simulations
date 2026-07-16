#pragma once

#include "container/Container.h"
#include "container/SoAContainer.h"
// #include "alg/NeighborSearch.h"

class SortedSoAContainer:public SoAContainer{


    public:

    void Build();
    void Rebuild();
    void SortIntoCells();
    
    protected:
    

    // NeighborSearch neighbor_search;
};