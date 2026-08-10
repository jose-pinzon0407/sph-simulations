
template<class NeighborSearchType>
GenericNeighborLists<NeighborSearchType>::GenericNeighborLists(double r):h2{r*r}{



}

template<class NeighborSearchType>
void GenericNeighborLists<NeighborSearchType>::AssembleNeighborList(SoAContainer& container, GenericNeighborSearch<NeighborSearchType>& neighbor_search){

    ResetNeighborList();

    size_t N = container.Size();
    Resize(N);

    //prepare codifier
    // encoder.ResetCodifier();
    // encoder.Init(N);

    #pragma omp parallel for schedule(runtime)
    for(size_t i=0;i<N;++i){
        uint cell_i = container.GetParticleCellIndex(i);

        Vec3i m_idx_i = neighbor_search.Unflatten(cell_i);

        for(int c1=-1;c1<=1;++c1){
            for(int c2=-1;c2<=1;++c2){
                
                Vec3i m_idx_j = m_idx_i;
                
                m_idx_j[0] += c1;
                
                if(m_idx_j[0] < 0 || m_idx_j[0] >= neighbor_search.GetCellsInDirectionX()) continue;
                
                m_idx_j[1] += c2;
                
                if(m_idx_j[1] < 0|| m_idx_j[1] >= neighbor_search.GetCellsInDirectionY()) continue;
                

                int cell_j = neighbor_search.Flatten(m_idx_j);


                int j_start = neighbor_search.GetCellStart(cell_j);
                int j_end = neighbor_search.GetCellEnd(cell_j);

                if(j_start < 0 ||  j_end<0){
                    continue;
                }
                
                for(int pj=j_start;pj<j_end;++pj){

                    Real dx = container.GetParticles().rx[i]-container.GetParticles().rx[pj];
                    Real dy = container.GetParticles().ry[i]-container.GetParticles().ry[pj];
                    Real dz = container.GetParticles().rz[i]-container.GetParticles().rz[pj];

                    Real r2 = dx*dx+dy*dy+dz*dz; 

                    if(r2>h2) continue;

                    //add neighbor
                    neighbor_lists[i].emplace_back(pj);



                }
            }
        }

        // std::sort(neighbor_lists[i].begin(),neighbor_lists[i].end());

        //Encode into differences
        // encoder.Encode(i,neighbor_lists[i]);

    }
}

template<class NeighborSearchType>
neighbors& GenericNeighborLists<NeighborSearchType>::GetNeighborList(){
    return neighbor_lists;
}

template<class NeighborSearchType>
void GenericNeighborLists<NeighborSearchType>::ResetNeighborList(){
    size_t N = neighbor_lists.size();
    #pragma omp parallel for schedule(runtime)
    for(int i=0;i<N;++i){
        neighbor_lists[i].clear();
    }

}

template<class NeighborSearchType>
void GenericNeighborLists<NeighborSearchType>::Resize(int N){
    neighbor_lists.resize(N);
}