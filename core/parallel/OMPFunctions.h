#pragma once

#include "OpenMPSetup.h"

#ifdef OPENMP
void OMPConfiguration(){
    #pragma omp parallel
    {
        #pragma omp master
        {
        std::cout<<"[OpenMP] Parallelization using "<<omp_get_num_threads()<<" threads"<<std::endl;
        int affinity = omp_get_proc_bind();
        std::string affinity_name;
        if(affinity==0)
        affinity_name = "false";
        if(affinity==1)
        affinity_name = "true";
        if(affinity==2)
        affinity_name = "master";
        if(affinity==3)
        affinity_name = "close";
        if(affinity==4)
        affinity_name = "spread";
        std::cout<<"[OpenMP] Thread affinity policy '"<<affinity_name<<"' used"<<std::endl;
        }
    }
    #pragma omp parallel 
    {
        int ID= omp_get_thread_num();
        #pragma omp critical
        std::cout<<"Thread "<<ID<<" from "<<omp_get_num_threads()<<" running on CPU: "<<sched_getcpu()<<std::endl;
    }

    #pragma omp parallel
    {
        #pragma omp master
        std::cout<<"[OpenMP] Cores: "<<omp_get_num_threads()<<std::endl;
    }
    // omp_sched_t scheduling_type;
    //     int chunk;
    //         omp_get_schedule(&scheduling_type,&chunk);
        
        omp_sched_t scheduling_type;
        int chunk;
        omp_get_schedule(&scheduling_type,&chunk);
        omp_sched_t base = (omp_sched_t)(scheduling_type & 0xF);
        std::cout<<"[OpenMP] Scheduling type: ";
        switch(base){
            case omp_sched_static:
                std::cout<<"static\n";
                break;
            case omp_sched_dynamic:
                std::cout<<"dynamic\n";
                break;
            case omp_sched_guided:
                std::cout<<"guided\n";
                break;
            case omp_sched_auto:
                std::cout<<"auto\n";
                break;
            default:
			    printf("other (implementation specific)\n");
			    break;
        }
    }
#endif