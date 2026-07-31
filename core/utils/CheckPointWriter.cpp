#include "CheckPointWriter.h"
// #include "ensemble/Ensemble.h"

CheckPointWriter::CheckPointWriter(const char* pref, int freq):prefix{pref},write_frequency{freq}{
}


CheckPointWriter::CheckPointWriter(std::string pref, int freq):prefix{pref},write_frequency{freq}{

}

void CheckPointWriter::WriteCheckPoint(int step, Container& particles){

    if(step%write_frequency != 0) return;

    std::string file_name = prefix + std::to_string(step) + ".txt";
    this->CreateNewFile(file_name);

    int N = particles.Size();
    for(int i=0;i<N;++i){
        Particle particle = particles[i];
        file_stream<<particle.r[0]<<" "<<particle.r[1]<<" "<<particle.r[2]<<std::endl;
    }
    

    this->ReleaseStream();

}