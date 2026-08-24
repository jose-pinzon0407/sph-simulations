#include "TimerManager.h"

void TimerManager::OutputTimersToFile(int N, int steps){
    this->OpenInputAppend(file_name.c_str());
    std::ofstream& output = this->GetStream();

    int num_threads=1;
    #ifdef OPENMP
    #pragma omp parallel 
    {
        #pragma omp critical
        num_threads=omp_get_num_threads();
    }
    #endif
    output<<std::fixed<<std::setw(20)<<std::left<<num_threads<<"\t";
    output<<std::fixed<<std::setw(20)<<std::left<<steps<<"\t";
    output<<std::fixed<<std::setw(20)<<std::left<<N<<"\t";
    
    for(int i=0;i<timers.size();i++){
        output<<std::fixed<<std::setw(20)<<std::setprecision(8)<<std::left<<timers[i].GetTotalTimeInSeconds()<<"\t";
    }
    
    output<<std::endl;

    this->ReleaseStream();
}

void TimerManager::PrintTimers(){
    double time_sum=0.0;
    std::cout<<"[HighPerformanceClock] Measured times: \n";
    for(int i=0;i<timers.size();++i){
        std::cout<<timers[i].GetName()<<"="<<timers[i].GetTotalTimeInSeconds()<<"\n";
        time_sum += timers[i].GetTotalTimeInSeconds();
    }

    std::cout<<"Total time="<<time_sum-GetTimer("main_loop")->GetTotalTimeInSeconds()<<"\n";
}

void TimerManager::Init(const char* path){
    if(path){
        file_name=path;
    }

    bool file_was_created=this->CreateIfNotExists(file_name.c_str());

    if(file_was_created){
        std::cout<<"[TimerManager] Set header"<<std::endl;
        this->SetFileHeader();
    }
}

void TimerManager::SetFileHeader(){

    file_stream<<std::fixed<<std::setw(20)<<std::left<<"#threads"<<"\t";
    file_stream<<std::fixed<<std::setw(20)<<std::left<<"#steps"<<"\t";
    file_stream<<std::fixed<<std::setw(20)<<std::left<<"N"<<"\t";
    for(int i=0;i<timers.size();++i){
        file_stream<<std::fixed<<std::setw(20)<<std::left<<timers[i].GetName()<<"\t";
    }
    file_stream<<std::endl;
    this->ReleaseStream();
}

void TimerManager::AddTimer(const char* name){
    timers.emplace_back(Timer{name});
}

void TimerManager::Start(const char* name){
    auto timer = GetTimer(name);
    
    if(timer){
        timer->Start();
    }
}

void TimerManager::Stop(const char* name){
    auto timer = GetTimer(name);
    
    if(timer){
        timer->Stop();
    }
}

Timer* TimerManager::GetTimer(const char* name){

    auto timer = std::find_if(timers.begin(),timers.end(),
                              [&name](Timer& t){
                                return t.GetName()==name;
                              }  
                             );

    if(timer!=timers.end()){
        return &(*timer);
    }

    return nullptr;

}