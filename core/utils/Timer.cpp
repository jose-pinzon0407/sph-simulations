#include "Timer.h"

Timer::Timer(const char* n):name{n},total{0.0}{

}

void Timer::Start(){
    this->start = std::chrono::high_resolution_clock::now();
}

void Timer::Stop(){
    this->end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time = end -start;
    this->total += elapsed_time;
}

double Timer::GetTotalTimeInSeconds(){
    timelapse time = std::chrono::duration_cast<timelapse>(total);
    return time.count()/1000000.0;
}

std::string& Timer::GetName(){
    return this->name;
}