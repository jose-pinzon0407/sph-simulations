#pragma once

#include<string>
#include<chrono>
/**
 * \brief Measure time by wrapping around a routine(s)
 * 
 * 
 * 
 */

class Timer{
    using timelapse = std::chrono::microseconds;
    using clock_type = std::chrono::high_resolution_clock;
    private:
    std::string name;
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    std::chrono::duration<double> total;
    public:

    Timer(const char* name);
    void Start();
    void Stop();
    double GetTotalTimeInSeconds();
    std::string& GetName();
    

};