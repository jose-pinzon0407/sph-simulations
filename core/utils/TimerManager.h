#pragma once

#include<vector>
#include<algorithm>

#include "Timer.h"
#include "Printer.h"
#include "parallel/OpenMPSetup.h"
/**
 * Class that handles the Timer class
 */

 class TimerManager:public FileSystem{

    public:
    void Init(const char* path=nullptr);
    void AddTimer(const char* name);
    Timer* GetTimer(const char* name);
    void OutputTimersToFile(int N, int steps);
    void PrintTimers();

    void Start(const char* name);
    void Stop(const char* name);

    private:

    void SetFileHeader();

    private:
    std::string file_name = "timer_manager.data";
    std::vector<Timer> timers;

 };