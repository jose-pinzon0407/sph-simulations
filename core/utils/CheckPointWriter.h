#include "Printer.h"
#include"container/Container.h"


class CheckPointWriter:public FileSystem{
    protected:
    std::string prefix;
    int write_frequency;

    public:
    CheckPointWriter(const char* prefix, int freq=1);
    CheckPointWriter(std::string prefix, int freq=1);
    void WriteCheckPoint(int step, Container& container);

};