#pragma once

#include<memory>
#include<vector>
#include<iterator>

#include<iostream>
#include<iomanip>
#include<fstream>
#include<algorithm>
#include<filesystem>
 

class FileSystem{
    protected:

    std::ofstream file_stream;

    public:

    std::ofstream& GetStream();
    /**
     * \brief Checks if the file exists, clears if needed. 
     * 
     *  Uses the in flag to see if the file exists, then clears and creates if needed.
     */
    bool CreateNewFile(std::string name);
    bool CreateNewFile(const char* name);
    bool CreateIfNotExists(const char* name);
       /**
     * \brief Clears file and opens for writing
     * 
     * Uses the trunc and out flags
     */
    bool ExistsAndClear(const char* name);
    bool Exists(const char* path);

    bool OpenInputAppend(const char* name);

    void ReleaseStream();
};