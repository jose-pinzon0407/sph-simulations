#pragma once

#include<fstream>
#include<vector>
#include<string>
#include<iostream>
#include<map>
#include<regex>

#include "ConfigVariables.hpp"

/**
 * Class for parsing the simulation input parameters
 */

class InputParser{
    private:
    std::map<std::string,std::string> parameters;
    std::vector<std::string> input;
    std::string names;

    public:

    std::vector<std::string>& GetInput();
    void ReadInputFile(const char* filename);
    void ReadInputFile(const std::string& filename);

    std::string GetParameter(const char* name) const;
    std::string GetParameter(std::string& name) const;

    int GetIntParameter(const char* name) const;
    Real GetRealParameter(const char* name) const;
    private:
    /**
     * Find the numerical value associated to the given parameter name
     */
    std::string GetValue(std::string& input);
    std::string TrimLine(const std::string& line);
};