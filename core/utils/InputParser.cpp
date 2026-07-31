#include "InputParser.h"
#include<regex>


void InputParser::ReadInputFile(const char* filename){
    std::string name{filename};
    ReadInputFile(name);
}

void InputParser::ReadInputFile(const std::string& filename){
    std::cout<<"Configuration folder: "<<CONFIG_DIR<<std::endl;
    std::string line;
    std::ifstream file{CONFIG_DIR  +filename};
    // file.open(filename);

    if(!file.is_open()){
        std::cerr<<"[Error]Cannot find configuration file"<<std::endl;
        return;
    }

    while(std::getline(file,line)){
        size_t pos_equal = line.find("=");

        if(pos_equal != std::string::npos){
            std::string parameter = TrimLine(line.substr(0,pos_equal));
            std::string value = TrimLine(line.substr(pos_equal+1));

            //add to map
            parameters.insert({parameter,value});
        }
    }

}

std::vector<std::string>& InputParser::GetInput(){
    return input;
}

std::string InputParser::GetParameter(const char* name) const{
    if(parameters.find(name)==parameters.end()){
        std::string dummy_value{"0"};
        return dummy_value;
    }
    else{
       return parameters.at(name); 
    }
    
    // return parameters[name];

}

Real InputParser::GetRealParameter(const char* name) const{
    if(parameters.find(name)==parameters.end()){
        std::string dummy_value{"0"};
        return std::stod(dummy_value);
    }
    else{
       return std::stod(parameters.at(name)); 
    }
    
    // return parameters[name];

}

int InputParser::GetIntParameter(const char* name) const{
    if(parameters.find(name)==parameters.end()){
        std::string dummy_value{"0"};
        return std::stoi(dummy_value);
    }
    else{
       return std::stoi(parameters.at(name)); 
    }
    
    // return parameters[name];

}

std::string InputParser::GetParameter(std::string& name) const{
    return this->GetParameter(name.c_str());
}

std::string InputParser::GetValue(std::string& input){
    std::string output = std::regex_replace(input,
                                            std::regex(R"(\w+\s*=\s*([-+]?\d*\.?\d+([eE][-+]?\d+)?))"),
                                            "$1"
                                            );
    return output;

}

std::string InputParser::TrimLine(const std::string& line){
    size_t start = line.find_first_not_of(" \t");
    size_t end = line.find_last_not_of(" \t");
    return (start == std::string::npos) ?"":line.substr(start, end - start +1);
}