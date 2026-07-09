#include "Printer.h"

std::ofstream& FileSystem::GetStream(){
    return this->file_stream;
}

bool FileSystem::CreateIfNotExists(const char* name){
    
    std::ifstream checking{name};
    if(checking.good()){
        return false;
    }
    
    CreateNewFile(name);
    return true;

}

bool FileSystem::Exists(const char* name){
    std::ifstream checking{name};
    if(checking.good()){
        return true;
    }else{
        return false;
    }

}

bool FileSystem::CreateNewFile(std::string name){
    file_stream.open(name,std::ios::in);
    if(file_stream.good()){
        file_stream.close();
        //std::cout<<"The file already exists, rewriting\n";
        file_stream.open(name,std::ios::trunc|std::ios::out);
        return false;
    }
    file_stream.open(name, std::ios::out);
    return true;
}

bool FileSystem::CreateNewFile(const char* name){
    std::string file_name{name};
    return CreateNewFile(file_name);
}

bool FileSystem::ExistsAndClear(const char* name){
    file_stream.open(name,std::ios::trunc|std::ios::out);

    if(file_stream.is_open()){
        file_stream.close();//release for next file
        return true;
    }
    std::cout<<"Could not open and clear file\n";
    return false;
}

bool FileSystem::OpenInputAppend(const char* name){
    file_stream.open(name, std::ios::app|std::ios::out);
    if(file_stream.is_open()){
        return true;
    }
    std::cout<<"Could not open file\n";
    return false;
}

void FileSystem::ReleaseStream(){
    file_stream.close();
}