#pragma once
#include<cmath>
#include<iostream>
#include<array>
#include<algorithm>

template<typename T>
std::ostream& operator<<(std::ostream& os, std::array<T,3> const& v ){
    // os<<"(";
    for(auto it=v.begin();it!=v.end();++it){
        os<<*it<<" ";
    }
    // os<<")";
    return os;
}


template<typename T>
 std::array<T,3>& operator+=(std::array<T,3>& v1, std::array<T,3> const& v2){
    
    if(v1.size() != v2.size()){
        std::exit(EXIT_FAILURE);
    }

    // if(v1 == v2){
    //     return v1;
    // }
    
    for(size_t i=0;i<v1.size();++i){
        v1[i] += v2[i];
    }
    return v1;
}

template<typename T>
 std::array<T,3>& operator +=(std::array<T,3>& V, const T cst){
    for(size_t i=0;i<V.size();++i){
        V[i] += cst;
    }
    return V;
}

template<typename T>
 std::array<T,3> operator+(std::array<T,3> const& v1, std::array<T,3> const& v2){


    if(v1.size() != v2.size()){
        std::exit(EXIT_FAILURE);
    }
    
    std::array<T,3> v3=v1;

    for(int i=0;i<v3.size();++i){
        v3[i] += v2[i];
    }

    // v3 += v2;
    return v3;
}

template<typename T>
 std::array<T,3>& operator-=(std::array<T,3>& v1, std::array<T,3> const& v2){
    
    if(v1.size() != v2.size()){
        std::exit(EXIT_FAILURE);
    }

    // if(v1 == v2){
    //     return v1;
    // }
    
    for(size_t i=0;i<v1.size();++i){
        v1[i] -= v2[i];
    }
    return v1;
}

template<typename T>
 std::array<T,3> operator-(std::array<T,3> const& v1, std::array<T,3> const& v2){
    std::array<T,3> v3=v1;

    for(int i=0;i<v3.size();++i){
        v3[i] -= v2[i];
    }

    return v3;
}

template<typename T, typename U>
 std::array<T,3>& operator*=(std::array<T,3>& v1, const U val){
    for(size_t i=0;i<v1.size();++i){
        v1[i] = (static_cast<T>(val)*v1[i]);
    }
    return v1;
}

template<typename T, typename U>
 std::array<T,3> operator*(std::array<T,3> const& v1, const U val){
    std::array<T,3> scaled_vec = v1;
    for(size_t i=0;i<v1.size();++i){
        scaled_vec[i] = (static_cast<T>(val)*scaled_vec[i]);
    }
    return scaled_vec;
}

template<typename T, typename U>
 std::array<T,3> operator*(const U val, std::array<T,3> const& v1){
    std::array<T,3> scaled_vec = v1;
    for(size_t i=0;i<v1.size();++i){
        scaled_vec[i] = (static_cast<T>(val)*scaled_vec[i]);
    }
    return scaled_vec;
}

template<typename T, typename U>
 std::array<T,3>& operator/=(std::array<T,3>& v1, U val){
    if(val == 0.0){
        std::exit(EXIT_FAILURE);
    }
    for(size_t i=0;i<v1.size();++i){
        v1[i] = (v1[i]/static_cast<T>(val));
    }
    return v1;
}

template<typename T>
double SquaredNorm(std::array<T,3> const& V){
    double n2 = 0.0;    

    for(size_t i=0;i<V.size();++i){
        n2 += V[i]*V[i];
    }

    return n2;
}

template<typename T>
double VecNorm(std::array<T,3> const& V){
    return std::sqrt(SquaredNorm(V));
}


template<typename T>
 double DotProduct(std::array<T,3> const& v1, std::array<T,3> const& v2){
    
    double dot=0;
    
    for(size_t i=0;i<v1.size();++i){
        dot += v1[i]*v2[i];
    }

    return dot;
}

template<typename T>
 std::array<T,3>& VectorNormalize(std::array<T,3>& V){
    double norm = VecNorm(V);
    V /= norm;

    return V;
}

template<typename T>
 std::array<T,3> VectorMultiplyEntries(std::array<T,3> const& v1,std::array<T,3> const& v2){
    std::array<T,3> output;
    for(int i=0;i<output.size();++i){
        output[i] = v1[i]*v2[i];
    }

    return output;
}

template<typename T>
Mat3d OuterProduct(std::array<T,3> const& v1,std::array<T,3> const& v2){
    Mat3d mat;

    for(int i=0;i<3;++i){
        for(int j=0;j<3;++j){
            mat[i][j] = v1[i]*v2[j];
        }
    }

    return mat;

}