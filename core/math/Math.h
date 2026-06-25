#pragma once
#include<math.h>
#include<iostream>
#include<array>
#include<algorithm>
#include<vector>
#include "ConfigVariables.hpp"

using Vec3d = std::array<Real,3>;
using Vec3i = std::array<int,3>;
using Mat3d = std::array<std::array<Real,3>,3>;
using ScalarField = std::vector<Real>;
using VectorField = std::vector<Vec3d>;



template<typename T>
std::ostream& operator<<(std::ostream& os, std::array<T,3> const& v );

template<typename T>
inline std::array<T,3>& operator+=(std::array<T,3>& v1, std::array<T,3> const& v2);

template<typename T>
inline std::array<T,3>& operator +=(std::array<T,3>& V, const T cst);

template<typename T>
inline std::array<T,3> operator+(std::array<T,3> const& v1, std::array<T,3> const& v2);

template<typename T>
inline std::array<T,3>& operator-=(std::array<T,3>& v1, std::array<T,3> const& v2);

template<typename T>
inline std::array<T,3> operator-(std::array<T,3> const& v1, std::array<T,3> const& v2);

template<typename T, typename U>
inline std::array<T,3>& operator*=(std::array<T,3>& v1, const U val);

template<typename T, typename U>
inline std::array<T,3> operator*(std::array<T,3> const& v1, const U val);

template<typename T, typename U>
inline std::array<T,3> operator*(U val, std::array<T,3> const& v1);

template<typename T, typename U>
inline std::array<T,3>& operator/=(std::array<T,3>& v1, U val);

template<typename T>
inline double SquaredNorm(std::array<T,3> const& V);

template<typename T>
inline double DotProduct(std::array<T,3> const& v1, std::array<T,3> const& v2);

template<typename T>
inline std::array<T,3>& VectorNormalize(std::array<T,3>& V);

template<typename T>
inline std::array<T,3> VectorMultiplyEntries(std::array<T,3> const& v1,std::array<T,3> const& v2);

template<typename T>
double VecNorm(std::array<T,3> const& V);

/**
 * \brief Tensor product between two vectors
 */
template<typename T>
std::array<std::array<double,3>,3> OuterProduct(std::array<T,3> const& v1,std::array<T,3> const& v2);

#include "Math_impl.h"