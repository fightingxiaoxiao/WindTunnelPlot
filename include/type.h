#ifndef TYPE_H
#define TYPE_H

#define EIGEN_VECTORIZE
#define USE_INTEL_MKL_ALL
#include <Eigen/Dense>
#include <map>
#include <string>
#include <vector>

typedef Eigen::Vector3d Vector;

typedef double      Scalar;
typedef long        Label;
typedef std::string Word;

template <class T> using List = std::vector<T>;

template <class KEY, class VALUE> using Dict = std::map<KEY, VALUE>;

#endif