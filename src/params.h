//! R-exposed builder class for LSHConstructionParameters

#ifndef FALCONNR_PARAMS_H
#define FALCONNR_PARAMS_H

#include <string>
#include <map>

#include <Rcpp.h>
#include "falconn/lsh_nn_table.h"

class LshParameterSetter {
  public:
    typedef std::map<std::string, falconn::DistanceFunction> distancesMap;
    typedef std::map<std::string, falconn::StorageHashTable> storageTypesMap;
    typedef std::map<std::string, falconn::LSHFamily>        familiesMap;

    static const familiesMap     families;
    static const distancesMap    distances;
    static const storageTypesMap storageTypes;

    LshParameterSetter(int n, int d);

    falconn::LSHConstructionParameters params() const;

    LshParameterSetter& withDefaults(std::string distance  = "euclidean_squared");
    LshParameterSetter& distance(std::string distance);
    LshParameterSetter& numHashFunctions(int funcs);
    LshParameterSetter& numHashTables(int tables);
    LshParameterSetter& storage(std::string storage);
    LshParameterSetter& family(std::string family);
    LshParameterSetter& rotations(int numRotations);
    Rcpp::List asList();

  private:

    int                                _n;
    int                                _d;
    falconn::LSHConstructionParameters _p;
};

#endif

// Local Variables:
// mode: c++
// End:
