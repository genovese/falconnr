#include <Rcpp.h>
#include "falconn/lsh_nn_table.h"
#include <string>
#include <map>

using namespace Rcpp;

using falconn::construct_table;
using falconn::compute_number_of_hash_functions;
using falconn::DenseVector;
using falconn::DistanceFunction;
using falconn::LSHConstructionParameters;
using falconn::LSHFamily;
using falconn::LSHNearestNeighborTable;
using falconn::StorageHashTable;
using falconn::get_default_parameters;

typedef DenseVector<double> Point;

typedef std::map<std::string, LSHFamily>        familiesMap;
typedef std::map<std::string, DistanceFunction> distancesMap;
typedef std::map<std::string, StorageHashTable> storageTypesMap;

template <typename V>
V get( const std::map<std::string,V> &m, const std::string &key, const V &def) {
    typename std::map<std::string,V>::const_iterator it = m.find( key );
    if ( it == m.end() ) {
        return def;
    } else {
        return it->second;
    }
}

class LshParameterSetter {
  public:
    static const familiesMap families;
    static const distancesMap distances;
    static const storageTypesMap storageTypes;

    LshParameterSetter(int n, int d) : _n(n), _d(d) {
        setDefaults();
    }

    LSHConstructionParameters params() { return _p; }

    LshParameterSetter& setDefaults(std::string distance="euclidean_squared") {
        _p =
            get_default_parameters<Point>(_n, _d,
                                          get<DistanceFunction>(distances,
                                                                distance,
                                                                DistanceFunction::Unknown),
                                          true);
        return *this;
    }

    LshParameterSetter&  setDistance(std::string distance) {
        _p.distance_function = get<DistanceFunction>(distances,
                                                     distance,
                                                     DistanceFunction::Unknown);
        return *this;
    }

    LshParameterSetter& setHashes(int funcs, int tables, std::string storage) {
        _p.k = funcs;
        _p.l = tables;
        _p.storage_hash_table =
            get<StorageHashTable>(storageTypes, storage, StorageHashTable::Unknown);
        return *this;
    }

    LshParameterSetter& setFamily(std::string family) {
        _p.lsh_family = get<LSHFamily>(families, family, LSHFamily::Unknown);
        return *this;
    }

    void dump() {
        Rcout << "Sizes: " << _n << ", " << _d <<
            " (" << _p.dimension << ")" << std::endl;
        Rcout << "Hash k,l: " << _p.k << ", " << _p.l << std::endl;
        Rcout << "Distance: ";
        if ( _p.distance_function == DistanceFunction::EuclideanSquared ) {
            Rcout << "Euclidean Squared";
        } else if ( _p.distance_function == DistanceFunction::NegativeInnerProduct ) {
            Rcout << "Negative Inner Product";
        } else {
            Rcout << "unknown";
        }
        Rcout << std::endl;
            
        Rcout << "Family: ";
        if ( _p.lsh_family == LSHFamily::CrossPolytope ) {
            Rcout << "Cross Polytope";
        } else if ( _p.lsh_family == LSHFamily::Hyperplane ) {
            Rcout << "Hyperplane";
        } else {
            Rcout << "unknown";
        }
        Rcout << std::endl;
    }

  private:

    int                       _n;
    int                       _d;
    LSHConstructionParameters _p;
};

const familiesMap LshParameterSetter::families = {
    {"unknown",        LSHFamily::Unknown},
    {"hyperplane",     LSHFamily::Hyperplane},
    {"cross_polytope", LSHFamily::CrossPolytope}
};

const distancesMap LshParameterSetter::distances = {
    {"unknown",                DistanceFunction::Unknown},
    {"negative_inner_product", DistanceFunction::NegativeInnerProduct},
    {"euclidean_squared",      DistanceFunction::EuclideanSquared}
};

const storageTypesMap LshParameterSetter::storageTypes = {
    {"unknown",                    StorageHashTable::Unknown},
    {"flat_hash_table",            StorageHashTable::FlatHashTable},
    {"bit_packed_flat_hash_table", StorageHashTable::BitPackedFlatHashTable},
    {"stl_hash_table",             StorageHashTable::STLHashTable},
    {"linear_probing_hash_table",  StorageHashTable::LinearProbingHashTable}
};



//// 
////
//// [[Rcpp::export]]
//SEXP lsh_parameters_default(SEXP n_, SEXP d_) {
//    LshParameterSetter* p = new LshParameterSetter(as<int>(n_), as<int>(d_));
// 
//    return Rcpp::XPtr<LshParameterSetter> ptr(p, true);
//}
// 
//// 
////
//// [[Rcpp::export]]
//SEXP lsh_parameters_(SEXP n_, SEXP d_) {
//    LshParameterSetter* p = new LshParameterSetter(as<int>(n_), as<int>(d_));
// 
//    return Rcpp::XPtr<LshParameterSetter> ptr(p, true);
//}


RCPP_EXPOSED_CLASS(LshParameterSetter)

RCPP_MODULE(mod_params) {
    class_<LshParameterSetter>("LshParameterSetter")

    .constructor<int,int>()

    .method("setDefaults", &LshParameterSetter::setDefaults,
            "Fill with defaults")
    .method("setDistance", &LshParameterSetter::setDistance,
            "Set distance metric")
    .method("setHashes",   &LshParameterSetter::setHashes,
            "Set numbers of hash functions and tables")
    .method("setFamily",   &LshParameterSetter::setFamily,
            "Set LSH Family")
    .method("dump",        &LshParameterSetter::dump,
            "Check the values")
   ;
}
