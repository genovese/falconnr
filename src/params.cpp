#include "falconnr.h"
#include "params.h"

using namespace Rcpp;

using falconn::DistanceFunction;
using falconn::LSHFamily;
using falconn::StorageHashTable;
using falconn::LSHConstructionParameters;
using falconn::get_default_parameters;


/// Get value from string-keyed map with default value if key is missing
///  
/// @param m    map from strings to V values
/// @param key  string key
/// @param def  default value to be returned when key is not in map
///
/// @return key's entry in map or def if not
/// 
template <typename V>
V get( const std::map<std::string,V> &m, const std::string &key, const V &def) {
    typename std::map<std::string,V>::const_iterator it = m.find( key );
    if ( it == m.end() ) {
        return def;
    } else {
        return it->second;
    }
}

/// Get key from string-keyed map matching given value, or given label if none
///
/// Assumes implicitly that the map represents a one-to-one mapping.
///  
/// @param m     map from strings to V values
/// @param value a value in the map
/// @param label string to return if no matching value is found
///
/// @return value's key in map or label if no match
/// 
template <typename V>
std::string
invert( const std::map<std::string,V> &m, const V &value, std::string label) {
    using map_it = typename std::map<std::string,V>::const_iterator;

    for ( map_it it = m.begin(); it != m.end(); ++it ) {
        if ( it->second == value ) {
            return it->first;
        }
    }

    return label;    
}


// Maps taking strings from R to FALCONN enumerated types

const LshParameterSetter::distancesMap LshParameterSetter::distances = {
    {"unknown",                DistanceFunction::Unknown},
    {"negative_inner_product", DistanceFunction::NegativeInnerProduct},
    {"euclidean_squared",      DistanceFunction::EuclideanSquared}
};

const LshParameterSetter::storageTypesMap LshParameterSetter::storageTypes = {
    {"unknown",                    StorageHashTable::Unknown},
    {"flat_hash_table",            StorageHashTable::FlatHashTable},
    {"bit_packed_flat_hash_table", StorageHashTable::BitPackedFlatHashTable},
    {"stl_hash_table",             StorageHashTable::STLHashTable},
    {"linear_probing_hash_table",  StorageHashTable::LinearProbingHashTable}
};

const LshParameterSetter::familiesMap LshParameterSetter::families = {
    {"unknown",        LSHFamily::Unknown},
    {"hyperplane",     LSHFamily::Hyperplane},
    {"cross_polytope", LSHFamily::CrossPolytope}
};


//' Construct an instance of an \code{LshParameterSetter}
//'
//' @param n  -- number of data points
//' @param d  -- dimension of the data points
//'
LshParameterSetter::LshParameterSetter(int n, int d) : _n(n), _d(d) {
    withDefaults();
}


//' Returns a copy of the underlying parameters structure
//'
//' @seealso \code{\link{falconn/lsh_nn_table.h}}, specifically
//' \code{struct LSHConstructionParameters} for details of the returned 
//' structure.
//'
//' @return instance of \code{struct LSHConstructionParameters} representing
//'         the underlying LSH parameters
//'
LSHConstructionParameters LshParameterSetter::params() const {
    return _p;
}

//' Sets all parameters to their default values based on data size
//'
//' @param distance -- one of the strings: "negative_inner_product", 
//'                    "euclidean_squared", or "unknown"; all other values
//'                    lead to a setting of "unknown".
//'
//' @return a reference to the original object, enabling chaining
//'
LshParameterSetter& LshParameterSetter::withDefaults(std::string distance) {
    _p = get_default_parameters<Point>(_n, _d,
                                       get<DistanceFunction>(distances,
                                                             distance,
                                                             DistanceFunction::Unknown),
                                       true);
    return *this;
}

//' Sets the distance function used in similarity search
//'
//' @param distance -- one of the strings: "negative_inner_product", 
//'                    "euclidean_squared", or "unknown"; all other values
//'                    lead to a setting of "unknown".
//'
//' @return a reference to the original object, enabling chaining
//'
LshParameterSetter&  LshParameterSetter::distance(std::string distance) {
    _p.distance_function = get<DistanceFunction>(distances,
                                                 distance,
                                                 DistanceFunction::Unknown);
    return *this;
}

//' Sets the number of hash functions used in locality-sensitive hashing
//'
//' @param tables -- a positive integer, the number of hash functions
//'                  to use in locality-sensitive hashing
//'
//' @return a reference to the original object, enabling chaining
//'
LshParameterSetter& LshParameterSetter::numHashFunctions(int funcs) {
    _p.k = funcs;
    return *this;
}

//' Sets the number of hash tables used in locality-sensitive hashing
//'
//' @param tables -- a positive integer, the number of hash tables
//'                  to use in locality-sensitive hashing
//' @return a reference to the original object, enabling chaining
//'
LshParameterSetter& LshParameterSetter::numHashTables(int tables) {
    _p.l = tables;
    return *this;
}

//' Sets the storage mode used in similarity search
//'
//' @param storage -- one of the strings: "flat_hash_table",
//'                   "bit_packed_flat_hash_table", "stl_hash_table",
//'                   "linear_probing_hash_table", or "unknown"; 
//'                   all other values lead to a setting of "unknown".
//'
//' @return a reference to the original object, enabling chaining
//'
LshParameterSetter& LshParameterSetter::storage(std::string storage) {
    _p.storage_hash_table =
        get<StorageHashTable>(storageTypes, storage, StorageHashTable::Unknown);
    return *this;
}

//' Sets the family used for locality-sensitive hashing
//'
//' @param family -- one of the strings: "hyperplane", "cross_polytope",
//'                  or "unknown"; all other values lead to a setting
//'                  of "unknown".
//'
//' @return a reference to the original object, enabling chaining
//'
LshParameterSetter& LshParameterSetter::family(std::string family) {
    _p.lsh_family = get<LSHFamily>(families, family, LSHFamily::Unknown);
    return *this;
}

//' Sets the number of pseudo-rotations used with cross-polytope hash
//'
//' For sparse data, a value of 2 is recommended; for dense data, 1.
//' 
//' @param rotations -- number of rotations to use
//'                  or "unknown"; all other values lead to a setting
//'                  of "unknown".
//'
//' @return a reference to the original object, enabling chaining
//'
LshParameterSetter& LshParameterSetter::rotations(int numRotations) {
    _p.num_rotations = numRotations;
    return *this;
}


//' Represents parameters as an R list
//'
//' @return R-list with names corresponding to the parameters
//'
Rcpp::List LshParameterSetter::asList() {
    return List::create(_["points"] = _n,
                        _["dimension"] = _d,
                        _["hashFunctions"] = _p.k,
                        _["hashTables"] = _p.l,
                        _["seed"] = _p.seed,
                        _["lshFamily"] = invert(families, _p.lsh_family,
                                                "unknown"),
                        _["distance"] = invert(distances, _p.distance_function,
                                               "unknown"),
                        _["storage"] = invert(storageTypes, _p.storage_hash_table,
                                              "unknown"),
                        _["rotations"] = _p.num_rotations,
                        _["threads"] = _p.num_setup_threads,
                        _["last_cp_dimension"] = _p.last_cp_dimension,
                        _["feature_hashing_dimension"] = _p.feature_hashing_dimension);
}


// Module mod_params exposes the LshParameterSetter class

/* RCPP_EXPOSED_CLASS(LshParameterSetter) // Rcpp updates broke module import. Drop this? */

RCPP_MODULE(mod_params) {
    class_<LshParameterSetter>("LshParameterSetter")

    .constructor<int,int>()

    .method("withDefaults", &LshParameterSetter::withDefaults,
            "Fill with defaults")
    .method("distance", &LshParameterSetter::distance,
            "Set distance metric")
    .method("numHashFunctions",   &LshParameterSetter::numHashFunctions,
            "Set numbers of hash functions and tables")
    .method("numHashTables",   &LshParameterSetter::numHashTables,
            "Set numbers of hash functions and tables")
    .method("storage",   &LshParameterSetter::storage,
            "Set LSH Hash Table Storage type")
    .method("family",   &LshParameterSetter::family,
            "Set LSH Family")
    .method("rotations",   &LshParameterSetter::rotations,
            "Number of rotations")
    .method("asList",        &LshParameterSetter::asList,
            "List of parameter values by name")
   ;
}
