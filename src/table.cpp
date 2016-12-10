//! Interface to LSH Nearest Neighbor Table for a data set

// [[Rcpp::depends(RcppEigen)]]


#include <iterator>

#include "table.h"
#include <RcppEigen.h>

using namespace Rcpp;

using falconn::construct_table;
using falconn::compute_number_of_hash_functions;
using falconn::LSHNearestNeighborTable;

// typedef Eigen::MatrixXd  Md;
// typedef Eigen::VectorXd  vd;
// typedef Eigen::Map<Md>   MMd;
// typedef Eigen::Map<Vd>   MVd;

typedef Eigen::Map<Point>    EPoint;
typedef std::vector<int32_t> KeyVector;

LshNnTable::LshNnTable(const NumericMatrix tDataMatrix,
                       const LshParameterSetter& params) {
    int n_points = tDataMatrix.ncol();
    // int dimension = tDataMatrix.nrow();
    // assert(params.dimension == dimension)

    std::vector<double> asvec = Rcpp::as<std::vector<double> >(tDataMatrix);
    const double* data = asvec.data();

    _params = params.params();
    _data = {data, static_cast<int_fast32_t>(n_points), _params.dimension};
    _table = construct_table<Point,int32_t,DataPoints>(_data, _params);
}

LshNnTable::LshNnTable(LshNnTable &&o) : _table(std::move(o._table)) {}

LshNnTable& LshNnTable::operator=(LshNnTable &&o) {
    if ( this != &o ) {
        _table = std::move(o._table);
    }
    return *this;
}


IntegerVector LshNnTable::find_nearest_neighbor(const NumericVector& q) {
    EPoint query(Rcpp::as<EPoint>(q));
    int    nearest_index;

    nearest_index = _table->find_nearest_neighbor(query);

    return IntegerVector::create(nearest_index);
}

IntegerVector LshNnTable::find_k_nearest_neighbors(const NumericVector& q,
                                                   int k) {
    EPoint        query(Rcpp::as<EPoint>(q));
    KeyVector     nearest_indices(k);
    IntegerVector nearest_indices_r;

    _table->find_k_nearest_neighbors(query, k, &nearest_indices);
    nearest_indices_r.assign(nearest_indices.begin(), nearest_indices.end());

    return nearest_indices_r;
}

IntegerVector LshNnTable::find_near_neighbors(const NumericVector& q,
                                              double radius) {
    EPoint        query(Rcpp::as<EPoint>(q));
    KeyVector     nearest_indices;
    IntegerVector nearest_indices_r;

    _table->find_near_neighbors(query, radius, &nearest_indices);
    nearest_indices_r.assign(nearest_indices.begin(), nearest_indices.end());

    return nearest_indices_r;
}

IntegerVector LshNnTable::get_candidates(const NumericVector& q) {
    EPoint        query(Rcpp::as<EPoint>(q));
    KeyVector     nearest_indices;
    IntegerVector nearest_indices_r;

    _table->get_candidates_with_duplicates(query, &nearest_indices);
    nearest_indices_r.assign(nearest_indices.begin(), nearest_indices.end());

    return nearest_indices_r;
}

IntegerVector LshNnTable::get_unique_candidates(const NumericVector& q) {
    EPoint        query(Rcpp::as<EPoint>(q));
    KeyVector     nearest_indices;
    IntegerVector nearest_indices_r;

    _table->get_unique_candidates(query, &nearest_indices);
    nearest_indices_r.assign(nearest_indices.begin(), nearest_indices.end());

    return nearest_indices_r;
}

LshNnTable&   LshNnTable::setNumProbes(int num_probes) {
    _table->set_num_probes(num_probes);
    return std::move(*this);
}

int           LshNnTable::getNumProbes() const {
    return _table->get_num_probes();
}

LshNnTable&   LshNnTable::setMaxNumCandidates(int num_candidates) {
    _table->set_max_num_candidates(num_candidates);
    return std::move(*this);
}

int           LshNnTable::getMaxNumCandidates() const {
    return _table->get_max_num_candidates();
}


double        LshNnTable::computeProbePrecision(const NumericMatrix queries,
                                                IntegerVector answers,
                                                int num_probes) {
    int num_cols = queries.ncol();
    int answer_index = 0;
    int num_matches = 0;
    std::vector<int32_t> candidates;

    _table->set_num_probes(num_probes);

    for ( int column = 0; column < num_cols; ++column ) {
        NumericVector col = queries(_, column); // ATTN: as fails with queries
        EPoint query = Rcpp::as<EPoint>(col);   // directly in place of col
        _table->get_candidates_with_duplicates(query, &candidates);
        for ( auto candidate : candidates ) {
            if ( candidate == answers[answer_index] ) {
                ++num_matches;
                break;
            }
        }
        ++answer_index;
    }
    return (num_matches + 0.0) / (num_cols + 0.0);
}

int           LshNnTable::tuneNumProbes(const NumericMatrix queries,
                                        IntegerVector answers,
                                        double target_precision,
                                        int init_num_probes,
                                        int max_iterations) {
  int    num_probes = init_num_probes;
  
  for ( ;; ) {
      double precision = computeProbePrecision(queries, answers, num_probes);
      if (precision >= target_precision) {
          break;
      }
      num_probes *= 2;
  }

  for ( int mid, lo = num_probes/2; num_probes - lo > 1; ) {
      mid = (num_probes + lo) / 2;
      double precision = computeProbePrecision(queries, answers, num_probes);

      if (precision >= target_precision) {
          num_probes = mid;
      } else {
          lo = mid;
      }
  }

  return num_probes;
}


RCPP_EXPOSED_CLASS(LshNnTable)
RCPP_EXPOSED_CLASS(LshParameterSetter)

RCPP_MODULE(mod_table) {
    class_<LshNnTable>("LshNnTable")

    .constructor<const NumericMatrix, const LshParameterSetter&>()

    .method("find_nearest_neighbor", &LshNnTable::find_nearest_neighbor,
            "Returns index of the (approximate) nearest neighbor to a given query point")
    .method("find_k_nearest_neighbors", &LshNnTable::find_k_nearest_neighbors,
            "Returns indices of the (approximate) k nearest neighbors to a given query point")
    .method("find_near_neighbors", &LshNnTable::find_near_neighbors,
            "Returns indices of (approximate) neighbors within a given radius of query point")

    .method("getNumProbes", &LshNnTable::getNumProbes,
            "Returns number of probes used for multi-probe LSH")
    .method("setNumProbes", &LshNnTable::setNumProbes,
            "Sets number of probes used for multi-probe LSH and returns self")

    .method("getMaxNumCandidates", &LshNnTable::getMaxNumCandidates,
            "Returns maximum number of candidates to consider in each query")
    .method("setNumProbes", &LshNnTable::setNumProbes,
            "Sets maximum number of candidates to consider in each query and returns self")
    .method("tuneNumProbes", &LshNnTable::tuneNumProbes,
            "Trains number of probes to target specified precision, returns number of probes")
    ;
}
