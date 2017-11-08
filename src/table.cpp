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


//' Constructs a LSH search table for a specified data set
//'
//' Note: The FALCONN LSH nearest-neighbor search uses a static table,
//' as contructed here. To add to or change the data set, a new
//' object needs to be constructed.
//'
//' @param tDataMatrix -- the data, an R numeric matrix where
//'                       each \emph{column} is a data point.
//'                       (Thus, pass the transpose of a typical data matrix.)
//'
//' @param params -- the LSH configuration parameters
//'
LshNnTable::LshNnTable(const NumericMatrix tDataMatrix,
                       const LshParameterSetter& params) {
    _params = params.params();
    _n_points = tDataMatrix.ncol();

    if ( tDataMatrix.nrow() != _params.dimension ) {
        stop("dimension mismatch between data matrix and LshTable parameters");
    }

    std::vector<double> asvec = Rcpp::as<std::vector<double> >(tDataMatrix);
    const double* data = asvec.data();

    DataPoints _data = {data,
                        static_cast<int_fast32_t>(_n_points),
                        _params.dimension};
    _table = std::shared_ptr<FnnTable>(construct_table<Point,int32_t,DataPoints>
                                       (_data, _params).release());
}

//' The dimension of the data points
//'
//' @return the dimension of points in the data matrix
int LshNnTable::dimension() const {
    return _params.dimension;
}

//' The number of data points
//'
//' @return the number points in the data matrix
int LshNnTable::size() const {
    return _n_points;
}

//' Find the data point nearest to the given query point
//'
//' Searches in the encapsulated data set with Locality-Sensitive Hashing.
//' This is an approximate nearest neighbor search.
//'
//' @param q -- query point, an R numeric vector of dimension d 
//'
//' @return index of the data point nearest to \code{q}, wrapped in
//'         an R integer vector
//'
IntegerVector LshNnTable::find_nearest_neighbor(const NumericVector& q) {
    EPoint query(Rcpp::as<EPoint>(q));
    int    nearest_index;

    nearest_index = _table->find_nearest_neighbor(query) + 1; // 1-indexed output

    return IntegerVector::create(nearest_index);
}

//' Find the data points nearest to the given query point
//'
//' Searches in the encapsulated data set with Locality-Sensitive Hashing.
//' This is an approximate nearest neighbor search.
//'
//' @param q -- query point, an R numeric vector of dimension d 
//' 
//' @param k -- the number of nearest neighbors to return
//'
//' @return vector of indicies of the data points nearest to \code{q},
//'         wrapped in an R integer vector
//'
IntegerVector LshNnTable::find_k_nearest_neighbors(const NumericVector& q,
                                                   int k) {
    EPoint        query(Rcpp::as<EPoint>(q));
    KeyVector     nearest_indices(k);
    IntegerVector nearest_indices_r;

    _table->find_k_nearest_neighbors(query, k, &nearest_indices);
    nearest_indices_r.assign(nearest_indices.begin(), nearest_indices.end());
    std::transform(nearest_indices_r.begin(),
                   nearest_indices_r.end(),
                   nearest_indices_r.begin(),
                   [](int index) { return index + 1; });
    return nearest_indices_r;
}

//' Find the data points within a specified radius of the given query point
//'
//' Searches in the encapsulated data set with Locality-Sensitive Hashing.
//' This is an approximate nearest neighbor search.
//'
//' @param q -- query point, an R numeric vector of dimension d 
//' 
//' @param radius -- radius around the query point in which to search
//'
//' @return vector of indicies of the data points within radius of \code{q},
//'         wrapped in an R integer vector
//'
IntegerVector LshNnTable::find_near_neighbors(const NumericVector& q,
                                              double radius) {
    EPoint        query(Rcpp::as<EPoint>(q));
    KeyVector     nearest_indices;
    IntegerVector nearest_indices_r;

    _table->find_near_neighbors(query, radius, &nearest_indices);
    nearest_indices_r.assign(nearest_indices.begin(), nearest_indices.end());
    std::transform(nearest_indices_r.begin(),
                   nearest_indices_r.end(),
                   nearest_indices_r.begin(),
                   [](int index) { return index + 1; });
    return nearest_indices_r;
}

//' Find all data points found in a single probing sequence
//'
//' This is a low-level operation. Note that a single data point might
//' appear multiple times in the result if several hash tables are used.
//' These duplicates are included in the returned value.
//'
//' @param q -- query point, an R numeric vector of dimension d 
//' 
//' @return vector of indicies of the data points found in the search,
//'         which may include duplicates, wrapped in an R integer vector
//'
IntegerVector LshNnTable::get_candidates(const NumericVector& q) {
    EPoint        query(Rcpp::as<EPoint>(q));
    KeyVector     nearest_indices;
    IntegerVector nearest_indices_r;

    _table->get_candidates_with_duplicates(query, &nearest_indices);
    nearest_indices_r.assign(nearest_indices.begin(), nearest_indices.end());

    return nearest_indices_r;
}

//' Find all data points found in a single probing sequence, duplicates removed
//'
//' This is a low-level operation. Note that a single data point might
//' appear multiple times in the result if several hash tables are used.
//' These duplicates are removed from the returned value.
//'
//' @param q -- query point, an R numeric vector of dimension d 
//' 
//' @return vector of indicies of the data points found in the search,
//'         which will not include duplicates, wrapped in an R integer vector
//'
IntegerVector LshNnTable::get_unique_candidates(const NumericVector& q) {
    EPoint        query(Rcpp::as<EPoint>(q));
    KeyVector     nearest_indices;
    IntegerVector nearest_indices_r;

    _table->get_unique_candidates(query, &nearest_indices);
    nearest_indices_r.assign(nearest_indices.begin(), nearest_indices.end());

    return nearest_indices_r;
}

//' Set the number of probes used in multi-probe LSH
//' 
//' Note that this is not a costly operation and can be
//' done repeatedly once the table has been constructed.
//'
//' @param num_probes -- the number of probes to use
//'
//' @return a reference to the table object to enable chaining
//' 
LshNnTable&   LshNnTable::setNumProbes(int num_probes) {
    _table->set_num_probes(num_probes);
    return *this;
}

//' Returns the number of probes currently being used for multi-probe LSH
//' 
//' @return the number of probes 
int           LshNnTable::getNumProbes() const {
    return _table->get_num_probes();
}

//' Set the maximum number of candidates to consider during similarity search
//' 
//' @param num_candidates -- the maximum number of candidates
//'
//' @return a reference to the table object to enable chaining
//' 
LshNnTable&   LshNnTable::setMaxNumCandidates(int num_candidates) {
    _table->set_max_num_candidates(num_candidates);
    return *this;
}

//' Returns the maximum number of candidates considered for each search
//' 
//' @return the maximum number of candidates
int           LshNnTable::getMaxNumCandidates() const {
    return _table->get_max_num_candidates();
}

//' Calculates the search accuracy for a specified number of probes
//'
//' @param queries -- matrix of queries, one per column
//'                   @seealso \code{\link{LshNnTable::tuneNumProbes}}
//'
//' @param answers -- vector of data indices,
//'                   @seealso \code{\link{LshNnTable::tuneNumProbes}}
//'
//' @param num_probes -- the number of probes to use for this check
//'
//' @return search accuracy for the given number of probes
//'
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

//' Find number of probes to achieve target precision on training data
//'
//' During the search for the optimal number of probes, this changes
//' the set number of probes via \code{setNumProbes()}. However, the
//' number originally set on entry is restored on exit.
//'
//' @param queries          -- matrix of query points, one point per
//'                            \emph{column} (pass transpose if necessary)
//' @param answers          -- vector of indices into data table representing
//'                            the nearest neighbor to the corresponding query
//' @param target_precision -- minimal correctness probability to achieve
//' @param init_num_probes  -- number of probes at which to start search
//' @param max_iterations   -- maximum number of iterations to consider;
//'                            default of -1 means no limit
//'                            
//' @return number of probes to use to achieve
//'
int           LshNnTable::tuneNumProbes(const NumericMatrix queries,
                                        IntegerVector answers,
                                        double target_precision,
                                        int init_num_probes,
                                        int max_iterations) {
    int original_num_probes = getNumProbes();
    int num_probes = init_num_probes;
    int iter = max_iterations;
  
    for ( ; iter != 0; --iter ) {
        double precision = computeProbePrecision(queries, answers, num_probes);
        if (precision >= target_precision) {
            break;
        }
        num_probes *= 2;
    }

    if ( iter == 0 ) {
        stop("maximum iterations exceeded while tuning number of probes");
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

    setNumProbes(original_num_probes);

    return num_probes;
}

/* RCPP_EXPOSED_CLASS(LshNnTable)
RCPP_EXPOSED_CLASS(LshParameterSetter) */

// Module mod_table exposes the LshNnTable class to R

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
