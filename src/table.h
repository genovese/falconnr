//! R-exposed proxy for the LSH Nearest-Neighbor table

#ifndef FALCONNR_TABLE_H
#define FALCONNR_TABLE_H

#include "falconnr.h"
#include "params.h"

using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Rcpp::IntegerVector;

using falconn::LSHConstructionParameters;
using falconn::PlainArrayPointSet;

class LshNnTable {
  public:
    typedef falconn::LSHNearestNeighborTable<Point> FnnTable;
    typedef std::shared_ptr<FnnTable>               FnnTablePtr;
    typedef PlainArrayPointSet<double>              DataPoints;

    LshNnTable(const NumericMatrix tDataMatrix,
               const LshParameterSetter& params);

    int dimension() const;
    int size() const;

    IntegerVector find_nearest_neighbor(const NumericVector& q);
    IntegerVector find_k_nearest_neighbors(const NumericVector& q, int k);
    IntegerVector find_near_neighbors(const NumericVector& q, double radius);

    IntegerVector get_candidates(const NumericVector& q);
    IntegerVector get_unique_candidates(const NumericVector& q);

    LshNnTable& setNumProbes(int num_probes);
    int         getNumProbes() const;
    int         tuneNumProbes(const NumericMatrix queries,
                              IntegerVector answers,
                              double target_precision,
                              int init_num_probes = 1,
                              int max_iterations = -1);

    LshNnTable& setMaxNumCandidates(int num_candidates = FnnTable::kNoMaxNumCandidates);
    int         getMaxNumCandidates() const;

  private:
    FnnTablePtr                 _table;
    LSHConstructionParameters   _params;
    int                         _n_points;

    double      computeProbePrecision(const NumericMatrix queries,
                                      IntegerVector answers,
                                      int num_probes);
};

#endif

// Local Variables:
// mode: c++
// End:
