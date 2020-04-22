// Enable C++11
// [[Rcpp::plugins(cpp11)]]

#ifndef FASTSMOOTHV2_H
#define FASTSMOOTHV2_H

#include <Rcpp.h>
#include <cstddef>
#include <cassert>
#include <utility>
#include <iterator>
// #include <vector>
// #include <iostream>
#include <stdexcept>
#include <limits>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include <math.h>

using namespace Rcpp;

class OverlappingSegments {
  // references the two sequences
  const NumericVector& loc1;
  const NumericVector& loc2;
  const double binWidth1, binWidth2;

public:
  OverlappingSegments(NumericVector& loc1_, NumericVector& loc2_,
                      double binWidth1_ = 1.0, double binWidth2_ = 1.0)
    : loc1(loc1_), loc2(loc2_), binWidth1(binWidth1_), binWidth2(binWidth2_)
  {
    // check the requirement that loc1 and loc2 are nonempty and
    // (strictly) sorted

    if (loc1.size() == 0 || loc2.size() == 0)
      throw std::invalid_argument("Input vectors must be nonempty");
    
    for (int i = 0; i < loc1.size() - 1; i++)
      if (loc1[i] > loc1[i + 1])
        throw std::invalid_argument("Input vectors must be sorted in strict ascending order");

    for (int i = 0; i < loc2.size() - 1; i++)
      if (loc2[i] > loc2[i + 1])
        throw std::invalid_argument("Input vectors must be sorted in strict ascending order");
  }

  class iterator {
    const OverlappingSegments * const segs;

    long Nx, Ny;

    // the minimum and maximum contained values;
    double minloc, maxloc;

    // the current iteration state
    std::pair<long, long> idx;

  public:
    typedef std::pair<long, long> value_type;
    typedef void difference_type;
    typedef value_type* pointer;
    typedef value_type& reference;
    typedef std::input_iterator_tag iterator_category;

    // Iterate over pairs of indices (i,j) into the sequences loc1 and
    // loc2, where the intervals [loc1[i], loc1[i+1]] and [loc2[j],
    // loc2[j+1]] overlap.
    //
    // These indices are returned from the iterator as
    //   std::pair<long, long>.
    //
    // A sequence has an implicit segment from minloc (with index -1)
    // to its zeroth element.  The elements loc1[0] and loc2[0] are
    // compared to determine whether, for either sequence, this
    // initial implicit segment overlaps the zeroth segment of the
    // other one.  If both sequences start with the same value, the
    // iteration starts at (0,0).
    //
    explicit iterator(OverlappingSegments *segs_)
      : segs(segs_), Nx(segs_->loc1.size()), Ny(segs_->loc2.size())
    {
      // The initial state of the iterator
      
      if (segs->loc1[0] < segs->loc2[0]) {
        minloc = segs->loc1[0];
        idx.first = 0;
        idx.second = -1;
      }
      else if (segs->loc1[0] == segs->loc2[0]) {
        minloc = segs->loc1[0]; // == segs->loc2[0]
        idx.first = 0;
        idx.second = 0;
      }
      else {
        minloc = segs->loc2[0];
        idx.first = -1;
        idx.second = 0;
      }

      // maxloc is used to signal the end of the current segment:
      // since we consider only left-hand endpoints, the only
      // requirement is that it compares greater than any real loc in
      // either sequence - in reality it will be the largest loc plus
      // the corresponding bin width
      // maxloc = std::numeric_limits<double>::infinity();
      maxloc = std::max(segs->loc1[segs->loc1.size() - 1] + segs->binWidth1,
                        segs->loc2[segs->loc2.size() - 1] + segs->binWidth2);
    }

    double get_loc1(long i) const {
      if (i < 0) return minloc;
      else if (i >= Nx) return maxloc;
      else return segs->loc1[i];
    }

    double get_loc2(long i) const {
      if (i < 0) return minloc;
      else if (i >= Ny) return maxloc;
      else return segs->loc2[i];
    }

    // Does interval i (from the first collection of segments) overlap
    // interval j (from the second)?
    bool intervals_overlap(long i, long j) const {
      return (get_loc1(i) < get_loc2(j + 1) && get_loc2(j) < get_loc1(i + 1));
    }

    bool at_end() const { return idx.first == Nx && idx.second == Ny - 1; }

    iterator& advance_to_end() {
      idx.first = Nx;
      idx.second = Ny - 1;
      return *this;
    }
    
    value_type operator*() const { return idx; }

    const value_type *operator->() const { return &idx; }

    iterator& operator++() {

#if !NDEBUG
      // Verify precondition
      if (!intervals_overlap(idx.first, idx.second)) {
        throw std::logic_error("Iterator precondition not satisfied: "
                               "current intervals do not overlap");
      }
#endif

      // Advance the second segment if it would still overlap the first
      // 
      // The condition below is equivalent to
      //   idx.second < Ny - 1 && intervals_overlap(idx.first, idx.second + 1)
      //
      // For the right-hand condition, we know that (by precondition)
      //   get_loc1(idx.first) < get_loc2(idx.second + 1)
      //
      // and therefore that
      //   get_loc1(idx.first) < get_loc2(idx.second + 2),
      //
      // so this inequality is all that remains to be checked.
      //
      if (idx.second < Ny - 1
          && get_loc2(idx.second + 1) < get_loc1(idx.first + 1)) {
        idx.second++;
      }
      // Could not advance the second segment above: advance the first instead,
      // and the second as well if they share an endpoint
      else {
        if (idx.second < Ny - 1
            && get_loc1(idx.first + 1) == get_loc2(idx.second + 1)) {
          idx.second++;
        }
        idx.first++;
      }

#if !NDEBUG
      // Verify postcondition
      if (!(at_end() || intervals_overlap(idx.first, idx.second))) {
        throw std::logic_error("Iterator postcondition not satisfied: "
                               "current intervals do not overlap (not at end)");
      }
#endif

      return *this;
    }

    iterator operator++(int) {
      iterator res = *this;
      operator++();
      return res;
    }

    friend bool operator==(const iterator& lhs, const iterator& rhs) {
      return lhs.idx == rhs.idx;
    }

    friend bool operator!=(const iterator& lhs, const iterator& rhs) {
      return !(lhs == rhs);
    }
  };

  iterator begin() { return iterator(this); };
  iterator end() { return iterator(this).advance_to_end(); };
};


double bowtie_area(double length, double val1_start, double val1_end,
                          double val2_start, double val2_end);
double get_segment(double start, double end, double val1_start,
                          double val1_end, double val2_start, double val2_end);
double get_segment_constrained(double seg1L1, double seg1L2,
                                      double seg2L1, double seg2L2,
                                      double seg1V1, double seg1V2,
                                      double seg2V1, double seg2V2);

double get_double_segment_constrained(
  double seg1Loc1, double seg1Loc2, double seg1Loc3,
  double seg1Val1, double seg1Val2,
  double seg2Loc1, double seg2Loc2, double seg2Loc3,
  double seg2Val1, double seg2Val2);

inline double leftmost_segments(const NumericVector& loc1,
                                const NumericVector& loc2,
                                const NumericVector& val1,
                                const NumericVector& val2,
                                double binWidth1,
                                double maxLoc);

double NetEmdSmoothV2(NumericVector loc1, NumericVector val1, double binWidth1,
                      NumericVector loc2, NumericVector val2, double binWidth2);

#endif // FASTSMOOTHV2_H
