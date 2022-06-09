// Enable C++11
// [[Rcpp::plugins(cpp11)]]

#ifndef FASTSMOOTHV2_H
#define FASTSMOOTHV2_H

#include <Rcpp.h>
#include <cstddef>
#include <cassert>
#include <utility>
#include <iterator>
#include <stdexcept>
#include <limits>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include <math.h>

using namespace Rcpp;

class OverlappingSegments {

  // The two sequences of left-hand segment endpoints
  const NumericVector& loc1;
  const NumericVector& loc2;
  const double binWidth1, binWidth2;

  // shorter names for loc1.size() and loc2.size()
  const long N1, N2;

  // the minimum and maximum values over both sequences;
  double minloc, maxloc;

public:
  OverlappingSegments(NumericVector& loc1_, NumericVector& loc2_,
                      double binWidth1_ = 1.0, double binWidth2_ = 1.0)
    : loc1(loc1_), loc2(loc2_),
      binWidth1(binWidth1_), binWidth2(binWidth2_),
      N1(loc1.size()), N2(loc2.size())
  {
    if (N1 == 0 || N2 == 0)
      throw std::invalid_argument("Input vectors must be nonempty");

    for (int i = 0; i < N1 - 1; i++)
      if (loc1[i] > loc1[i + 1] + binWidth1)
        throw std::invalid_argument(
          "Elements of loc1 must be sorted in ascending order, "
          "with elements separated by at least binWidth1");

    for (int i = 0; i < N2 - 1; i++)
      if (loc2[i] > loc2[i + 1] + binWidth2)
        throw std::invalid_argument(
          "Elements of loc2 must be sorted in ascending order, "
          "with elements separated by at least binWidth2");

    minloc = std::min(loc1[0], loc2[0]);
    maxloc = std::max(loc1[N1 - 1] + binWidth1, loc2[N2 - 1] + binWidth2);
  }

  // left, mid and right locations of the segments, for loc1 and loc2
  double loc1_left(long i) const {
    return (i >= 0) ? loc1[i] : minloc;
  }

  double loc1_mid(long i) const {
    return (i >= 0) ? (loc1[i] + binWidth1) : minloc;
  }

  double loc1_right(long i) const {
    return (i + 1 < N1) ? loc1[i + 1] : maxloc;
  }

  double loc2_left(long i) const {
    return (i >= 0) ? loc2[i] : minloc;
  }

  double loc2_mid(long i) const {
    return (i >= 0) ? (loc2[i] + binWidth2) : minloc;
  }

  double loc2_right(long i) const {
    return (i + 1 < N2) ? loc2[i + 1] : maxloc;
  }

  // Does interval i (from the first collection of segments) overlap
  // interval j (from the second)?
  bool intervals_overlap(long i, long j) const {
    return (loc1_left(i) < loc2_right(j) && loc2_left(j) < loc1_right(i));
  }

  // OverlappingSegments iterator
  class iterator {
    const OverlappingSegments& segs;

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
    explicit iterator(const OverlappingSegments& segs_)
      : segs(segs_)
    {
      if (segs.loc1[0] < segs.loc2[0]) {
        idx.first = 0;
        idx.second = -1;
      }
      else if (segs.loc1[0] == segs.loc2[0]) {
        idx.first = 0;
        idx.second = 0;
      }
      else {
        idx.first = -1;
        idx.second = 0;
      }
    }

    // Is the current iterator at one-past-the-end?  Equivalent to an
    // equality comparison with segs.end().
    bool at_end() const {
      return idx.first == segs.N1 && idx.second == segs.N2 - 1;
    }

    // Update the current iterator to point to one-past-the-end
    iterator& advance_to_end() {
      idx.first = segs.N1;
      idx.second = segs.N2 - 1;
      return *this;
    }

    iterator& operator++() {
#if !NDEBUG
      // Verify precondition
      if (!segs.intervals_overlap(idx.first, idx.second)) {
        throw std::logic_error("Iterator precondition not satisfied: "
                               "current intervals do not overlap");
      }
#endif

      // Advance the second segment if it would still overlap the first
      //
      // The condition below is equivalent to
      //   idx.second < N2 - 1 && intervals_overlap(idx.first, idx.second + 1)
      // given that we know (by the precondition) that
      //   loc1_left(idx.first) < loc2_right(idx.second)
      // and therefore that
      //   loc1_left(idx.first) < loc2_right(idx.second + 1),
      //
      if (idx.second < segs.N2 - 1
          && segs.loc2_left(idx.second + 1) < segs.loc1_right(idx.first)) {
        idx.second++;
      }
      // Could not advance the second segment above: advance the first instead,
      // and the second as well if they share an endpoint
      else {
        if (idx.second < segs.N2 - 1
            && segs.loc2_left(idx.second + 1) == segs.loc1_right(idx.first)) {
          idx.second++;
        }
        idx.first++;
      }

#if !NDEBUG
      // Verify postcondition
      if (!(at_end() || segs.intervals_overlap(idx.first, idx.second))) {
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

    value_type operator*() const { return idx; }

    const value_type *operator->() const { return &idx; }

    friend bool operator==(const iterator& lhs, const iterator& rhs) {
      return lhs.idx == rhs.idx;
    }

    friend bool operator!=(const iterator& lhs, const iterator& rhs) {
      return !(lhs == rhs);
    }
  };

  iterator begin() { return iterator(*this); }
  iterator end() { return iterator(*this).advance_to_end(); }  
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

double netemd_smooth(NumericVector loc1, NumericVector val1, double binWidth1,
                      NumericVector loc2, NumericVector val2, double binWidth2);

#endif // FASTSMOOTHV2_H
