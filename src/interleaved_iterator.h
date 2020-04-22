#ifndef OVERLAPPING_SEGMENTS_H
#define OVERLAPPING_SEGMENTS_H

#include <cstddef>
#include <cassert>
#include <utility>
#include <iterator>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <limits>

// Class used to support iteration (with OverlappingSegments::iterator)
// over pairs of overlapping segments within two sequences.
//
template <typename Container>
class OverlappingSegments {
  // short alias for the value type of the container
  typedef typename std::remove_reference<typename Container::value_type>::type CvalT;
  typedef long int IndexT;

  // references the two sequences
  const Container& loc1;
  const Container& loc2;

public:
  OverlappingSegments(Container& loc1_, Container& loc2_)
    : loc1(loc1_), loc2(loc2_) {
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

    IndexT Nx, Ny;

    // the minimum and maximum contained values;
    CvalT minloc, maxloc;

    // the current iteration state
    std::pair<IndexT, IndexT> idx;

  public:
    typedef std::pair<IndexT, IndexT> value_type;
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
      // the corresponding binwidth
      maxloc = std::numeric_limits<double>::infinity();
    }

    CvalT get_loc1(IndexT i) const {
      if (i < 0) return minloc;
      else if (i >= Nx) return maxloc;
      else return segs->loc1[i];
    }

    CvalT get_loc2(IndexT i) const {
      if (i < 0) return minloc;
      else if (i >= Ny) return maxloc;
      else return segs->loc2[i];
    }

    // Does interval i (from the first collection of segments) overlap
    // interval j (from the second)?
    bool intervals_overlap(IndexT i, IndexT j) const {
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

    iterator& operator++(int) {
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

#endif // OVERLAPPING_SEGMENTS_H
