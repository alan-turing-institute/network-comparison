#ifndef INTERLEAVED_ITERATOR_H
#define INTERLEAVED_ITERATOR_H

#include <cstddef>
#include <utility>
#include <iterator>
#include <vector>
#include <iostream>

// Interleave two containers of ordered elements in a particular way.
//
// Supports iteration with Interleave::iterator
template <typename Container>
class Interleaved {
  // short alias for the value type of the container
  typedef typename Container::value_type CvalT;
  typedef long int IndexT;

  Container& xs;
  Container& ys;

public:
  Interleaved(Container& xs_, Container& ys_) : xs(xs_), ys(ys_) { }
  
  class iterator {
    const Interleaved * const xys;

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
          
    explicit iterator(Interleaved *xys_)
      : xys(xys_), Nx(xys_->xs.size()), Ny(xys_->ys.size())
    {
      //// this slightly simpler version assumes both xs and ys are
      //// nonempty:
      
      // if (xys->xs[0] <= xys->ys[0]) {
      //   minloc = xys->xs[0];
      //   idx.first = 0;
      //   idx.second = -1;
      // }
      // else {
      //   minloc = xys->ys[0];
      //   idx.first = -1;
      //   idx.second = 0;
      // }
      // maxloc = std::max(xys->xs.back(), xys->ys.back());

      //// no such assumption:
      
      if (Nx != 0 && Ny != 0) {
        maxloc = std::max(xys->xs.back(), xys->ys.back());
      }
      else if (Nx != 0) {
        maxloc = xys->xs.back();
      }
      else if (Ny != 0) {
        maxloc = xys->ys.back();
      }
      else {
        maxloc = 0;
        minloc = 0;
        idx.first = -1;
        idx.second = -1;
      }
      
      if (Nx != 0 && (Ny == 0 || Ny != 0 && xys->xs[0] <= xys->ys[0])) {
        minloc = xys->xs[0];
        idx.first = 0;
        idx.second = -1;
      }
      else if (Ny != 0) {
        minloc = xys->ys[0];
        idx.first = -1;
        idx.second = 0;
      }
    }

    CvalT get_x(IndexT i) const {
      if (i < 0) return minloc;
      else if (i >= Nx) return maxloc;
      else return xys->xs[i];
    }

    CvalT get_y(IndexT i) const {
      if (i < 0) return minloc;
      else if (i >= Ny) return maxloc;
      else return xys->ys[i];
    }
    
    std::pair<CvalT, CvalT> get_loc(const std::pair<IndexT, IndexT>& i) const {
      CvalT x(get_x(i.first));
      CvalT y(get_y(i.second));

      return std::make_pair(x, y);
    }

    std::pair<CvalT, CvalT> get_current_loc() const {
      return get_loc(idx);
    }

    iterator& advance_to_end() {
      idx.first = Nx;
      idx.second = Ny;
      return *this;
    }

    value_type operator*() const {
      return idx;
    }

    const value_type *operator->() const {
      return &idx;
    }

    iterator& operator++() {
      std::pair<CvalT, CvalT> loc(get_current_loc());
      while (idx.second < Ny && get_y(++idx.second) < loc.first) { }

      if (get_y(idx.second + 1) > get_x(idx.first + 1) || idx.second >= Ny) {
        if (idx.first < Nx) idx.first++;
      } else {
        if (idx.second < Ny) idx.second++;
      }
      
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

#endif // INTERLEAVED_ITERATOR_H
