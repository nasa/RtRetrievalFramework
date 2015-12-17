#ifndef BIN_MAP_H
#define BIN_MAP_H
#include <map>

namespace FullPhysics {
/****************************************************************//**
  This is a map that takes and index and returns a value. The index is
  a double, and we find the bin that it maps to.  This extrapolates,
  so if a index is greater the last bin than we return the last bin.
*******************************************************************/
template<class T> class BinMap {
public:
  typedef typename std::map<double, T>::const_iterator const_iterator;
  typedef typename std::map<double, T>::iterator iterator;
  const_iterator begin() const {return val.begin();}
  iterator begin() {return val.begin();}
  const_iterator end() const {return val.end();}
  iterator end() {return val.end();}
  
  BinMap() {}

//-----------------------------------------------------------------------
/// Constructor. This takes an iterator to fill in the boundaries of
/// the bins.
//-----------------------------------------------------------------------
  template<class It> BinMap(It xstart, It xend, const T& Initial_value)
  {
    if(xstart == xend)		// Handle degenerate case
      return;
    xstart++;			// Don't care what the lower edge is,
				// since we extrapolate anyways.
    while(xstart != xend)
      val.insert(std::pair<double, T>(*xstart++, Initial_value));
  }

//-----------------------------------------------------------------------
/// Constructor. This takes an iterator to fill in the boundaries of
/// the bins, and a function that is used to create the initial value.
//-----------------------------------------------------------------------
  template<class It, class C> 
  BinMap(It xstart, It xend, const C& Initial_value_creator)
  {
    if(xstart == xend)		// Handle degenerate case
      return;
    xstart++;			// Don't care what the lower edge is,
				// since we extrapolate anyways.
    while(xstart != xend)
      val.insert(std::pair<double, T>(*xstart++, Initial_value_creator()));
  }

//-----------------------------------------------------------------------
/// Return bin that "x" falls into.
//-----------------------------------------------------------------------

  const T& operator[](double x) const
  {
    typedef typename 
      std::map<double, T>::const_iterator Itype;
    Itype i = val.lower_bound(x);
    if(i == val.end())
      return val.rbegin()->second;
    else
      return i->second;
  }
  T& operator[](double x)
  {
    typedef typename 
      std::map<double, T>::iterator Itype;
    Itype i = val.lower_bound(x);
    if(i == val.end())
      return val.rbegin()->second;
    else
      return i->second;
  }

//-----------------------------------------------------------------------
// Return all the values. This is ordered by the bins.
//-----------------------------------------------------------------------

  std::vector<T> value() const
  {
    std::vector<T> res;
    typedef typename std::map<double, T>::value_type vt;
    for(typename std::map<double, T>::const_iterator i = val.begin();
	i != val.end(); ++i)
      res.push_back(i->second);
    return res;
  }
private:
  std::map<double, T> val;
};
}
#endif
