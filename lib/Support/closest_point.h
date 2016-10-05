#ifndef CLOSEST_POINT
#define CLOSEST_POINT

#include <blitz/array.h>

namespace FullPhysics {


//  Does not assume sorted array.
//  Returns the index of the element in
//  target that is closest to data.
template <class T1, class T2> 
int closest_point(
  const blitz::Array<T1, 1>& target,
  T2 data)
{
  return minIndex(abs(target-data))(0);
}


//  Does not assume sorted array.
//  Returns the indexes of the element in
//  target that are closest to the elements
//  in data.
template <class T1, class T2> 
blitz::Array<int, 1> closest_point(
  const blitz::Array<T1, 1>& target,
  const blitz::Array<T2, 1>& data)
{
  blitz::Array<int, 1> a(data.shape());
  for(int i=a.base(blitz::firstDim); i<a.base(blitz::firstDim)+a.rows(); i++)
    a(i) = closest_point(target, data(i));
    
  return a;
}


}

#endif
