#ifndef LOG_INTERPOLATE_H
#define LOG_INTERPOLATE_H
#include "linear_interpolate.h"
#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/foreach.hpp>
#include <cmath>

namespace FullPhysics {
/****************************************************************//**
  Wrapper around LinearInterpolate that uses log(x) rather than x 
  in interpolating.
*******************************************************************/
template<class TX, class TY> class LogLinearInterpolate : 
    public Printable<LogLinearInterpolate<TX,TY> > {
public:
  typedef typename 
  LinearInterpolate<TX, TY>::BehaviorOutOfRange BehaviorOutOfRange;

//-----------------------------------------------------------------------
/// Constructor. This takes iterators to give x and y, along with the 
/// behavior when called with out of range data.
//-----------------------------------------------------------------------

  LogLinearInterpolate() {};
  // Gcc 4.2 does not like having template default arugment. As an
  // easy workaround, use the hardcoded value for extrapolate
  // template<class I1, class I2> LogLinearInterpolate(I1 xstart, I1 xend,
  //   I2 ystart, BehaviorOutOfRange Out_of_range = 
  //   LinearInterpolate<TX, TY>::OUT_OF_RANGE_EXTRAPOLATE)
  template<class I1, class I2> LogLinearInterpolate(I1 xstart, I1 xend,
    I2 ystart, BehaviorOutOfRange Out_of_range = BehaviorOutOfRange(0))
  {
    std::vector<TX> xval(xstart, xend);
    BOOST_FOREACH(TX &x ,xval)
      x = std::log(x);
    interp = LinearInterpolate<TX,TY>(xval.begin(), xval.end(), ystart, 
				      Out_of_range);
  }
  TY operator()(const TX& x) const
  { 
    return interp(std::log(x));
  }
  void print(std::ostream& Os) const { Os << "LogLinearInterpolate"; }
private:
  LinearInterpolate<TX,TY> interp;
};

/****************************************************************//**
  Wrapper around LinearInterpolate that uses log(x) and log(y) rather
  than x and y when interpolating.
*******************************************************************/
template<class TX, class TY> class LogLogInterpolate : 
    public Printable<LogLogInterpolate<TX, TY> > {
public:
  typedef typename 
  LinearInterpolate<TX, TY>::BehaviorOutOfRange BehaviorOutOfRange;
  LogLogInterpolate() {};
  // Gcc 4.2 does not like having template default arugment. As an
  // easy workaround, use the hardcoded value for extrapolate
  // template<class I1, class I2> LogLogInterpolate(I1 xstart, I1 xend,
  //    I2 ystart, BehaviorOutOfRange Out_of_range = 
  //    LinearInterpolate<TX, TY>::OUT_OF_RANGE_EXTRAPOLATE)
  template<class I1, class I2> LogLogInterpolate(I1 xstart, I1 xend,
     I2 ystart, BehaviorOutOfRange Out_of_range = BehaviorOutOfRange(0))
  {
    std::vector<TX> xval;
    std::vector<TY> yval;
    for(; xstart != xend; ++xstart, ++ystart) {
      xval.push_back(std::log(*xstart));
      yval.push_back(std::log(*ystart));
    }
    interp = LinearInterpolate<TX, TY>(xval.begin(), xval.end(), yval.begin(),
				       Out_of_range);
  }
  TY operator()(const TX& x) const
  { 
    return std::exp(interp(std::log(x)));
  }
  void print(std::ostream& Os) const { Os << "LogLogInterpolate"; }
private:
  LinearInterpolate<TX,TY> interp;
};

/****************************************************************//**
  Specialization for AutoDerivative<double> types
*******************************************************************/
template<> class LogLogInterpolate<AutoDerivative<double>, AutoDerivative<double> >: 
  public Printable<LogLogInterpolate<AutoDerivative<double>, AutoDerivative<double> > > {
public:
  typedef LinearInterpolate<AutoDerivative<double>, AutoDerivative<double> >::BehaviorOutOfRange BehaviorOutOfRange;
  LogLogInterpolate() {};
  // Gcc 4.2 does not like having template default arugment. As an
  // easy workaround, use the hardcoded value for extrapolate
  // template<class I1, class I2> LogLogInterpolate(I1 xstart, I1 xend,
  //    I2 ystart, BehaviorOutOfRange Out_of_range = 
  //    LinearInterpolate<AutoDerivative<double>, AutoDerivative<double> >::OUT_OF_RANGE_EXTRAPOLATE)
  template<class I1, class I2> LogLogInterpolate(I1 xstart, I1 xend,
     I2 ystart, BehaviorOutOfRange Out_of_range = BehaviorOutOfRange(0))
  {
    std::vector<AutoDerivative<double> > xval;
    std::vector<AutoDerivative<double> > yval;
    for(; xstart != xend; ++xstart, ++ystart) {
      xval.push_back(std::log(*xstart));
      yval.push_back(std::log(*ystart));
    }
    interp = LinearInterpolate<AutoDerivative<double>, AutoDerivative<double> >(xval.begin(), xval.end(), yval.begin(),
				       Out_of_range);
  }
  AutoDerivative<double> operator()(const AutoDerivative<double> & x) const
  { 
    return std::exp(interp(std::log(x)));
  }

  void interpolate(const AutoDerivative<double>& x, 
		   const AutoDerivativeRef<double>& res) const
  {
    interp.interpolate(std::log(x), res);
    res = std::exp(res);
  }

  void print(std::ostream& Os) const { Os << "LogLogInterpolate"; }
private:
  LinearInterpolate<AutoDerivative<double>, AutoDerivative<double> > interp;
};

/****************************************************************//**
  Wrapper around LinearInterpolate that uses log(y) rather
  than x and y when interpolating.
*******************************************************************/
template<class TX, class TY> class LinearLogInterpolate : 
    public Printable<LinearLogInterpolate<TX, TY> > {
public:
  typedef typename 
  LinearInterpolate<TX, TY>::BehaviorOutOfRange BehaviorOutOfRange;
  LinearLogInterpolate() {};
  // Gcc 4.2 does not like having template default arugment. As an
  // easy workaround, use the hardcoded value for extrapolate
  // template<class I1, class I2> LinearLogInterpolate(I1 xstart, I1 xend,
  //    I2 ystart, BehaviorOutOfRange Out_of_range = 
  //    LinearInterpolate<TX, TY>::OUT_OF_RANGE_EXTRAPOLATE)
  template<class I1, class I2> LinearLogInterpolate(I1 xstart, I1 xend,
     I2 ystart, BehaviorOutOfRange Out_of_range = BehaviorOutOfRange(0))
  {
    std::vector<TX> xval;
    std::vector<TY> yval;
    for(; xstart != xend; ++xstart, ++ystart) {
      xval.push_back(*xstart);
      yval.push_back(std::log(*ystart));
    }
    interp = LinearInterpolate<TX, TY>(xval.begin(), xval.end(), yval.begin(),
				       Out_of_range);
  }
  TY operator()(const TX& x) const
  { 
    return std::exp(interp(x));
  }
  void print(std::ostream& Os) const { Os << "LinearLogInterpolate"; }
private:
  LinearInterpolate<TX,TY> interp;
};

/****************************************************************//**
  Specialization for AutoDerivative<double> types
*******************************************************************/
template<> class LinearLogInterpolate<AutoDerivative<double>, AutoDerivative<double> >: 
  public Printable<LinearLogInterpolate<AutoDerivative<double>, AutoDerivative<double> > > {
public:
  typedef LinearInterpolate<AutoDerivative<double>, AutoDerivative<double> >::BehaviorOutOfRange BehaviorOutOfRange;
  LinearLogInterpolate() {};
  // Gcc 4.2 does not like having template default arugment. As an
  // easy workaround, use the hardcoded value for extrapolate
  // template<class I1, class I2> LinearLogInterpolate(I1 xstart, I1 xend,
  //    I2 ystart, BehaviorOutOfRange Out_of_range = 
  //    LinearInterpolate<AutoDerivative<double>, AutoDerivative<double> >::OUT_OF_RANGE_EXTRAPOLATE)
  template<class I1, class I2> LinearLogInterpolate(I1 xstart, I1 xend,
     I2 ystart, BehaviorOutOfRange Out_of_range = BehaviorOutOfRange(0))
  {
    std::vector<AutoDerivative<double> > xval;
    std::vector<AutoDerivative<double> > yval;
    for(; xstart != xend; ++xstart, ++ystart) {
      xval.push_back(*xstart);
      yval.push_back(std::log(*ystart));
    }
    interp = LinearInterpolate<AutoDerivative<double>, AutoDerivative<double> >(xval.begin(), xval.end(), yval.begin(),
										Out_of_range);
  }
  AutoDerivative<double> operator()(const AutoDerivative<double> & x) const
  { 
    return std::exp(interp(x));
  }

  void interpolate(const AutoDerivative<double>& x, 
		   const AutoDerivativeRef<double>& res) const
  {
    interp.interpolate(x, res);
    res = std::exp(res);
  }

  void print(std::ostream& Os) const { Os << "LinearLogInterpolate"; }
private:
  LinearInterpolate<AutoDerivative<double>, AutoDerivative<double> > interp;
};


}
#endif
