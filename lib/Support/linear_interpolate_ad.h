#ifndef LINEAR_INTERPOLATE_AD_H
#define LINEAR_INTERPOLATE_AD_H
#include "linear_interpolate.h"
#include "auto_derivative.h"
#include <iostream>

namespace FullPhysics {

template<> class InterpolatePoint<AutoDerivative<double>, 
				  AutoDerivative<double> > {
public:
  virtual const AutoDerivative<double>& x_min() const = 0;
  virtual const AutoDerivative<double>& x_max() const = 0;
  virtual void interpolate(const AutoDerivative<double>& x,
			   const AutoDerivativeRef<double>& res) const = 0;
};

template<> class Return1Point<AutoDerivative<double>, 
			      AutoDerivative<double> > : 
    public InterpolatePoint<AutoDerivative<double>, 
			    AutoDerivative<double> >,
    public Printable<Return1Point<AutoDerivative<double>, 
				  AutoDerivative<double> > > {
public:
  Return1Point(const AutoDerivative<double>& x0, const AutoDerivative<double>& y0)
    : x0_(x0), y0_(y0) 
  {}
  const AutoDerivative<double>& x_min() const {return x0_; }
  const AutoDerivative<double>& x_max() const {return x0_; }
  void interpolate(const AutoDerivative<double>& x,
		   const AutoDerivativeRef<double>& res) const { res = y0_; }
  void print(std::ostream& Os) const { Os << "Return1Point"; }
private:
  AutoDerivative<double> x0_;
  AutoDerivative<double> y0_;
};

/****************************************************************//**
  This is a specialization of LinearInterpolate2Point for AutoDerivative.
  This specialization should run faster.
*******************************************************************/

template<> 
class LinearInterpolate2Point<AutoDerivative<double>, AutoDerivative<double> > :
    public InterpolatePoint<AutoDerivative<double>, AutoDerivative<double> >,
    public Printable<LinearInterpolate2Point<AutoDerivative<double>, 
					     AutoDerivative<double> > > {
public:
  LinearInterpolate2Point(const AutoDerivative<double>& x0, 
			  const AutoDerivative<double> & y0, 
			  const AutoDerivative<double>& x1, 
			  const AutoDerivative<double> & y1)
    : x0_(x0), x1_(x1), y0_(y0), delta_y0_((y1 - y0) / (x1 - x0))
  {
    range_max_check(x0, x1);
  }
  const AutoDerivative<double>& x_min() const {return x0_; }
  const AutoDerivative<double>& x_max() const {return x1_; }
  virtual void interpolate(const AutoDerivative<double>& x,
			   const AutoDerivativeRef<double>& res) const
  {
    res.value_ref() = y0_.value() + 
      delta_y0_.value() * (x.value() - x0_.value());
    if(x.is_constant() && y0_.is_constant())
      res.gradient_ref() = 0;
    else if(x.is_constant())
      res.gradient_ref() = y0_.gradient() + delta_y0_.gradient() * 
	(x.value() - x0_.value()) - delta_y0_.value() * x0_.gradient();
    else if(y0_.is_constant())
      res.gradient_ref() = delta_y0_.value() * x.gradient();
    else
      res.gradient_ref() = y0_.gradient() + delta_y0_.gradient() * 
	(x.value() - x0_.value()) + 
	delta_y0_.value() * (x.gradient() - x0_.gradient());
  }
  void print(std::ostream& Os) const { Os << "LinearInterpolate2Point"; }
private:
  AutoDerivative<double> x0_;
  AutoDerivative<double> x1_;
  AutoDerivative<double> y0_;
  AutoDerivative<double> delta_y0_;
};

/****************************************************************//**
  This is a specialization of LinearInterpolate for AutoDerivative.
  This specialization should run faster.
*******************************************************************/

template<> class LinearInterpolate<AutoDerivative<double>, 
				   AutoDerivative<double> > : 
    public Printable<LinearInterpolate<AutoDerivative<double>, 
				       AutoDerivative<double> > > {
public:
  enum BehaviorOutOfRange 
    {OUT_OF_RANGE_EXTRAPOLATE = 0, OUT_OF_RANGE_CLIP, OUT_OF_RANGE_ERROR};
  LinearInterpolate() {}

//-----------------------------------------------------------------------
/// Constructor. This takes iterators to give x and y, along with the 
/// behavior when called with out of range data.
//-----------------------------------------------------------------------

  template<class I1, class I2> LinearInterpolate(I1 xstart, I1 xend,
      I2 ystart, BehaviorOutOfRange Out_of_range = OUT_OF_RANGE_EXTRAPOLATE)
    : out_of_range(Out_of_range), nvar(0)

  {
    if(std::distance(xstart, xend) == 0) {
      // If no values then we can not do any interpolation, 
      // error in operator if attempted
      return;
    } else {
      nvar = std::max(xstart->number_variable(), ystart->number_variable());
      if(std::distance(xstart, xend) == 1) {
      // If only one value available then just return that one point
      // no need for interpolation
      AutoDerivative<double> x0 = *xstart;
      AutoDerivative<double> y0 = *ystart;
      if(nvar > 0 && x0.is_constant()) {
	x0.gradient().resize(nvar);
	x0.gradient() = 0;
      }
      if(nvar > 0 && y0.is_constant()) {
	y0.gradient().resize(nvar);
	y0.gradient() = 0;
      }
      if(x0.number_variable() != y0.number_variable())
	throw Exception("x0 and y0 need  to have the same number of variables");
      inter[x0] = boost::shared_ptr<InterpolatePoint<AutoDerivative<double>, 
						     AutoDerivative<double> > >
	(new Return1Point<AutoDerivative<double>, 
			  AutoDerivative<double> >(x0, y0));

      } else while(xstart != xend) {
	  // For two more points we can use linear interpolation
	  AutoDerivative<double> x0 = *xstart++;
	  AutoDerivative<double> y0 = *ystart++;
	  if(nvar > 0 && x0.is_constant()) {
	    x0.gradient().resize(nvar);
	    x0.gradient() = 0;
	  }
	  if(nvar > 0 && y0.is_constant()) {
	    y0.gradient().resize(nvar);
	    y0.gradient() = 0;
	  }
	  if(x0.number_variable() != nvar)
	    throw Exception("All X and Y need to have the same number of variables");
	  if(y0.number_variable() != nvar)
	    throw Exception("All X and Y need to have the same number of variables");
      
	  if(xstart == xend)
	    break;
	  AutoDerivative<double> x1 = *xstart;
	  AutoDerivative<double> y1 = *ystart;
	  if(nvar > 0 && x1.is_constant()) {
	    x1.gradient().resize(nvar);
	    x1.gradient() = 0;
	  }
	  if(nvar > 0 && y1.is_constant()) {
	    y1.gradient().resize(nvar);
	    y1.gradient() = 0;
	  }
	  if(x1.number_variable() != nvar)
	    throw Exception("All X and Y need to have the same number of variables");
	  if(y1.number_variable() != nvar)
	    throw Exception("All X and Y need to have the same number of variables");
	  if(x1 <= x0)
	    throw Exception("X needs to be sorted");
      
	  inter[x1] = boost::shared_ptr<InterpolatePoint<AutoDerivative<double>, 
							 AutoDerivative<double> > >
	    (new LinearInterpolate2Point<AutoDerivative<double>, 
				     AutoDerivative<double> >(x0, y0, x1, y1));
	}
    }
  }
  AutoDerivative<double>  operator()(const AutoDerivative<double>& x) const
  { AutoDerivative<double> res;
    res.gradient().resize(std::max(nvar, x.number_variable()));
    AutoDerivativeRef<double> res2(res.value(), res.gradient());
    interpolate(x, res2);
    return res;
  }
  // Version of operator() that avoids creating a new temporary for
  // returning. Instead we use an existing ArrayRef. Note that this
  // should already be the right size, if it isn't an error will
  // occur.
  void interpolate(const AutoDerivative<double>& x, 
		   const AutoDerivativeRef<double>& res) const
  { 
    if(inter.empty()) {
      throw Exception("Can not interpolate since interpolation map is empty.");
    }

    typedef 
      std::map<AutoDerivative<double>, boost::shared_ptr<InterpolatePoint<AutoDerivative<double>, AutoDerivative<double> > > >::
      const_iterator Itype;
    Itype i = inter.lower_bound(x);
    if(i == inter.end()) { // Are we past the end?
      switch(out_of_range) {
      case OUT_OF_RANGE_ERROR: 
	{
	  std::stringstream err_msg;
	  err_msg << "Linear interpolation AD point "
		  << x.value() << " is past the end of interpolation range: "
		  << "[" << inter.begin()->first.value()
		  << ", " << inter.rbegin()->first.value() << "]";
	  throw Exception(err_msg.str());
	}
      case OUT_OF_RANGE_EXTRAPOLATE:
	inter.rbegin()->second->interpolate(x, res);
	break;
      case OUT_OF_RANGE_CLIP:
	inter.rbegin()->second->interpolate(inter.rbegin()->second->x_max(), 
					    res);
	break;
      default:
	throw Exception("Unknown out_of_range value");
      }
    } else {
      if(x >= i->second->x_min())
	i->second->interpolate(x, res);
      else 
	switch(out_of_range) {
	case OUT_OF_RANGE_ERROR:
	  {
	    std::stringstream err_msg;
	    err_msg << "Linear interpolation AD point "
		    << x.value() << " is outside of interpolation range: "
		    << "[" << inter.begin()->first.value()
		    << ", " << inter.rbegin()->first.value() << "]";
	    throw Exception(err_msg.str());
	  }
	case OUT_OF_RANGE_EXTRAPOLATE:
	  i->second->interpolate(x, res);
	  break;
	case OUT_OF_RANGE_CLIP:
	  i->second->interpolate(i->second->x_min(), res);
	  break;
	default:
	  throw Exception("Unknown out_of_range value");
	}
    }
  }
  void print(std::ostream& Os) const { Os << "LinearInterpolate"; }
private:
  std::map<AutoDerivative<double>, boost::shared_ptr<InterpolatePoint<AutoDerivative<double>, AutoDerivative<double> > > > inter;
  BehaviorOutOfRange out_of_range;
  int nvar;
};
}
#endif
