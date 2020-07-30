#include "ils_convolution.h"
#include "hdf_file.h"
#include "ostream_pad.h"
#include <boost/lexical_cast.hpp>
using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(IlsConvolution, Ils)
.def(luabind::constructor<const boost::shared_ptr<Dispersion>&,
			  const boost::shared_ptr<IlsFunction>&>())
.def(luabind::constructor<const boost::shared_ptr<Dispersion>&,
			  const boost::shared_ptr<IlsFunction>&, DoubleWithUnit&>())
REGISTER_LUA_END()
#endif

// See base class for description.
blitz::Array<double, 1> IlsConvolution::apply_ils
(const blitz::Array<double, 1>& Hres_wn,
 const blitz::Array<double, 1>& Hres_rad,
 const std::vector<int>& Pixel_list) const
{
  if(Hres_wn.rows() != Hres_rad.rows())
    throw Exception("wave_number and radiance need to be the same size");
  Array<double, 1> disp_wn(disp->pixel_grid().data());
  ArrayAd<double, 1> response;
  Array<double,1> res((int) Pixel_list.size());
  for(int i = 0; i < res.rows(); ++i) {
    // Find the range of Hres_wn that lie within
    // wn_center +- ils_half_width

    double wn_center = disp_wn(Pixel_list[i]);
    Array<double, 1>::const_iterator itmin = 
      std::lower_bound(Hres_wn.begin(), Hres_wn.end(),
		       wn_center - ils_half_width_.value);
    Array<double, 1>::const_iterator itmax = 
      std::lower_bound(Hres_wn.begin(), Hres_wn.end(),
		       wn_center + ils_half_width_.value);
    int jmin = (itmin == Hres_wn.end() ? Hres_wn.rows() - 1 : itmin.position()(0));
    int jmax = (itmax == Hres_wn.end() ? Hres_wn.rows() - 1 : itmax.position()(0));
    Range r(jmin, jmax);

    // Convolve with response

    ils_func->ils(wn_center, Hres_wn(r), response);
    Array<double, 1> conv(response.value() * Hres_rad(r));

    // And integrate to get pixel value.

    res(i) = integrate(Hres_wn(r), conv) / 
      integrate(Hres_wn(r), response.value());
  }
  return res;
}


//-----------------------------------------------------------------------
/// Simple trapezoid rule for integration.
//-----------------------------------------------------------------------

double IlsConvolution::integrate(const blitz::Array<double, 1>& x, 
				 const blitz::Array<double, 1>& y) const
{
  double res = 0;
  for(int i = 1; i < x.rows(); ++i)
    res += (x(i) - x(i - 1)) * (y(i) + y(i - 1));
  return res / 2.0;
}

//-----------------------------------------------------------------------
/// Simple trapezoid rule for integration.
//-----------------------------------------------------------------------

AutoDerivative<double> IlsConvolution::integrate(const blitz::Array<double, 1>& x, 
				 const ArrayAd<double, 1>& y) const
{
  double resv = 0;
  Array<double, 1> resgrad(y.number_variable());
  resgrad = 0;
  for(int i = 1; i < x.rows(); ++i) {
    double t = x(i) - x(i - 1);
    resv +=  t * (y.value()(i) + y.value()(i - 1));
    resgrad += t * (y.jacobian()(i, Range::all()) + 
		    y.jacobian()(i - 1, Range::all()));
  }
  resv /= 2.0;
  resgrad /= 2.0;
  return AutoDerivative<double>(resv, resgrad);
}


// See base class for description.
ArrayAd<double, 1> IlsConvolution::apply_ils
(const blitz::Array<double, 1>& Hres_wn,
 const ArrayAd<double, 1>& Hres_rad,
 const std::vector<int>& Pixel_list) const
{
  firstIndex i1; secondIndex i2;
  if(Hres_wn.rows() != Hres_rad.rows())
    throw Exception("wave_number and radiance need to be the same size");
  ArrayAd<double, 1> disp_wn(disp->pixel_grid().data_ad());
  ArrayAd<double,1> res((int) Pixel_list.size(), 
			std::max(disp_wn.number_variable(), 
				 Hres_rad.number_variable()));
  // A few scratch variables, defined outside of loop so we don't keep
  // recreating them.
  Array<double, 1> normfact_grad(disp_wn.number_variable());
  ArrayAd<double, 1> conv(1, res.number_variable());
  ArrayAd<double, 1> response;
  for(int i = 0; i < res.rows(); ++i) {
    // Find the range of Hres_wn that lie within
    // wn_center +- ils_half_width
    AutoDerivative<double> wn_center = disp_wn(Pixel_list[i]);
    Array<double, 1>::const_iterator itmin = 
      std::lower_bound(Hres_wn.begin(), Hres_wn.end(),
		       wn_center.value() - ils_half_width_.value);
    Array<double, 1>::const_iterator itmax = 
      std::lower_bound(Hres_wn.begin(), Hres_wn.end(),
		       wn_center.value() + ils_half_width_.value);
    int jmin = (itmin == Hres_wn.end() ? Hres_wn.rows() - 1 : itmin.position()(0));
    int jmax = (itmax == Hres_wn.end() ? Hres_wn.rows() - 1 : itmax.position()(0));
    Range r(jmin, jmax);

    // Convolve with response

    // For speed, we just calculate dresponse/d[wn_center, coeff]. This along
    // with the chain rule is enough to calculate the Jacobian, and is
    // faster. 
    AutoDerivative<double> wn_center2(wn_center.value(), 0, 2);
    ils_func->ils(wn_center2, Hres_wn(r), response, true);

    // Note that this is currently hardcoded to assume at most one
    // coefficent, (like in IlsTable with 1 scale, of IlsGaussian with
    // 0 coefficients)
    Array<double, 1> coeff_grad;
    if(ils_func->used_flag_func().rows() > 1)
      throw Exception("IlsConvolution is currently hardwired to assume at most one coefficient");
    if(ils_func->used_flag_func().rows() == 1 && ils_func->used_flag_func()(0))
      coeff_grad.reference(ils_func->coeff_func().jacobian()(0,Range::all()));

    // Response may not be normalized, so calculate normalization
    // factor.
    AutoDerivative<double> normfact_dwn = integrate(Hres_wn(r), response);
    // Convert gradient from dwn_center to dState
    normfact_grad = normfact_dwn.gradient()(0) * wn_center.gradient();
    if(coeff_grad.rows() > 0)
      normfact_grad += normfact_dwn.gradient()(1) * coeff_grad;
    AutoDerivative<double> normfact(normfact_dwn.value(), normfact_grad);

    
    conv.resize(response.rows(), res.number_variable());
    conv.value() = response.value() * Hres_rad(r).value();

    if (!Hres_rad.is_constant() && !wn_center.is_constant()) {
      conv.jacobian()  = response.value()(i1) * Hres_rad(r).jacobian()(i1, i2) +
        response.jacobian()(Range::all(), 0)(i1) * wn_center.gradient()(i2) * 
        Hres_rad(r).value()(i1);
      if (coeff_grad.rows() > 0)
        conv.jacobian() += response.jacobian()(Range::all(), 1)(i1) *
          coeff_grad(i2) * Hres_rad(r).value()(i1);
    }

    else if (!wn_center.is_constant()) {
      conv.jacobian()  = response.jacobian()(Range::all(), 0)(i1) * 
        wn_center.gradient()(i2) * Hres_rad(r).value()(i1);
      if (coeff_grad.rows() > 0)
        conv.jacobian() += response.jacobian()(Range::all(), 1)(i1) *
          coeff_grad(i2) * Hres_rad(r).value()(i1);
    }

    else if (!Hres_rad.is_constant()) {
      conv.jacobian()  = response.value()(i1) * Hres_rad(r).jacobian()(i1, i2);
      if (coeff_grad.rows() > 0)
        conv.jacobian() += response.jacobian()(Range::all(), 1)(i1) *
          coeff_grad(i2) * Hres_rad(r).value()(i1);
    }

    res(i) = integrate(Hres_wn(r), conv) / normfact;
  }

  return res;
}

boost::shared_ptr<Ils> IlsConvolution::clone() const
{
  return boost::shared_ptr<Ils>(new IlsConvolution(disp->clone(),
			   ils_func, ils_half_width_));
}

void IlsConvolution::print(std::ostream& Os) const 
{ 
  Os << "IlsConvolution:\n";
  OstreamPad opad(Os, "  ");
  opad << *disp << "\n" << *ils_func << "\n"
       << "Half Width: " << ils_half_width_ << "\n";
  opad.strict_sync();
}
