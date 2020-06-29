#ifndef ILS_CONVOLUTION_H
#define ILS_CONVOLUTION_H
#include "ils.h"
#include "ils_function.h"
#include "dispersion.h"

namespace FullPhysics {
  class HdfFile;
/****************************************************************//**
  This is a ILS where we use a Dispersion object to determine the
  wavenumbers of each pixel, and convolve against a IlsFunction.
*******************************************************************/

class IlsConvolution : public Ils, public Observer<Dispersion>,
                       public Observer<IlsFunction> {
public:
//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------
  IlsConvolution(const boost::shared_ptr<Dispersion>& Disp,
		 const boost::shared_ptr<IlsFunction>& Ils_func,
		 const DoubleWithUnit&
		 Ils_half_width = DoubleWithUnit(20, units::inv_cm))
    : disp(Disp), ils_func(Ils_func), ils_half_width_(Ils_half_width) 
  {
    disp->add_observer(*this);
    ils_func->add_observer(*this);
  }

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------

  IlsConvolution(const boost::shared_ptr<Dispersion>& Disp,
		 const boost::shared_ptr<IlsFunction>& Ils_func,
		 double Ils_half_width)
    : disp(Disp), ils_func(Ils_func), 
      ils_half_width_(Ils_half_width, units::inv_cm)
  {
    disp->add_observer(*this);
    ils_func->add_observer(*this);
  }
  virtual ~IlsConvolution() {}
  virtual void notify_update(const StateVector& Sv) 
  { /* Nothing to do, we have everything done by attached Disp. */ }
  virtual void notify_add(StateVector& Sv) { Sv.add_observer(*disp);
                                             Sv.add_observer(*ils_func); }
  virtual void notify_remove(StateVector& Sv) { Sv.remove_observer(*disp);
                                                Sv.remove_observer(*ils_func);}
  virtual void notify_update(const Dispersion& D)
  {
    notify_update_do(*this);
  }
  virtual blitz::Array<double, 1> apply_ils
  (const blitz::Array<double, 1>& High_resolution_wave_number,
   const blitz::Array<double, 1>& High_resolution_radiance,
   const std::vector<int>& Pixel_list) const;
  virtual ArrayAd<double, 1> apply_ils
  (const blitz::Array<double, 1>& High_resolution_wave_number,
   const ArrayAd<double, 1>& High_resolution_radiance,
   const std::vector<int>& Pixel_list) const;
  virtual void print(std::ostream& Os) const;
  virtual std::string band_name() const {return ils_func->band_name(); }
  virtual std::string hdf_band_name() const 
  { return ils_func->hdf_band_name();}
  virtual SpectralDomain pixel_grid() const
  { return disp->pixel_grid(); }
  virtual DoubleWithUnit ils_half_width() const
  {return ils_half_width_;}
  virtual void ils_half_width(const DoubleWithUnit& half_width) 
  {ils_half_width_ = half_width;}
  virtual boost::shared_ptr<Ils> clone() const;

//-----------------------------------------------------------------------
/// Underlying IlsFunction
//-----------------------------------------------------------------------
  boost::shared_ptr<IlsFunction> ils_function() const {return ils_func; }

//-----------------------------------------------------------------------
/// Underlying dispersion.
//-----------------------------------------------------------------------
  boost::shared_ptr<Dispersion> dispersion() const {return disp; }

private:
  boost::shared_ptr<Dispersion> disp;
  boost::shared_ptr<IlsFunction> ils_func;
  DoubleWithUnit ils_half_width_;
  double integrate(const blitz::Array<double, 1>& x, 
		   const blitz::Array<double, 1>& y) const;
  AutoDerivative<double> integrate(const blitz::Array<double, 1>& x, 
				   const ArrayAd<double, 1>& y) const;
};
}
#endif
