#ifndef ILS_FTS_H
#define ILS_FTS_H
#include "ils.h"
#include <vector>
#include "level_1b_fts.h"
#include "dispersion_polynomial.h"

namespace FullPhysics {

/****************************************************************//**
  This does an ILS convolution for FTS data. This wraps some old
  Fortran code that models and sinc + box-car. We may replace this
  at some point.
*******************************************************************/

class IlsFts : public Ils, public Observer<Dispersion> {
public:
  IlsFts(const boost::shared_ptr<DispersionPolynomial>& Disp,
	 const blitz::Array<double, 2>& Dispersion_perturb,
	 const boost::shared_ptr<Level1bFts>& Level_1b,
	 int Spec_index,
	 const std::string& Band_name,
	 const std::string& Hdf_band_name,
	 const DoubleWithUnit&
	 Ils_half_width = DoubleWithUnit((2000 + 1) * 1e-2, units::inv_cm));
  virtual void notify_update(const StateVector& Sv) 
  { /* Nothing to do, we have everything done by attached Disp. */ }
  virtual void notify_add(StateVector& Sv) { Sv.add_observer(*disp); }
  virtual void notify_remove(StateVector& Sv) { Sv.remove_observer(*disp);}
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
  virtual std::string band_name() const
  { return band_name_; }
  virtual std::string hdf_band_name() const
  { return hdf_band_name_; }
  virtual boost::shared_ptr<Ils> clone() const;
  virtual SpectralDomain pixel_grid() const
  {return disp->pixel_grid(); }
  virtual DoubleWithUnit ils_half_width() const 
  {return ils_half_width_;}
  virtual void ils_half_width(const DoubleWithUnit& half_width)
  {ils_half_width_ = half_width;}
   virtual void print(std::ostream& Os) const;
private:
  std::string band_name_;
  std::string hdf_band_name_;
  boost::shared_ptr<DispersionPolynomial> disp;
  blitz::Array<double, 2> disp_perturb;
  boost::shared_ptr<Level1bFts> level_1b;
  int spec_index;
  DoubleWithUnit ils_half_width_;
};
}
#endif
