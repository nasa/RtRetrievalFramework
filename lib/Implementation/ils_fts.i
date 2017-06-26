// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "ils_fts.h"
%}
%base_import(ils)
%base_import(observer)
%base_import(dispersion)
%import "dispersion_polynomial.i"
%import "level_1b_fts.i"
%import "state_vector.i"
%import "array_ad.i"
%fp_shared_ptr(FullPhysics::IlsFts);

namespace FullPhysics {
class IlsFts : public Ils, public Observer<Dispersion> {
public:
  IlsFts(const boost::shared_ptr<DispersionPolynomial>& Disp,
	 const blitz::Array<double, 2>& Dispersion_perturb,
	 const boost::shared_ptr<Level1bFts>& Level_1b,
	 int Spec_index,
	 const std::string& Band_name,
	 const std::string& Hdf_band_name);
  virtual void notify_update(const StateVector& Sv);
  virtual void notify_add(StateVector& Sv);
  virtual void notify_remove(StateVector& Sv);
  virtual void notify_update(const Dispersion& D);
  virtual blitz::Array<double, 1> apply_ils
  (const blitz::Array<double, 1>& High_resolution_wave_number,
   const blitz::Array<double, 1>& High_resolution_radiance,
   const std::vector<int>& Pixel_list) const;
  virtual ArrayAd<double, 1> apply_ils
  (const blitz::Array<double, 1>& High_resolution_wave_number,
   const ArrayAd<double, 1>& High_resolution_radiance,
   const std::vector<int>& Pixel_list) const;
  virtual boost::shared_ptr<Ils> clone() const;
};
}
