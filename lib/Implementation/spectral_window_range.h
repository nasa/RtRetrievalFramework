#ifndef SPECTRAL_WINDOW_RANGE_H
#define SPECTRAL_WINDOW_RANGE_H
#include "heritage_file.h"
#include "spectral_window.h"
#include "array_with_unit.h"

namespace FullPhysics {
  class Dispersion;

/****************************************************************//**
  This is an implementation of a SpectralWindow that covers a fixed
  window. The window can be made up of multiple microwindow ranges if
  desired. It can also include a bad pixel mask.

  Note that there are a few closely related classes, with similar 
  sounding names. See \ref spectrum_doxygen for a description of each
  of these.
*******************************************************************/

class SpectralWindowRange : public SpectralWindow {
public:
  SpectralWindowRange(const ArrayWithUnit<double, 3>& Microwindow_ranges);

  template<class T>
  SpectralWindowRange(const ArrayWithUnit<double, 3>& Microwindow_ranges,
		      const blitz::Array<T, 2>& Bad_sample_mask);
  virtual ~SpectralWindowRange() {}

//-----------------------------------------------------------------------
/// The Dispersion is optional. If supplied, we can convert from
/// sample_index to wavelength/wavenumber.
//-----------------------------------------------------------------------
  
  const std::vector<boost::shared_ptr<Dispersion> >& dispersion() const 
  {return disp_;}
  void dispersion(const std::vector<boost::shared_ptr<Dispersion> >& D)
  { disp_ = D; }
  
  virtual int number_spectrometer() const {return range_.value.rows();}
  virtual SpectralBound spectral_bound() const;
  virtual std::vector<int> grid_indexes(const SpectralDomain& Grid, 
					int Spec_index) const;
  void print(std::ostream& Os) const;
  const ArrayWithUnit<double, 3>& range_array() const {return range_;}
  void range_array(const ArrayWithUnit<double, 3>& Ran) 
  { range_ = Ran;}
  const blitz::Array<bool, 2>& bad_sample_mask() const 
  { return bad_sample_mask_; }
  void bad_sample_mask(const blitz::Array<bool, 2>& M) 
  { bad_sample_mask_.reference(M.copy()); }
private:
  ArrayWithUnit<double, 3> range_;
  // Mask of bad samples, True for a bad sample, False for a good one
  blitz::Array<bool, 2> bad_sample_mask_;
  std::vector<boost::shared_ptr<Dispersion> > disp_;
};
}
#endif
