#include "forward_model_spectral_grid.h"
#include "linear_interpolate.h"
using namespace FullPhysics;
using namespace blitz;

//-----------------------------------------------------------------------
/// The low resolution grid.
//-----------------------------------------------------------------------

  const SpectralDomain ForwardModelSpectralGrid::low_resolution_grid(int Spec_index) const
  { 
    range_check(Spec_index, 0, number_spectrometer());

    std::vector<int> plist = pixel_list(Spec_index);
    if (plist.size() > 0) {
        return spectral_window->apply(inst->pixel_spectral_domain(Spec_index), Spec_index);
    } else {
        return SpectralDomain(blitz::Array<double, 1>(0));
    }
  }

//-----------------------------------------------------------------------
/// The high resolution grid, possibly nonuniform
//-----------------------------------------------------------------------

  const SpectralDomain ForwardModelSpectralGrid::high_resolution_grid(int Spec_index) const
  { 
    range_check(Spec_index, 0, number_spectrometer());

    std::vector<int> plist = pixel_list(Spec_index);
    if (plist.size() > 0) {
        return spectrum_sampling->spectral_domain(Spec_index, low_resolution_grid(Spec_index), inst->ils_half_width(Spec_index));
    } else {
        return SpectralDomain(blitz::Array<double, 1>(0));
    }
  }

//-----------------------------------------------------------------------
/// The high resolution grid, interpolated to be uniform.
//-----------------------------------------------------------------------

  const SpectralDomain ForwardModelSpectralGrid::high_resolution_interpolated_grid(int Spec_index) const
  { 
    range_check(Spec_index, 0, number_spectrometer());

    std::vector<int> plist = pixel_list(Spec_index);
    if (plist.size() > 0) {
        return spectrum_sampling->spectral_domain_interpolated(Spec_index, low_resolution_grid(Spec_index), inst->ils_half_width(Spec_index));
    } else {
        return SpectralDomain(blitz::Array<double, 1>(0));
    }
  }

//-----------------------------------------------------------------------
/// Pixel indexes to use for low resolution grid.
//-----------------------------------------------------------------------
  const std::vector<int> ForwardModelSpectralGrid::pixel_list(int Spec_index) const
  {
    range_check(Spec_index, 0, number_spectrometer());
    return spectral_window->grid_indexes(inst->pixel_spectral_domain(Spec_index), Spec_index);
  }

//-----------------------------------------------------------------------
/// Interpolate a spectrum to the high_resolution_interpolated_grid()
/// sampling. 
//-----------------------------------------------------------------------
Spectrum ForwardModelSpectralGrid::interpolate_spectrum
(const Spectrum& Spec_in, int Spec_index) const
{ 
  range_check(Spec_index, 0, number_spectrometer());
  if(!spectrum_sampling->need_interpolation(Spec_index))
    return Spec_in;
  Range ra(Range::all());
  SpectralDomain hgrid_inter = high_resolution_interpolated_grid(Spec_index);
  ArrayAd<double, 1> res(hgrid_inter.data().rows(),
             Spec_in.spectral_range().data_ad().number_variable());
  Array<double, 1> spec_in_dom = Spec_in.spectral_domain().data();
  Array<double, 1> ispec_dom = 
    hgrid_inter.convert_wave(Spec_in.spectral_domain().units());
  LinearInterpolate<double, double> 
    interp(spec_in_dom.begin(),
       spec_in_dom.end(),
       Spec_in.spectral_range().data().begin(),
       LinearInterpolate<double, double>::OUT_OF_RANGE_ERROR);
  for(int i = 0; i< ispec_dom.rows(); ++i)
    res.value()(i) = interp(ispec_dom(i));
  if(res.number_variable() > 0) {
    std::vector<Array<double, 1> > yvec;
    for(int i = 0; i< spec_in_dom.rows(); ++i)
      yvec.push_back(Spec_in.spectral_range().data_ad().jacobian()(i, ra));
    LinearInterpolate<double, Array<double, 1> > 
      jacinterp(spec_in_dom.begin(), spec_in_dom.end(),
        yvec.begin(),
    LinearInterpolate<double, Array<double, 1> >::OUT_OF_RANGE_ERROR);
    for(int i = 0; i< ispec_dom.rows(); ++i)
      res.jacobian()(i,ra) = jacinterp(ispec_dom(i));
  }
  return Spectrum(hgrid_inter.clone(), 
          SpectralRange(res, Spec_in.spectral_range().units()));
}
