#include "forward_model_spectral_grid.h"
#include "spectral_domain.h"
#include "instrument.h"
#include "spectral_window.h"
#include "spectrum_sampling.h"
#include "linear_interpolate.h"
using namespace FullPhysics;
using namespace blitz;

ForwardModelSpectralGrid::ForwardModelSpectralGrid
(const Instrument& Inst,
 const SpectralWindow& Spectral_window,
 const SpectrumSampling& Spectrum_sampling)
{
    for(int i = 0; i < Spectral_window.number_spectrometer(); ++i) {
        std::vector<int> spec_plist = Spectral_window.grid_indexes(Inst.pixel_spectral_domain(i), i);
        if (spec_plist.size() > 0) {
            plist.push_back(spec_plist);
            lgrid.push_back(Spectral_window.apply(Inst.pixel_spectral_domain(i), i));
            hgrid.push_back(Spectrum_sampling.spectral_domain(i, lgrid[i], 
                         Inst.ils_half_width(i)));
            hgrid_inter.push_back(Spectrum_sampling.
                    spectral_domain_interpolated(i, lgrid[i],
                    Inst.ils_half_width(i)));
            need_interpolation.push_back(Spectrum_sampling.need_interpolation(i));
        } else {
            plist.push_back(std::vector<int> ());
            SpectralDomain empty_domain = SpectralDomain(blitz::Array<double, 1>(0));
            lgrid.push_back(empty_domain);
            hgrid.push_back(empty_domain);
            hgrid_inter.push_back(empty_domain);
            need_interpolation.push_back(false);
        }
    }
}

//-----------------------------------------------------------------------
/// Interpolate a spectrum to the high_resolution_interpolated_grid()
/// sampling. 
//-----------------------------------------------------------------------
Spectrum ForwardModelSpectralGrid::interpolate_spectrum
(const Spectrum& Spec_in, int Spec_index) const
{ 
  range_check(Spec_index, 0, number_spectrometer());
  if(!need_interpolation[Spec_index])
    return Spec_in;
  Range ra(Range::all());
  ArrayAd<double, 1> res(hgrid_inter[Spec_index].data().rows(),
			 Spec_in.spectral_range().data_ad().number_variable());
  Array<double, 1> spec_in_dom = Spec_in.spectral_domain().data();
  Array<double, 1> ispec_dom = 
    hgrid_inter[Spec_index].convert_wave(Spec_in.spectral_domain().units());
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
  return Spectrum(hgrid_inter[Spec_index].clone(), 
		  SpectralRange(res, Spec_in.spectral_range().units()));
}
