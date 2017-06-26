// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "lsi_rt.h"
#include "sub_state_vector_array.h"
#include "pressure.h"
%}
%base_import(radiative_transfer_fixed_stokes_coefficient)
%import "radiative_transfer_single_wn.i"
%import "spectral_domain.i"
%import "hdf_file.i"

%fp_shared_ptr(FullPhysics::LsiRt);

namespace FullPhysics {
class LsiRt : public RadiativeTransferFixedStokesCoefficient {
public:
  LsiRt(const boost::shared_ptr<RadiativeTransferSingleWn>& Low_stream_rt,
	const boost::shared_ptr<RadiativeTransferSingleWn>& High_stream_rt,
	const std::string& Lsi_fname, 
	double Water_vapor_fraction_threshold = 0.8);
  LsiRt(const boost::shared_ptr<RadiativeTransferSingleWn>& Low_stream_rt,
	const boost::shared_ptr<RadiativeTransferSingleWn>& High_stream_rt,
	const HdfFile& Config_file,
	const std::string& Lsi_group = "LSI",
	double Water_vapor_fraction_threshold = 0.8);
  %python_attribute(number_stokes, virtual int)
  virtual blitz::Array<double, 2> stokes(const SpectralDomain& Spec_domain,
					 int Spec_index) const;
  virtual ArrayAd<double, 2> stokes_and_jacobian
  (const SpectralDomain& Spec_domain, int Spec_index) const;
  virtual ArrayAd<double, 2> correction_only
  (const SpectralDomain& Spec_domain, int Spec_index) const;
  %python_attribute(low_stream_radiative_transfer, boost::shared_ptr<RadiativeTransfer>)
  %python_attribute(high_stream_radiative_transfer, boost::shared_ptr<RadiativeTransfer>)
};
}

