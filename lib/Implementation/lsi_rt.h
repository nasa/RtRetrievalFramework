#ifndef LSI_RT_H
#define LSI_RT_H
#include "radiative_transfer_single_wn.h"
#include "hdf_file.h"
#include "log_interpolate.h"

namespace FullPhysics {
/****************************************************************//**
  This does a Low Stream Interpolator correction to another
  RadiativeTransfer object. There is a paper in the doc directory "LSI
  Paper.pdf" which describes this algorithm. Note that the paper
  describes an improved version, the code here is for an older
  implementation (we will update the papers version eventually).

  There is a configuration file that gives the optical depth
  boundaries to use in the LSI binning. This can either be read from
  an HDF file (the preferred way), or for backwards compatibility from
  an ASCII file. 
*******************************************************************/
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
  virtual ~LsiRt() {}
  virtual int number_stokes() const
  { return high_stream_rt->number_stokes(); }
  virtual blitz::Array<double, 2> stokes(const SpectralDomain& Spec_domain,
					 int Spec_index) const;
  virtual ArrayAd<double, 2> stokes_and_jacobian
  (const SpectralDomain& Spec_domain, int Spec_index) const;
  virtual ArrayAd<double, 2> correction_only
  (const SpectralDomain& Spec_domain, int Spec_index) const;
  virtual void print(std::ostream& Os, bool Short_form = false) const;
  boost::shared_ptr<RadiativeTransfer>
  low_stream_radiative_transfer() const { return low_stream_rt; }
  boost::shared_ptr<RadiativeTransfer>
  high_stream_radiative_transfer() const { return high_stream_rt; }
private:
  void calc_correction(const SpectralDomain& Spec_domain,
		       int Spec_index, bool Calc_jacobian, 
		       bool Skip_stokes_calc) const;
  boost::shared_ptr<RadiativeTransferSingleWn> low_stream_rt, high_stream_rt;
  std::vector<std::vector<double> > optical_depth_boundary;
  std::vector<std::string> main_gas;
  double wv_threshold;
  mutable blitz::Array<double, 2> stokes_only;
  mutable ArrayAd<double, 2> stokes_and_jac;
  mutable ArrayAd<double, 1> gas_opd;
  mutable blitz::Array<int, 1> wv_index;
  typedef LinearInterpolate<AutoDerivative<double>, 
       blitz::Array<AutoDerivative<double>, 1> > linear_interp_type;
  typedef LogLinearInterpolate<AutoDerivative<double>, 
       blitz::Array<AutoDerivative<double>, 1> > log_linear_interp_type;
  mutable std::vector<linear_interp_type> linear_interp;
  mutable std::vector<log_linear_interp_type> log_linear_interp;
};
}
#endif
