#include "ils_fts.h"
#include "linear_algebra.h"
#include "ostream_pad.h"
using namespace FullPhysics;
using namespace blitz;

boost::shared_ptr<Ils> ils_fts_create
(const boost::shared_ptr<Dispersion>& Disp,
 const blitz::Array<double, 2>& Disp_pert,
 const boost::shared_ptr<Level1b>& L1b,
 int Spec_index,
 const std::string& Band_name, const std::string& Hdf_band_name)
{
  const boost::shared_ptr<DispersionPolynomial> Disp_pol(boost::dynamic_pointer_cast<DispersionPolynomial>(Disp));
  const boost::shared_ptr<Level1bFts> L1b_fts(boost::dynamic_pointer_cast<Level1bFts>(L1b));
  return boost::shared_ptr<Ils>(new IlsFts(Disp_pol, Disp_pert,
                                           L1b_fts, Spec_index,
                                           Band_name, Hdf_band_name));
}

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(IlsFts, Ils)
.scope
[
 luabind::def("create", &ils_fts_create)
 ]
REGISTER_LUA_END()
#endif

extern "C" {

  void iof_model_instrument_fts(
    const double *ils_1,
    const double *ils_2,
    const double *ils_3,
    const double *dispersion_2,
    const double* disp_wn, 
    const int* disp_wn_size, 
    const double* wn, 
    const int* wn_size,
    const double* rad_calc, 
    const int *num_rad,
    const double *start,
    const double  *end,
    const double* rad_conv);
}

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------

IlsFts::IlsFts
(const boost::shared_ptr<DispersionPolynomial>& Disp,
 const blitz::Array<double, 2>& Dispersion_perturb,
 const boost::shared_ptr<Level1bFts>& Level_1b,
 int Spec_index,
 const std::string& Band_name,
 const std::string& Hdf_band_name,
 const DoubleWithUnit& Ils_half_width)
  : band_name_(Band_name), hdf_band_name_(Hdf_band_name),
    disp(Disp), level_1b(Level_1b), 
    spec_index(Spec_index), ils_half_width_(Ils_half_width)
{ 
  disp_perturb.reference(Dispersion_perturb.copy());
  disp->add_observer(*this); 
}

// See base class for description
blitz::Array<double, 1> IlsFts::apply_ils
(const blitz::Array<double, 1>& High_resolution_wave_number,
 const blitz::Array<double, 1>& High_resolution_radiance,
 const std::vector<int>& Pixel_list) const
{
  if(High_resolution_wave_number.rows() != High_resolution_radiance.rows())
    throw Exception("High_resolution_wave_number and High_resolution_radiance need to be the same size");
 
  // Apply convolution

  // wave_number might not be contiguous. Since the Fortran code
  // requires that it is, we make a copy if needed.
  Array<double, 1> wave_number_contiguous = 
    to_fortran_const(High_resolution_wave_number);
  int number_wavelength = wave_number_contiguous.rows();

  Array<double, 1> radiance_contiguous = 
    to_fortran_const(High_resolution_radiance);
  int nrad=radiance_contiguous.rows();

  Array<double, 1> disp_wn(pixel_grid().wavenumber());
  int dwn_size = disp_wn.rows();

  // Result array
  Array<double,1> res_full(dwn_size);

  double spacing = disp_wn(1) - disp_wn(0);
  double opd = level_1b->run_log(spec_index).optical_path_difference;
  double fov = level_1b->run_log(spec_index).internal_fov;
  double ang_misalignment = level_1b->run_log(spec_index).angular_misalignment;
  
  double start = level_1b->frequency_start(spec_index);
  double end   = level_1b->frequency_end(spec_index);

  iof_model_instrument_fts(&opd,&fov,&ang_misalignment,&spacing,
                           disp_wn.data(),&dwn_size,
                           wave_number_contiguous.data(),
                           &number_wavelength,radiance_contiguous.data(),
                           &nrad,&start,&end,res_full.dataFirst());

  Array<double,1> res(Pixel_list.size());
  for (int i=0; i< (int) Pixel_list.size(); ++i) {
     res(i) = res_full(Pixel_list[i]);
  }
  
  return res;
}

// See base class for description
ArrayAd<double, 1> IlsFts::apply_ils
(const blitz::Array<double, 1>& High_resolution_wave_number,
 const ArrayAd<double, 1>& High_resolution_radiance,
 const std::vector<int>& Pixel_list) const
{
  ArrayAd<double, 1> res((int) Pixel_list.size(), 
                         High_resolution_radiance.number_variable());
  res.value() = apply_ils(High_resolution_wave_number, 
                          High_resolution_radiance.value(), 
                          Pixel_list);
  // Apply instrument model to jacobians
  for(int i = 0; i < res.number_variable(); ++i)
    res.jacobian()(Range::all(), i) = 
      apply_ils(High_resolution_wave_number, 
                High_resolution_radiance.jacobian()(Range::all(),i), 
                Pixel_list);

  // This will go away shortly.
  // svi = -1 in FM only mode
  int svi = disp->state_vector_start_index();
  if (svi > 0) {
    Array<double, 1> coeff(disp->coefficient().value());
    for(int i = 0; i < disp->used_flag_value().rows(); ++i)
      if(disp->used_flag_value()(i)) {
        double old_disp = coeff(i);
        coeff(i) += disp_perturb(spec_index, i);
        res.jacobian()(Range::all(), svi) +=
            (apply_ils(High_resolution_wave_number, 
                       High_resolution_radiance.value(), 
                       Pixel_list) -
             res.value()) / disp_perturb(spec_index, i);
        coeff(i) = old_disp;
        ++svi;
      }
  }
  return res;
}

// See base class for description
boost::shared_ptr<Ils> IlsFts::clone() const
{
  return boost::shared_ptr<Ils>
    (new IlsFts
     (boost::dynamic_pointer_cast<DispersionPolynomial>(disp->clone()),
      disp_perturb, level_1b, spec_index, band_name_, hdf_band_name_, 
      ils_half_width_));
}

// See base class for description
void IlsFts::print(std::ostream& Os) const 
{ 
  Os << "IlsFts:\n";
  Os << "Band: " << band_name_ << "\n";
  OstreamPad opad(Os, "  ");
  opad << *disp << "\n";
  opad << *level_1b << "\n";
  opad.strict_sync();
}


