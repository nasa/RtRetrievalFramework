#ifndef FORWARD_MODEL_H
#define FORWARD_MODEL_H
#include "printable.h"
#include "state_vector.h"
#include "auto_derivative.h"
#include "spectrum.h"
#include <blitz/array.h>
#include <boost/optional.hpp>

namespace FullPhysics {
/****************************************************************//**
  This is the Forward Model.

  Note that the forward model depends on the value of the StateVector.
  This is set separate from getting the radiance or jacobian
  values.
*******************************************************************/
class ForwardModel : public Printable<ForwardModel> {
public:
  virtual ~ForwardModel() {}
  virtual void print(std::ostream& Os) const {Os << "ForwardModel";}

//-----------------------------------------------------------------------
/// The grid that the forward model runs on may depend on the state
/// vector. This notifies the forward model that it should setup the
/// grid
//-----------------------------------------------------------------------

  virtual void setup_grid() = 0;

//-----------------------------------------------------------------------
/// The state vector associated with the forward model.
//-----------------------------------------------------------------------

  virtual boost::shared_ptr<StateVector> state_vector() const = 0;

//-----------------------------------------------------------------------
/// The number of spectral bands associated with forward model.
//-----------------------------------------------------------------------

  virtual int number_spectrometer() const = 0;

//-----------------------------------------------------------------------
/// The HDF field name to use for a particular band (e.g., "weak_co2")
//-----------------------------------------------------------------------

  virtual std::string hdf_band_name(int Spec_index) const = 0;

//-----------------------------------------------------------------------
/// Spectral domain for the given spectral band. Note that this may be
/// empty. 
//-----------------------------------------------------------------------

  virtual SpectralDomain spectral_domain(int Spec_index) const = 0;

//-----------------------------------------------------------------------
/// Spectrum for the given spectral band. Note that this may be empty.
///
/// \param Spec_index Band to give value for
/// \param Skip_jacobian If true, don't do the Jacobian
/// calculation. Often this is significantly faster to calculate.
/// \return The set of radiances, possibly empty.
//-----------------------------------------------------------------------
 
  virtual Spectrum radiance(int Spec_index, bool Skip_jacobian = false) 
    const = 0;

//-----------------------------------------------------------------------
/// Measured spectrum for the given spectral band. Note that this may be empty.
///
/// \param Spec_index Band to give value for
/// \return The set of radiances, possibly empty.
//-----------------------------------------------------------------------
 
  virtual Spectrum measured_radiance(int Spec_index) 
    const = 0;

//-----------------------------------------------------------------------
/// Type preference for spectral domain. This may seem an odd thing to
/// have a function for, but this is needed by ForwardModelOutput.
//-----------------------------------------------------------------------

  virtual SpectralDomain::TypePreference 
  spectral_domain_type_preference() const = 0;

  Spectrum radiance_all(bool Skip_jacobian = false) const;
  Spectrum measured_radiance_all() const;
  boost::optional<blitz::Range> pixel_range(int Spec_index) const;

//-----------------------------------------------------------------------
/// Description of the list of input files, can be using in print out
//-----------------------------------------------------------------------
  std::string input_file_description() const { return input_file_description_; }
  void input_file_description(const std::string& V) 
  { input_file_description_ = V; }
private:
  std::string input_file_description_;
};
}
#endif
