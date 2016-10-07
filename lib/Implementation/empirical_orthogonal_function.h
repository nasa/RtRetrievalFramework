#ifndef EMPIRICAL_ORTHOGONAL_FUNCTION_H
#define EMPIRICAL_ORTHOGONAL_FUNCTION_H
#include "dispersion.h"
#include "instrument_correction.h"
#include "hdf_file.h"
#include "sub_state_vector_array.h"
#include "array_with_unit.h"

namespace FullPhysics {
/****************************************************************//**
  This class applies a empirical orthogonal function (EOF) correction
  to instrument data. We use a supplied waveform, multiple by a single
  scale factor given by the state vector, and add this to the radiance
  calculated in InstrumentIls.

  Note that other than what we call this and there various metadata
  fields, this is the same thing as the ZeroOffsetWaveform.
*******************************************************************/

class EmpiricalOrthogonalFunction : 
public SubStateVectorArray<InstrumentCorrection> {
public:
//-----------------------------------------------------------------------
/// Constructor.
///
/// \param Coeff - Initial value of scale factor
/// \param Used_flag - If true, we update scale factor by values in
///    StateVector. If false, we hold this fixed and just used the
///    initial value.
/// \param Eof_waveform - Offset to add for zero level. This
///    is indexed by the instrument pixel.
/// \param Order - Order of the eigenvector (e.g., first order
///    correction, second order correction, etc.)
/// \param Band_name - Name of band
/// \param Hdf_group - HDF group name, if this was read from an HDF
///    file. This is only used in the diagnostic print out of this
///    object, we aren't actually reading from HDF in this particular
///    constructor. If the Eof_waveform does not come from and HDF
///    file, you can just leave this as the default "N/A" value.
/// \param Sounding_number - The sounding number
/// \param Eof_depend_on_sounding_number - True if the EOF depends
///    on the sounding number, false otherwise      
//-----------------------------------------------------------------------
  
  EmpiricalOrthogonalFunction(double Coeff, 
			      bool Used_flag,
			      const ArrayWithUnit<double, 1>& Eof_waveform,
			      int Order,
			      const std::string& Band_name,
			      const std::string& Hdf_group = "N/A",
			      int Sounding_number = 0,
			      bool Eof_depend_on_sounding_number = false)
    : SubStateVectorArray<InstrumentCorrection>(Coeff, Used_flag),
      band_name(Band_name),
      hdf_group(Hdf_group),
      order_(Order),
      sounding_number_(Sounding_number),
      eof_depend_on_sounding_number_(Eof_depend_on_sounding_number),
      eof_(Eof_waveform)
  { }

  EmpiricalOrthogonalFunction(double Coeff, 
			      bool Used_flag,
			      const Dispersion& Disp,
			      const HdfFile& Hdf_static_input,
			      int Spec_index,
			      int Sounding_number,
			      int Order,
			      const std::string& Band_name,
			      const std::string& Hdf_group = 
			      "Instrument/EmpiricalOrthogonalFunction");
  EmpiricalOrthogonalFunction(double Coeff, 
			      bool Used_flag,
			      const HdfFile& Hdf_static_input,
			      int Spec_index,
			      int Sounding_number,
			      int Order,
			      const std::string& Band_name,
			      const std::string& Hdf_group = 
			      "Instrument/EmpiricalOrthogonalFunction");
  EmpiricalOrthogonalFunction(double Coeff, 
			      bool Used_flag,
			      const HdfFile& Hdf_static_input,
			      const ArrayWithUnit<double, 1>& Uncertainty,
			      int Spec_index,
			      int Sounding_number,
			      int Order,
			      const std::string& Band_name,
			      const std::string& Hdf_group = 
			      "Instrument/EmpiricalOrthogonalFunction",
			      double Scale_to_stddev = 1e19);
  virtual ~EmpiricalOrthogonalFunction() {}
  virtual std::string state_vector_name_i(int i) const
  { return "EOF order " + boost::lexical_cast<std::string>(order_) +
      " scale factor " + band_name; }
  virtual boost::shared_ptr<InstrumentCorrection> clone() const;
  virtual void apply_correction
  (const SpectralDomain& Pixel_grid,
   const std::vector<int>& Pixel_list,
   SpectralRange& Radiance) const;
  virtual void print(std::ostream& Os) const;

//-----------------------------------------------------------------------
/// Current value of empirical orthogonal function, for each pixel number.
//-----------------------------------------------------------------------
  ArrayWithUnit<double, 1> eof() const
  {
    return ArrayWithUnit<double, 1>
      (blitz::Array<double, 1>(coeff(0).value() * eof_.value), 
       eof_.units);
  }

//-----------------------------------------------------------------------
/// Order of the empirical orthogonal function (e.g., first order,
/// second order, etc.)
//-----------------------------------------------------------------------

  int order() const {return order_;}
//-----------------------------------------------------------------------
/// Scale. This is just coeff(0), but we wrap this for use
/// by EofOutput
//-----------------------------------------------------------------------

  double scale() const { return coeff.value()(0); }

//-----------------------------------------------------------------------
/// Scale uncertainty. This is just sqrt(Cov(0,0)), but we
/// wrap this for use  by EofOutput
//-----------------------------------------------------------------------

  double scale_uncertainty() const 
  { 
    if(sv_cov_sub.rows() < 1)
      return 0;
    double t = sv_cov_sub(0,0);
    return (t < 0 ? 0 : sqrt(t)); 
  }

private:
  std::string band_name;
  std::string hdf_group;
  int order_;
  int sounding_number_;
  bool eof_depend_on_sounding_number_;
  bool eof_scale_uncertainty_;
  double scale_to_stddev_;
  ArrayWithUnit<double, 1> eof_;
};
}
#endif
