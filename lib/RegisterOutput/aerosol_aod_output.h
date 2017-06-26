#ifndef AEROSOL_AOD_OUTPUT_H
#define AEROSOL_AOD_OUTPUT_H
#include "register_output_base.h"
#include "aerosol_optical.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the Aerosol class that should be
  written as output.

  See the discussion in RegisterOutputBase why this isn't just part of
  the Aerosol class.
*******************************************************************/
class AerosolAodOutput : public RegisterOutputBase {
public:
//-----------------------------------------------------------------------
/// Constructor. We either name the output fields by the aerosol name
/// (e.g., aerosol_ice_aod), or we name it by the index number
/// (e.g., aerosol_1_aod). This is controlled by the
/// Number_instead_of_name field.
//-----------------------------------------------------------------------
  AerosolAodOutput(const boost::shared_ptr<Aerosol>& A, 
		   bool Number_instead_of_name = false);
  virtual ~AerosolAodOutput() {}
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  static const double low_boundary;
  static const double high_boundary;
private:
  boost::shared_ptr<AerosolOptical> a;
  int number_instead_of_name;
};
}
#endif
