#ifndef CO2_PROFILE_PRIOR_H
#include "pressure.h"
#include "hdf_file.h"
#include "oco_met_file.h"

namespace FullPhysics {
/****************************************************************//**
  This class is used to create the CO2 apriori from a profile file
  (L2CPr).
*******************************************************************/

class CO2ProfilePrior: public Printable<CO2ProfilePrior> {
public:
  CO2ProfilePrior(const OcoMetFile& Met_file,
		  const HdfFile& Profile_file);
  virtual ~CO2ProfilePrior() {}
  blitz::Array<double, 1> apriori_vmr(const Pressure& pressure) const;
  virtual void print(std::ostream& Os) const { Os << "CO2ProfilePrior";}
private:
  blitz::Array<double, 1> model_press;
  blitz::Array<double, 1> co2_vmr;
};
}
#endif
