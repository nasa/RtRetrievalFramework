#ifndef OCO_SIM_APRIORI_H
#define OCO_SIM_APRIORI_H
#include "pressure.h"
#include "hdf_file.h"
#include "oco_sounding_id.h"
#include "linear_interpolate.h"

namespace FullPhysics {
/****************************************************************//**
  This class is used to calculate a CO2 apriori using the scene
  file from a OCO simulator run.
*******************************************************************/

class OcoSimApriori : public Printable<OcoSimApriori> {
public:
  OcoSimApriori(const std::string& Oco_sim_scene,
		const HdfSoundingId& Sid);
  virtual ~OcoSimApriori() {}
  double co2_vmr(double P) const { return interp(P); }
  blitz::Array<double, 1> co2_vmr_grid(const Pressure& P) const;
  void print(std::ostream& Os) const { Os << "OcoSimApriori"; }
private:
  LinearInterpolate<double, double> interp;
};
}

#endif

