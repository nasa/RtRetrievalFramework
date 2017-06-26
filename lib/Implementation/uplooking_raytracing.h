#ifndef UPLOOKING_RAYTRACING_H
#define UPLOOKING_RAYTRACING_H

#include <vector>
#include <blitz/array.h>

#include "level_1b_fts.h"
#include "atmosphere_oco.h"

namespace FullPhysics {

// Use this struct as a namespace for these functions making it easy
// to set up for use by Lua
struct UplookingRaytracing : public Printable<UplookingRaytracing> {

  static ArrayWithUnit<double, 1> calculate_apparent_surface_sza
  (const boost::shared_ptr<Level1bFts>& l1b,
   const boost::shared_ptr<AtmosphereOco>& atm);

  static void tlpath(int nlev,double * z, double * t, double * p, double asza_in,
		     double fovr, double roc, double zobs, double wavtkr,
		     double wavmic, double * zmin_out, double * bend_out,
		     double * sp, int * ifail);
  virtual void print(std::ostream& Os) const {Os << "UplookingRaytracing";}
};
}

#endif
