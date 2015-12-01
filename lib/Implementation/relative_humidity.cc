#include "relative_humidity.h"
#include "ostream_pad.h"
using namespace FullPhysics;

//-----------------------------------------------------------------------
/// Print to a string.
//-----------------------------------------------------------------------
void RelativeHumidity::print(std::ostream& Os)
{
  OstreamPad opad(Os, "    ");
  Os << "RelativeHumidity:\n";
  Os << "  Absorber:\n";
  opad << *absorber << "\n";
  opad.strict_sync();
  Os << "  Temperature:\n";
  opad << *temp << "\n";
  opad.strict_sync();
  Os << "  Pressure:\n";
  opad << *press << "\n";
  opad.strict_sync();
}

//-----------------------------------------------------------------------
/// Calculate relative humidity.
//-----------------------------------------------------------------------

ArrayAd<double, 1> RelativeHumidity::relative_humidity_grid() const
{
}
