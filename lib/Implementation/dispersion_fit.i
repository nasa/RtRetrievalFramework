// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "dispersion_fit.h"
%}
%base_import(generic_object)
%import "level_1b.i"
%import "double_with_unit.i"
%import "array_with_unit.i"
%fp_shared_ptr(FullPhysics::DispersionFit);

namespace FullPhysics {
class DispersionFit : public GenericObject {
public:
  DispersionFit(const boost::shared_ptr<Level1b>& Level1b);
  blitz::Array<double, 2> fit(const blitz::Array<double, 2> disp_initial, 
			      const DoubleWithUnit& aband_solar_line_location,
			      const DoubleWithUnit& aband_solar_line_width,
			      const DoubleWithUnit& aband_search_width,
			      const DoubleWithUnit& aband_ils_offset,
			      const ArrayWithUnit<double, 1>& offset_scaling) const;
};
}
