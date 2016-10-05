#ifndef DISPERSION_FIT_H
#define DISPERSION_FIT_H
#include "level_1b.h"

namespace FullPhysics {
/****************************************************************//**
 Given a single frame of data, estimate the spectral shift in band 1P.
 This is accomplished by using the strong solar line at 12985.163 wavenumbers.
 This is in the a-band continuum and is really the only local feature there.

 Note: This routine fits for a global shift, which INCLUDES the instrument
 doppler shift.  If this is not desired, the user must subtract off
 the instrument doppler shift.
   
 Original Author: Chris Odell
 Converted from IDL
*******************************************************************/
  class DispersionFit : public Printable<DispersionFit> {
public:

    DispersionFit(const boost::shared_ptr<Level1b>& Level1b);
    blitz::Array<double, 2> fit(const blitz::Array<double, 2> disp_initial, 
				const DoubleWithUnit& aband_solar_line_location,
				const DoubleWithUnit& aband_solar_line_width,
				const DoubleWithUnit& aband_search_width,
				const DoubleWithUnit& aband_ils_offset,
				const ArrayWithUnit<double, 1>& offset_scaling) const;
    
    /// Return computed shift in each band
    blitz::Array<double, 1> shift() const { return shift_; };

    void print(std::ostream& Os) const;

private:

  boost::shared_ptr<Level1b> l1b;
  mutable blitz::Array<double, 1> shift_;

};
}
#endif


