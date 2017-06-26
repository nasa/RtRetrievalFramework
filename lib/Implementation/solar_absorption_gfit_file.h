#ifndef SOLAR_ABSORPTION_GFIT_FILE_H
#define SOLAR_ABSORPTION_GFIT_FILE_H
#include "solar_absorption_spectrum.h"

namespace FullPhysics {

// Note that we have the various equations here in Latex. This is
// fairly unreadable as text comments, but if you view the doxgyen
// documentation created by this it gives a nicely formatted output.

/****************************************************************//**
  This class calculates the solar absorption spectrum.

  This particular implementation reads the GFIT format absorption line
  list file. This is a fixed record file format that gives the
  absorption line list data.

  From the Fortran code:
 
  Calculates the solar optical thickness spectrum (SOT) at any wavelengths.
 
   Taking the exponential of SOT produces the solar spectrum
   as it would be observed at infinite spectral resolution.
 
   All solar lines are assumed to have a shape of the form
 
           \f[SOT = s  \exp\left(-\frac{x^2}{\sqrt{d^4+x^2 y^2}}\right)\f]
   where
           <DL>
           <DT>s</DT><DD>is the line-center optical thickness (dimensionless)</DD>
           <DT>x</DT><DD> is the frequency from line center (cm-1)</DD>
           <DT>y</DT><DD> is the 1/e folding width (cm-1)</DD>
           <DT>d</DT><DD> is the Doppler width (cm-1)</DD>
           </DL>
 
   In the doppler limit, i.e. \f$ d^2 \gg x.y \f$
          \f[SOT = s  \exp\left(-{\frac{x}{d}}^2\right)\f]
 
   In the far line wing limit, i.e. \f$x y \gg d^2\f$,  
          \f[SOT = s  \exp\left(- \left|\frac{x}{y}\right|\right)\f]
 
   So near the line center, the lineshape is Doppler, but in
   the line wings it decays exponentially (if y>0).
 
   This choice of lineshape has no physical basis. It just seems
   to give a reasonable representation is nearly all cases.
   The only cases in which this lineshape does not give an
   adequate representation of the absorption are the extremely
   broad lines of light atmos such as H (atomic hydrogen) or Mg.
   However, by representing the H absorptions as superpositions
   of two lines, one narrow and the other broad, adequate results
   were obtained.
 
   Molecular absorptions (e.g. CO, OH, NH, CN) tend to have narrow,
   Doppler lineshapes because they are confined to a relatively
   narrow layer in the cooler, upper, part of the solar atmosphere.
   In the hotter depths they are dissociated.
 
   Atomic transitions, on the other hand, are formed over a much
   wider range of solar altitudes, and hence temperatures. This
   gives rise to line shapes whose wings decay in an approximately
   exponential manner with the distance from line center. The line
   shape of equation (1) does a reasonable job in both cases.
 
   This subroutine also makes allowances for the effect of the
   finite FOV of the observing instrument, which gives rise to:
 
   -# broadening of the solar lines due to the linear variation
      of the Doppler shift from solar rotation across the solar disk.
   -# deepening of the solar lines due to limb darkening.
      It assumes that an instrument which observes the entire solar
      disk will observe lines which are, on average, twice the strength
      of an instrument just observing the center of the disk.

   \todo Note that the Fraction_solar_diameter parameter is completely
   ignored, the underlying Fortran code ignores this value and sets
   the variable "sld" to 1.0 regardless of the Fraction_solar_diameter
   passed in.
*******************************************************************/

class SolarAbsorptionGfitFile : public SolarAbsorptionSpectrum {  
public:
  SolarAbsorptionGfitFile(const std::string& Line_list_file,
                          double Fraction_solar_diameter = 1.0);
  virtual ~SolarAbsorptionGfitFile() {}
  virtual void print(std::ostream& Os) const;
  virtual Spectrum solar_absorption_spectrum(
     const SpectralDomain& spec_domain) const;

  //-----------------------------------------------------------------------
  /// Fraction of solar diameter that we view.
  //-----------------------------------------------------------------------
  double fraction_solar_diameter() const {return fraction_solar_diameter_;}

private:
  std::string line_list_file_;
  double fraction_solar_diameter_;
};
}
#endif
