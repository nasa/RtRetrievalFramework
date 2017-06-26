#include "rayleigh.h"

using namespace FullPhysics;
using namespace blitz;
extern "C" {
  void rayleigh_wrap(const double* wn,
		     const double* depolar_fact,
		     const double* a,
		     const double* b,
		     double* rayleigh_cross_section);
  /* depolar_fact depolarization factor for air is 0.02790 (young, 1980) */
  /* a and b are wavelength dependence coefficients for the refractive index
       (allen, 1964) (note: wavelengths must be in microns)
       For Earth, a = 2.871e-04, b = 5.67e-03 */
}

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------
Rayleigh::Rayleigh(const boost::shared_ptr<Pressure>& Pres, 
		   const std::vector<boost::shared_ptr<Altitude> >& Alt,
		   const Constant& C) 
  : pres(Pres), alt(Alt), cache_is_stale(true),
    a(C.rayleigh_a().value), b(C.rayleigh_b().value), 
    depolar_fact(C.rayleigh_depolarization_factor()),
    molar_weight_dry_air(C.molar_weight_dry_air().convert(Unit("g / mol")).value)
{ 
  pres->add_observer(*this);
  BOOST_FOREACH(boost::shared_ptr<Altitude>& a, alt)
    a->add_observer(*this);
}


//-----------------------------------------------------------------------
/// Calculate the rayleigh cross section for the given
/// wavenumber/wavelength. 
//-----------------------------------------------------------------------

DoubleWithUnit Rayleigh::cross_section(const DoubleWithUnit& W,
				       const Constant& C)
{
  double rayleigh_cross_section;
  double wn = W.convert_wave(units::inv_cm).value;
  double a = C.rayleigh_a().value;
  double b = C.rayleigh_b().value;
  double depolar_fact = C.rayleigh_depolarization_factor();
  rayleigh_wrap(&wn, &depolar_fact, &a, &b, &rayleigh_cross_section);
  return DoubleWithUnit(rayleigh_cross_section, Unit("m^2"));
}

//-----------------------------------------------------------------------
/// This gives the optical depth for each layer, for the given wave
/// number. Note this only includes the Rayleigh portion of this,
/// Atmosphere class combines this with Absorbers and Aerosol
/// scattering.
///
/// This has size of pres->number_active_layer().
///
/// \todo Determine what the mysterious a0 number is, and if it needs
/// to be changed.
//-----------------------------------------------------------------------

ArrayAd<double, 1> 
Rayleigh::optical_depth_each_layer(double wn, int spec_index) const
{
  range_check(spec_index, 0, (int) alt.size());
  fill_cache();
  double rayleigh_cross_section;

  rayleigh_wrap(&wn, &depolar_fact, &a, &b, &rayleigh_cross_section);
  ArrayAd<double, 1> res(part_independent_wn(spec_index, Range::all()).copy());
  res.value() *= rayleigh_cross_section;
  res.jacobian() *= rayleigh_cross_section;
  return res;
}

//-----------------------------------------------------------------------
/// Fill in cache, if needed.
//-----------------------------------------------------------------------
void Rayleigh::fill_cache() const
{
  if(!cache_is_stale)
    return;

  // A comment in the old fortran code indicates this is Avogadro's
  // number. However, this actually is not (Avogadro's number is
  // 6.02214179e23 /mol.  I'm not actually sure what this number is. 
  // You can look at the file exe/full_physics/src/rtmod/ray_tau.F90
  // in tagged B2.06.02_plus_2.07_backport version to find the
  // original code and comment
  const double a0 = 6.02297e26;

  int nvar = pres->pressure_grid().value.number_variable();
  for(int i = 0; i < (int) alt.size(); ++i)
    nvar = std::max(nvar, alt[i]->gravity(pres->pressure_grid()(0)).value.number_variable());
  part_independent_wn.resize((int) alt.size(), pres->number_layer(), nvar);
  for(int i = 0; i < part_independent_wn.cols(); ++i)
    for(int j = 0; j < part_independent_wn.rows(); ++j) {
      AutoDerivativeWithUnit<double> deltap = pres->pressure_grid()(i + 1) - 
	pres->pressure_grid()(i);
      AutoDerivativeWithUnit<double> play = 
	pres->pressure_grid()(i) + deltap / 2;
      part_independent_wn(j, i) = a0 * deltap.convert(units::Pa).value /
	(molar_weight_dry_air * alt[j]->gravity(play).convert(Unit("m/s^2")).value);
    }
  cache_is_stale = false;
}

