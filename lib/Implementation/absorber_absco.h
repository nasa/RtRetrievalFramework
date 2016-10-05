#ifndef ABSORBER_ABSCO_H
#define ABSORBER_ABSCO_H
#include "absorber.h"
#include "absorber_vmr.h"
#include "state_vector.h"
#include "pressure.h"
#include "pressure_level_input.h"
#include "absco.h"
#include "temperature.h"
#include "altitude.h"
#include "constant.h"
#include "fp_exception.h"
#include "fp_gsl_integrate.h"
#include <vector>

namespace FullPhysics {
/****************************************************************//**
  This class maintains the absorber portion of the state. This
  particular implementation uses the GasAbsorption classes for 
  calculating the gas absorption (e.g, the Absco tables).
*******************************************************************/

class AbsorberAbsco: public Absorber,
                     public Observer<AbsorberVmr>,
		     public Observer<Pressure>,
		     public Observer<Temperature>, 
		     public Observer<Altitude> {
public:
  AbsorberAbsco(const std::vector<boost::shared_ptr<AbsorberVmr> > Vmr,
		const boost::shared_ptr<Pressure>& Press,	   
		const boost::shared_ptr<Temperature>& Temp,
		const std::vector<boost::shared_ptr<Altitude> >& Alt,
		const std::vector<boost::shared_ptr<GasAbsorption> >& 
		Gas_absorption,
		const boost::shared_ptr<Constant>& C,
		int Nsub = 10);
  virtual ~AbsorberAbsco() {}
  virtual void notify_add(StateVector& Sv);
  virtual void notify_remove(StateVector& Sv);
  // We use attach_notify to directly attach the various object that
  // AbsorberAbsco contains. This means we don't need to do anything with
  // changes to the StateVector in this class, it is already handled
  // by the objects we contain.
  virtual void notify_update(const StateVector& Sv) 
  { notify_update_do(*this); }
  virtual int number_species() const {return (int) vmr.size(); }
  virtual int number_spectrometer() const {return (int) alt.size();}
  virtual int number_layer() const { return press->number_layer(); }
  virtual std::string gas_name(int Species_index) const 
  { range_check(Species_index, 0, number_species());
    return vmr[Species_index]->gas_name();
  }

//-----------------------------------------------------------------------
/// For performance, we cache some data as we calculate it. This
/// becomes stale when the pressure is changed, so we observe press
/// and mark the cache when it changes. 
//-----------------------------------------------------------------------

  virtual void notify_update(const Pressure& P) 
  { 
    cache_tau_gas_stale = true;
    notify_update_do(*this);
  }
  virtual void notify_update(const Temperature& T)
  {
    cache_tau_gas_stale = true;
    notify_update_do(*this);
  }
  virtual void notify_update(const Altitude& A)
  {
    cache_tau_gas_stale = true;
    notify_update_do(*this);
  }
  virtual void notify_update(const AbsorberVmr& A) 
  { 
    cache_tau_gas_stale = true;
    notify_update_do(*this);
  }
  virtual ArrayAd<double, 2> 
    optical_depth_each_layer(double wn, int spec_index) const;
  virtual blitz::Array<double, 2> 
    optical_depth_each_layer_nder(double wn, int spec_index) const;
  ArrayAdWithUnit<double, 1> specific_humidity_layer() const;
  ArrayAdWithUnit<double, 1> dry_air_molecular_density_layer() const;
  ArrayAdWithUnit<double, 1> dry_air_column_thickness_layer() const;
  ArrayAdWithUnit<double, 1> wet_air_column_thickness_layer() const;
  ArrayAd<double, 1> pressure_weighting_function_layer() const; 
  ArrayAd<double, 1> pressure_weighting_function_grid() const;
  ArrayAdWithUnit<double, 1> 
  gas_column_thickness_layer(const std::string& Gas_name) const;
  AutoDerivativeWithUnit<double> 
  gas_total_column_thickness(const std::string& Gas_name) const;

  virtual AutoDerivative<double> xgas(const std::string& Gas_name) const;

  AutoDerivative<double> total_number_density(const std::string& Gas_name) const;
  AutoDerivative<double> average_vmr(const std::string& Gas_name) const;

  virtual void print(std::ostream& Os) const;
  virtual boost::shared_ptr<Absorber> clone() const;
  virtual boost::shared_ptr<Absorber> clone
   (const boost::shared_ptr<Pressure>& Press,
    const boost::shared_ptr<Temperature>& Temp,
    const std::vector<boost::shared_ptr<Altitude> >& Alt) const;
  virtual boost::shared_ptr<AbsorberVmr> absorber_vmr(const std::string& Gas_name) const;
  const Pressure& pressure() const {return *press;}

  boost::shared_ptr<GasAbsorption> gas_absorption_ptr
  (const std::string& Gas_name) const;

//----------------------------------------------------------------
/// Return the pressure we use for each sublayer. This is meant for
/// diagnostic purposes.
//----------------------------------------------------------------

  ArrayWithUnit<double, 1> pressure_sublayer() const
  { fill_tau_gas_cache(); return psub;}

  ArrayAdWithUnit<double, 1> temperature_sublayer() const;
  ArrayAdWithUnit<double, 1> h2o_vmr_sublayer() const;
  ArrayAdWithUnit<double, 1> gravity_sublayer(int Spec_index) const;
  ArrayAdWithUnit<double, 1> vmr_sublayer(const std::string& Gas_name) const;
  AutoDerivativeWithUnit<double> integrand_independent_wn
  (int Spec_index, int Species_index, const DoubleWithUnit& P)
    const;
  double integrand
  (double wn, double p, int Spec_index, int Species_index) const;
  double optical_depth_each_layer_direct_integrate(double wn, int Spec_index,
						   int Species_index, int
						   Layer_index,
						   double eps_abs = 0,
						   double eps_rel = 1e-3) const;
  blitz::Array<double, 2> 
  optical_depth_each_layer_direct_integrate(double wn, int Spec_index,
					    double eps_abs = 0,
					    double eps_rel = 1e-3) const;
private:
  // Objects used to calculate the integrand
  boost::shared_ptr<Pressure> press;
  boost::shared_ptr<Temperature> temp;
  std::vector<boost::shared_ptr<Altitude> > alt;
  std::vector<boost::shared_ptr<AbsorberVmr> > vmr;
  std::vector<boost::shared_ptr<GasAbsorption> > gas_absorption;

  // Used for optical_depth_each_layer_direct_integrate. We don't use
  // this function when using AbsorberAbsco normally, but this is a
  // useful test function to see how well our approximate integral in
  // optical_depth_each_layer is working.
  GslIntegrate intg;

  // The species index that corresponds to water.
  int h2o_index;

  // Constants used to get things like avogadro_constant and
  // molar_weight_water. 
  boost::shared_ptr<Constant> c;

  // Number of sub layers to use in gas calculation.
  int nsub;

  // Cache a good portion of the tau_gas calculation, everything that
  // is not dependent on wn.

  mutable bool cache_tau_gas_stale; // Indicate if cache data is stale
  mutable std::vector<boost::shared_ptr<AbscoInterpolator> > absco_interp;
				// Used to quickly get absorption
				// coefficients for integral
  mutable ArrayAdWithUnit<double, 1> pgrid;
                                // Pressure grid on each level.
  mutable ArrayWithUnit<double, 1> psub;
				// Pressures that we calculate
				// integrand on
  mutable std::vector<blitz::Range> layer_range;
				// For each layer, give the Range in
				// psub that fits in that layer.
  mutable std::vector<ArrayWithUnit<double, 1> > weight;
                                // Weight for each layer.
  mutable ArrayAdWithUnit<double, 3> integrand_independent_wn_sub;
				// This is the integrand independent
				// of wavenumber wn at each psub.
				// This is indexed by Spec_index,
				// Species_index, sublayer index

  // The various jacobians we need to track tend to be fairly
  // sparse. To speed up calculation, we use a spare implementation of
  // these arrays.
  
  // For each layer, store the nonzero columns, and the gradient of p1
  // and p2 compressed to those nonzero columns.
  mutable std::vector<std::vector<int> > pressure_nonzero_column;
  mutable std::vector<blitz::Array<double, 1> > p1_grad;
  mutable std::vector<blitz::Array<double, 1> > p2_grad;

  // Same thing for integrand_independent_wn_sub
  mutable std::vector<std::vector<int> > integrand_nonzero_column;
  mutable std::vector<blitz::Array<double, 4> > integrand_jac;

  // And for Absco (only Absco is not by layer)
  mutable std::vector<int> absco_nonzero_column;


  // Scratch variable used to calculate taug. We keep this around so
  // we don't keep recreating the Array (so this is to improve performance)
  mutable ArrayAd<double, 2> taug;

  blitz::Array<double, 2> 
    tau_gas_nder(double wn, int spec_index) const;
  ArrayAd<double, 2> 
    tau_gas_der(double wn, int spec_index) const;
  void fill_tau_gas_cache() const;
  void create_sublayer() const;

  // Short functions to get the temperature, gravity, and mixing
  // levels at a particular pressure. Since we use these in a few
  // places, it is worth wrapping up into shorter functions rather
  // than repeating thing code in multiple places.
  AutoDerivativeWithUnit<double> 
  temp_func(const DoubleWithUnit& P) const
  { return 
      temp->temperature(AutoDerivativeWithUnit<double>(P.value, P.units)); }
  AutoDerivativeWithUnit<double>
  h2o_vmr_func(const DoubleWithUnit& P) const
  {     
    if(h2o_index < 0)
      return AutoDerivativeWithUnit<double>(0, units::dimensionless);
    else
      return AutoDerivativeWithUnit<double>(vmr[h2o_index]->volume_mixing_ratio(P.convert(Unit("Pa")).value), units::dimensionless);
  }
  AutoDerivativeWithUnit<double>
  vmr_func(int Species_index, const DoubleWithUnit& P) const
  { return AutoDerivativeWithUnit<double>(vmr[Species_index]->volume_mixing_ratio(P.convert(Unit("Pa")).value), units::dimensionless); }
  AutoDerivativeWithUnit<double>
  gravity_func(int Spec_index, const DoubleWithUnit& P) const
  { return alt[Spec_index]->gravity(AutoDerivativeWithUnit<double>(P.value, P.units)); }

};
}
#endif
