// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "absorber_absco.h"
#include "absco_hdf.h"
#include "sub_state_vector_array.h"
%}
%base_import(pressure)
%base_import(temperature)
%base_import(absorber_vmr)
%base_import(altitude)
%base_import(absorber)
%base_import(observer)
%import "gas_absorption.i"
%import "constant.i"
%fp_shared_ptr(FullPhysics::AbsorberAbsco);

%rename(gas_absorption) gas_absorption_ptr;

namespace FullPhysics {
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
	const std::vector<boost::shared_ptr<GasAbsorption> >& Gas_absorption,
		const boost::shared_ptr<Constant>& C,
		int Nsub = 10);
  AutoDerivativeWithUnit<double> integrand_independent_wn
  (int Spec_index, int Species_index, const DoubleWithUnit& P)
    const;
  double integrand
  (double wn, double p, int Spec_index, int Species_index) const;
  double optical_depth_each_layer_direct_integrate(double wn, int Spec_index,
						   int Species_index, int
						   Layer_index,
						   double eps_abs = 0,
						   double eps_rel = 1e-4) const;
  blitz::Array<double, 2> 
  optical_depth_each_layer_direct_integrate(double wn, int Spec_index,
					    double eps_abs = 0,
					    double eps_rel = 1e-4) const;
  virtual void notify_add(StateVector& Sv);
  virtual void notify_remove(StateVector& Sv);
  virtual void notify_update(const StateVector& Sv);
  %python_attribute_derived(number_species, int)
  %python_attribute_derived(number_spectrometer, int)
  %python_attribute_derived(number_layer, int)
  virtual std::string gas_name(int Species_index) const;
  virtual void notify_update(const Pressure& P);
  virtual void notify_update(const Temperature& T);
  virtual void notify_update(const Altitude& A);
  virtual void notify_update(const AbsorberVmr& A);
  virtual ArrayAd<double, 2> 
    optical_depth_each_layer(double wn, int spec_index) const;
  %python_attribute(specific_humidity_layer, ArrayAdWithUnit<double, 1>)
  %python_attribute(dry_air_molecular_density_layer, ArrayAdWithUnit<double, 1>)
  %python_attribute(dry_air_column_thickness_layer, ArrayAdWithUnit<double, 1>)
  %python_attribute(wet_air_column_thickness_layer, ArrayAdWithUnit<double, 1>)
  %python_attribute(pressure_weighting_function_layer, ArrayAd<double, 1>)
  %python_attribute(pressure_weighting_function_grid, ArrayAd<double, 1>)
  ArrayAdWithUnit<double, 1> 
  gas_column_thickness_layer(const std::string& Gas_name) const;
  AutoDerivativeWithUnit<double> 
  gas_total_column_thickness(const std::string& Gas_name) const;
  virtual AutoDerivative<double> xgas(const std::string& Gas_name) const;
  AutoDerivative<double> average_vmr(const std::string& Gas_name) const;
  virtual boost::shared_ptr<Absorber> clone() const;
  virtual boost::shared_ptr<Absorber> clone
   (const boost::shared_ptr<Pressure>& Press,
    const boost::shared_ptr<Temperature>& Temp,
    const std::vector<boost::shared_ptr<Altitude> >& Alt) const;
  virtual boost::shared_ptr<AbsorberVmr> absorber_vmr(const std::string& gas_name) const;
  virtual boost::shared_ptr<GasAbsorption> gas_absorption_ptr(const std::string& Gas_name) const;
  %python_attribute(pressure_sublayer, ArrayWithUnit<double, 1>)
  %python_attribute(temperature_sublayer, ArrayAdWithUnit<double, 1>)
  %python_attribute(h2o_vmr_sublayer, ArrayAdWithUnit<double, 1>)
  ArrayAdWithUnit<double, 1> vmr_sublayer(const std::string& Gas_name) const;
  ArrayAdWithUnit<double, 1> gravity_sublayer(int Spec_index) const;
};
}
