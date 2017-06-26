// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "tccon_apriori.h"
%}
%base_import(generic_object)
%import "meteorology.i"
%import "level_1b.i"
%import "fp_time.i"
%import "double_with_unit.i"
%import "pressure.i"
%fp_shared_ptr(FullPhysics::TcconApriori);
namespace FullPhysics {
class TcconApriori : public GenericObject {
public:
  TcconApriori(const boost::shared_ptr<Meteorology>& Met_file,
	       const boost::shared_ptr<Level1b>& L1b_file);
  TcconApriori(const boost::shared_ptr<Meteorology>& Met_file,
	       const boost::shared_ptr<Level1b>& L1b_file,
	       double Co2_ref);
  TcconApriori(const boost::shared_ptr<Meteorology>& Met_file,
	       const boost::shared_ptr<Level1b>& L1b_file,
	       double Co2_ref,
	       const Time& Ref_time);
  TcconApriori(const boost::shared_ptr<Meteorology>& Met_file,
	       const boost::shared_ptr<Level1b>& L1b_file,
	       double Co2_ref,
	       const Time& Ref_time,
	       DoubleWithUnit Rate_increase);
  TcconApriori(const boost::shared_ptr<Meteorology>& Met_file,
	       const boost::shared_ptr<Level1b>& L1b_file,
	       double Co2_ref,
	       const Time& Ref_time,
	       DoubleWithUnit Rate_increase,
	       DoubleWithUnit Age_air_pbl);
  TcconApriori(const boost::shared_ptr<Meteorology>& Met_file,
	       const boost::shared_ptr<Level1b>& L1b_file,
	       double Co2_ref,
	       const Time& Ref_time,
	       DoubleWithUnit Rate_increase,
	       DoubleWithUnit Age_air_pbl,
	       DoubleWithUnit Age_air_tropopause);
  TcconApriori(const boost::shared_ptr<Meteorology>& Met_file,
	       const boost::shared_ptr<Level1b>& L1b_file,
	       double Co2_ref,
	       const Time& Ref_time,
	       DoubleWithUnit Rate_increase,
	       DoubleWithUnit Age_air_pbl,
	       DoubleWithUnit Age_air_tropopause,
	       DoubleWithUnit Age_air_upper_stratosphere);
  double co2_vmr(double P) const;
  blitz::Array<double, 1> co2_vmr_grid(const Pressure& P) const;
  double tropopause_pressure() const;
  double planetary_boundary_layer_pressure() const;
  double fractional_amplitude_seasonal_cycle() const;
  DoubleWithUnit age_air(double P) const;
  std::string print_to_string() const;
};
}
