#include "absorber_absco_output.h"
#include "absorber_vmr_fixed_level.h"
#include "absorber_vmr_fixed_level_scaled.h"
#include "fill_value.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
// Lua doesn't know to cast a pointer type of base class to a derived class.
// Add a conversion routine.
boost::shared_ptr<RegisterOutputBase> abs_absco_create
(const boost::shared_ptr<Absorber>& A,
 const SpectralBound& Sb)
{
  return boost::shared_ptr<RegisterOutputBase>
    (new AbsorberAbscoOutput
     (boost::dynamic_pointer_cast<AbsorberAbsco>(A), Sb));
}
REGISTER_LUA_DERIVED_CLASS(AbsorberAbscoOutput, RegisterOutputBase)
.scope
[
 luabind::def("create", &abs_absco_create)
]
REGISTER_LUA_END()
#endif

// Helper class that gets Absco scaling information
class GasAbscoOutputHelper {
public:
  GasAbscoOutputHelper(const boost::shared_ptr<Absco>&
		       Ab,
		       const SpectralBound& Sb)
    : absco(Ab),
      sb(Sb)
  {}
  blitz::Array<double, 1> absco_scale() const
  { blitz::Array<double, 1> res(sb.number_spectrometer());
    for(int i = 0; i < res.rows(); ++i) {
      // Only convert non zero values, other wise we will get floating point errors
      // in single band mode
      if(abs(sb.upper_bound(i).value - sb.lower_bound(i).value) > 0) {
          DoubleWithUnit wn = sb.center(i, Unit("cm^-1"));
          res(i) = absco->table_scale(wn.value);
      } else {
          res(i) = 0.0;
      }
    }
    return res;
  }
private:
  boost::shared_ptr<Absco> absco;
  SpectralBound sb;
};

// Helper class that gets XCO2 value
class AbsorberAbscoOutputHelper {
public:
  AbsorberAbscoOutputHelper(const boost::shared_ptr<AbsorberAbsco>& Abs)
    : abs(Abs) {}
  double xco2_value() const { return abs->xgas("CO2").value(); }
  double total_column_thickness_o2_value() const { return abs->gas_total_column_thickness("O2").convert(Unit("m^-2")).value.value(); }
  double total_column_thickness_co2_value() const { return abs->gas_total_column_thickness("CO2").convert(Unit("m^-2")).value.value(); }
  double total_column_thickness_h2o_value() const { return abs->gas_total_column_thickness("H2O").convert(Unit("m^-2")).value.value(); }
  blitz::Array<double, 1>layer_column_thickness_h2o_value() const { return abs->gas_column_thickness_layer("H2O").convert(Unit("m^-2")).value.value(); }
  double average_vmr_o2_value() const { return abs->average_vmr("O2").value(); }
  blitz::Array<double, 1> pressure_weighting_function_grid_value() const
  { return abs->pressure_weighting_function_grid().value(); }
  blitz::Array<double, 1> dry_air_column_thickness_value() const
  { return abs->dry_air_column_thickness_layer().convert(Unit("m^-2")).value.value(); }
  blitz::Array<double, 1> wet_air_column_thickness_value() const
  { return abs->wet_air_column_thickness_layer().convert(Unit("m^-2")).value.value(); }
private:
  boost::shared_ptr<AbsorberAbsco> abs;
};


// See base class for description

void AbsorberAbscoOutput::register_output_apriori(const boost::shared_ptr<Output>& out) const
{
  // Freeze the pressure state
  boost::shared_ptr<AbsorberAbsco> afreeze = 
    boost::dynamic_pointer_cast<AbsorberAbsco>(a->clone());
  boost::shared_ptr<AbsorberAbscoOutputHelper> h(new AbsorberAbscoOutputHelper(afreeze));
  if(afreeze->gas_index("CO2") != -1)
    out->register_data_source("/RetrievalResults/xco2_apriori",
			     &AbsorberAbscoOutputHelper::xco2_value, h);

  if(afreeze->gas_index("O2") != -1)
    out->register_data_source("/RetrievalResults/apriori_o2_column",
			     &AbsorberAbscoOutputHelper::total_column_thickness_o2_value, h);

}

void AbsorberAbscoOutput::register_output(const boost::shared_ptr<Output>& out) const
{
  boost::shared_ptr<AbsorberAbscoOutputHelper> h(new AbsorberAbscoOutputHelper(a));
  if(a->gas_index("CO2") != -1) {
    out->register_data_source("/RetrievalResults/xco2",
			     &AbsorberAbscoOutputHelper::xco2_value, h);
    out->register_data_source_pad
      ("/RetrievalResults/xco2_pressure_weighting_function",
       &AbsorberAbscoOutputHelper::pressure_weighting_function_grid_value, h, 
       num_level, fill_value<double>());

    out->register_data_source("/RetrievalResults/retrieved_co2_column",
			     &AbsorberAbscoOutputHelper::total_column_thickness_co2_value, h);
  }

  if(a->gas_index("O2") != -1) {
    out->register_data_source("/RetrievalResults/retrieved_o2_column",
			     &AbsorberAbscoOutputHelper::total_column_thickness_o2_value, h);
 
    out->register_data_source("/Metadata/VMRO2", 
			     &AbsorberAbscoOutputHelper::average_vmr_o2_value, h);
  }

  if(a->gas_index("H2O") != -1) {
    out->register_data_source("/RetrievalResults/retrieved_h2o_column",
			     &AbsorberAbscoOutputHelper::total_column_thickness_h2o_value, h);

    out->register_data_source_pad("/RetrievalResults/retrieved_h2o_column_layer_thickness",
				 &AbsorberAbscoOutputHelper::layer_column_thickness_h2o_value, h, num_level - 1 , fill_value<double>());
  }

  // Wet and dry air mass, H2O used in calculation of results
  out->register_data_source_pad("/RetrievalResults/retrieved_dry_air_column_layer_thickness",
			       &AbsorberAbscoOutputHelper::dry_air_column_thickness_value, h, num_level - 1 , fill_value<double>());
    
  out->register_data_source_pad("/RetrievalResults/retrieved_wet_air_column_layer_thickness",
			       &AbsorberAbscoOutputHelper::wet_air_column_thickness_value, h, num_level - 1, fill_value<double>());
  for(int i = 0; i < a->number_species(); ++i) {
    std::string gas_name = a->gas_name(i);
    boost::shared_ptr<Absco> ab =
      boost::dynamic_pointer_cast<Absco>(a->gas_absorption_ptr(gas_name));
    if(ab) {
      boost::shared_ptr<GasAbscoOutputHelper> 
	h2(new GasAbscoOutputHelper(ab, sb));
      out->register_data_source("/Metadata/Absco" + gas_name + "Scale", 
				&GasAbscoOutputHelper::absco_scale, h2);
    }
  }
}
