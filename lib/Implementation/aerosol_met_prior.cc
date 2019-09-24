#include "aerosol_met_prior.h"
#include "aerosol_optical.h"
#include "aerosol_property_hdf.h"
#include "aerosol_property_rh_hdf.h"
#include "aerosol_shape_gaussian.h"
#include "initial_guess_value.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
// Luabind is restricted to 10 arguments in a function. The
// constructor for AerosolMetPrior takes too many. As an easy
// workaround, just stuff a number of the values into an array.
boost::shared_ptr<AerosolMetPrior> aerosol_met_prior_create(
const boost::shared_ptr<Meteorology>& Met_file,
const HdfFile& Aerosol_property,
const boost::shared_ptr<Pressure> &Press,
const boost::shared_ptr<RelativeHumidity> &Rh,
const blitz::Array<double, 2>& Aerosol_cov,
const blitz::Array<double, 1>& val)
{
  if(boost::dynamic_pointer_cast<OcoMetFile>(Met_file) == 0)
    throw Exception("Need a OcoMetFile");
  return boost::shared_ptr<AerosolMetPrior>
    (new AerosolMetPrior(*boost::dynamic_pointer_cast<OcoMetFile>(Met_file),
			 Aerosol_property, Press, Rh, Aerosol_cov,
			 val(0), static_cast<int>(val(1)),
			 static_cast<int>(val(2)), val(3),
			 static_cast<int>(val(4)) == 1,
			 static_cast<int>(val(5)) == 1));
}
REGISTER_LUA_CLASS(AerosolMetPrior)
.def("aerosol", &AerosolMetPrior::aerosol)
.def("initial_guess", &AerosolMetPrior::initial_guess)
.def("add_aerosol", &AerosolMetPrior::add_aerosol)
.def("number_merra_particle", &AerosolMetPrior::number_merra_particle)
.scope
[
 luabind::def("create", &aerosol_met_prior_create)
]
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
/// \param Met_file The Meteorological file. 
/// \param Aerosol_property The Aerosol property file
/// \param Press The Pressure object that gives the pressure grid.
/// \param Rh The RelativeHumidity object that gives the relative humidity.
/// \param Aerosol_cov The covariance matrix to use for each Aerosol.
/// \param Exp_aod Threshold for explained fraction of AOD
/// \param Min_types Minimum number of types to be selected
/// \param Max_types Maximum number of types to be selected
/// \param Linear_aod If true, then use linear aod rather than
///    logarithmic aod.
/// \param Relative_humidity_aerosol If true, then use AerosolPropertyRhHdf
///   instead of AerosolPropertyHdf which includes interpolation
///    by relative humidity.
/// \param Max_residual Maximum resdidual in total AOD (either >
///    threshold of fraction or < max residual will suffice 
/// \param Reference_wn The wavenumber that Aext is given for. This
///    is optional, the default value matches the reference band given
///    in the ATB.
//-----------------------------------------------------------------------

AerosolMetPrior::AerosolMetPrior
(const OcoMetFile& Met_file,
 const HdfFile& Aerosol_property,
 const boost::shared_ptr< Pressure > &Press, 
 const boost::shared_ptr<RelativeHumidity> &Rh,
 const blitz::Array<double, 2>& Aerosol_cov,
 double Exp_aod,
 int Min_types,
 int Max_types,
 bool Linear_aod,
 bool Relative_humidity_aerosol,
 double Max_residual,
 double Reference_wn)
: linear_aod(Linear_aod),
  rh_aerosol(Relative_humidity_aerosol),
  press(Press),
  rh(Rh),
  ref_wn(Reference_wn),
  met_fname(Met_file.file_name()),
  prop_fname(Aerosol_property.file_name())
{
  Array<double, 1> aod_frac = Met_file.read_array("/Aerosol/composite_aod_fraction_met");
  Array<double, 2> aod_gauss =  Met_file.read_array_2d("/Aerosol/composite_aod_gaussian_met");
  Array<int, 1> aod_sort_index = Met_file.read_array_int("/Aerosol/composite_aod_sort_index_met");
  Array<std::string, 1> comp_name = 
    Met_file.hdf_file().read_field<std::string, 1>("/Metadata/CompositeAerosolTypes");
  double total_aod = sum(aod_gauss(Range::all(), 3));
  ig.reset(new CompositeInitialGuess);
  // Loop over composite types until thresholds are reached:
  number_merra_particle_ = 0;
  for(int counter = 0; counter < aod_sort_index.shape()[0]; ++counter) {
    int i = aod_sort_index(counter);
    bool select;
    if(counter < Min_types)
      select = true;
    else if(aod_frac(counter) > Exp_aod ||
	    (1 - aod_frac(counter)) * total_aod < Max_residual ||
	    counter >= Max_types)
      select = false;
    else
      select = true;
    if(select) {
      ++number_merra_particle_;
      double aod = aod_gauss(i, 3);
      double peak_height = aod_gauss(i, 1);
      double peak_width = aod_gauss(i, 2);
      std::string aname = comp_name(i);
      Array<bool, 1> flag(3);
      Array<double, 1> coeff(3);
      flag = true, true, true;
      if(linear_aod)
	coeff = aod, peak_height, peak_width;
      else
	coeff = log(aod), peak_height, peak_width;
      aext.push_back
	(boost::shared_ptr<AerosolExtinction>
	 (new AerosolShapeGaussian(press, flag, coeff, aname, linear_aod)));
      if(rh_aerosol) 
	aprop.push_back
	  (boost::shared_ptr<AerosolProperty>
	   (new AerosolPropertyRhHdf(Aerosol_property, aname + "/Properties", 
				     press, rh)));
      else
	aprop.push_back
	  (boost::shared_ptr<AerosolProperty>
	   (new AerosolPropertyHdf(Aerosol_property, aname + "/Properties", 
				   press)));
      boost::shared_ptr<InitialGuessValue> nig(new InitialGuessValue);
      nig->apriori(coeff);
      nig->apriori_covariance(Aerosol_cov);
      ig->add_builder(nig);
    }
  }
}


//-----------------------------------------------------------------------
/// Return the aerosol setup generated from this class.
//-----------------------------------------------------------------------

boost::shared_ptr<Aerosol> AerosolMetPrior::aerosol() const 
{ 
  if(!aerosol_)
    aerosol_.reset(new AerosolOptical(aext, aprop, press, rh, ref_wn));
  return aerosol_;
}

//-----------------------------------------------------------------------
/// Add in an aerosol that comes from some source other than Dijian's
/// files (e.g., hardcoded aerosols to include at all times)
//-----------------------------------------------------------------------

void AerosolMetPrior::add_aerosol
(const boost::shared_ptr<AerosolExtinction>& Aext,
 const boost::shared_ptr<AerosolProperty>& Aprop,
 const boost::shared_ptr<InitialGuessBuilder>& Ig)
{
  aext.push_back(Aext);
  aprop.push_back(Aprop);
  ig->add_builder(Ig);
}

//-----------------------------------------------------------------------
/// Print to stream.
//-----------------------------------------------------------------------

void AerosolMetPrior::print(std::ostream& Os) const
{
  Os << "AerosolMetPrior: ("
     << (linear_aod ? "Linear, " : "Logarithmic, ")
     << (rh_aerosol ? "RH aerosol" : "Not RH aerosol")
     << ")\n"
     << "  Coefficient:\n"
     << "  Met file name:             " << met_fname << "\n"
     << "  Aerosol property file name: " << prop_fname << "\n";
}
