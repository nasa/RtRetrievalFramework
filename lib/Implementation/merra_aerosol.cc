#include "merra_aerosol.h"
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
// constructor for MerraAerosol takes too many. As an easy
// workaround, just stuff a number of the values into an array.
boost::shared_ptr<MerraAerosol> merra_aerosol_create(
const HdfFile& Merra_climatology,
const HdfFile& Aerosol_property,
DoubleWithUnit Latitude, DoubleWithUnit Longitude,
const boost::shared_ptr<Pressure> &Press,
const boost::shared_ptr<RelativeHumidity> &Rh,
const blitz::Array<double, 2>& Aerosol_cov,
const blitz::Array<double, 1>& val)
{
  return boost::shared_ptr<MerraAerosol>
    (new MerraAerosol(Merra_climatology, Aerosol_property,
		      Latitude, Longitude, Press, Rh, Aerosol_cov,
		      val(0), val(1), static_cast<int>(val(2)),
		      static_cast<int>(val(3)), val(4),
		      static_cast<int>(val(5)) == 1,
		      static_cast<int>(val(6)) == 1));
}
REGISTER_LUA_CLASS(MerraAerosol)
.def("aerosol", &MerraAerosol::aerosol)
.def("initial_guess", &MerraAerosol::initial_guess)
.def("add_aerosol", &MerraAerosol::add_aerosol)
.def("number_merra_particle", &MerraAerosol::number_merra_particle)
.scope
[
 luabind::def("create", &merra_aerosol_create)
]
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
/// \param Merra_climatology The Merra climatology file
/// \param Aerosol_property The Aerosol property file
/// \param Latitude The latitude of the ground point
/// \param Longitude The longitude of the ground point
/// \param Press The Pressure object that gives the pressure grid.
/// \param Rh The RelativeHumidity object that gives the relative humidity.
/// \param Aerosol_cov The covariance matrix to use for each Aerosol.
/// \param Max_aod Maximum AOD cap for each composite type
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

MerraAerosol::MerraAerosol
(const HdfFile& Merra_climatology,
 const HdfFile& Aerosol_property,
 DoubleWithUnit Latitude, DoubleWithUnit Longitude,
 const boost::shared_ptr< Pressure > &Press, 
 const boost::shared_ptr<RelativeHumidity> &Rh,
 const blitz::Array<double, 2>& Aerosol_cov,
 double Max_aod,
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
  merra_fname(Merra_climatology.file_name()),
  prop_fname(Aerosol_property.file_name()),
  lat(Latitude),
  lon(Longitude)
{
  int lat_ind = file_index(Merra_climatology, "Latitude", Latitude);
  int lon_ind = file_index(Merra_climatology, "Longitude", Longitude, true);
  TinyVector<int,3> shp = Merra_climatology.read_shape<3>("/AOD_Fraction");
  Array<double, 3> aod_frac = 
    Merra_climatology.read_field<double, 3>("/AOD_Fraction", 
					    shape(0,lat_ind,lon_ind),
					    shape(shp[0],1,1));
  TinyVector<int,4> shp2 = 
    Merra_climatology.read_shape<4>("/COMPOSITE_AOD_GAUSSIAN");
  Array<double, 4> aod_gauss = 
    Merra_climatology.read_field<double, 4>("/COMPOSITE_AOD_GAUSSIAN", 
					    shape(0,0,lat_ind,lon_ind),
					    shape(shp2[0],shp2[1], 1,1));
  shp = Merra_climatology.read_shape<3>("/COMPOSITE_AOD_SORT_Index");
  Array<int, 3> aod_sort_index =
    Merra_climatology.read_field<int, 3>("/COMPOSITE_AOD_SORT_Index", 
					shape(0,lat_ind,lon_ind),
					shape(shp[0],1,1));
  Array<std::string, 1> comp_name = 
    Merra_climatology.read_field<std::string, 1>("/COMPOSITE_NAME");
  double total_aod = sum(aod_gauss(3, Range::all(), 0, 0));
  ig.reset(new CompositeInitialGuess);
  // Loop over composite types until thresholds are reached:
  number_merra_particle_ = 0;
  for(int counter = 0; counter < aod_sort_index.shape()[0]; ++counter) {
    int i = aod_sort_index(counter);
    bool select;
    if(counter < Min_types)
      select = true;
    else if(aod_frac(counter, 0, 0) > Exp_aod ||
	    (1 - aod_frac(counter, 0, 0)) * total_aod < Max_residual ||
	    counter >= Max_types)
      select = false;
    else
      select = true;
    if(select) {
      ++number_merra_particle_;
      double aod = aod_gauss(3, i);
      double peak_height = aod_gauss(1, i);
      double peak_width = aod_gauss(2, i);
      if(aod > Max_aod)
	aod = Max_aod;
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

boost::shared_ptr<Aerosol> MerraAerosol::aerosol() const 
{ 
  if(!aerosol_)
    aerosol_.reset(new AerosolOptical(aext, aprop, press, rh, ref_wn));
  return aerosol_;
}

//-----------------------------------------------------------------------
/// Add in an aerosol that comes from some source other than Dijian's
/// files (e.g., hardcoded aerosols to include at all times)
//-----------------------------------------------------------------------

void MerraAerosol::add_aerosol
(const boost::shared_ptr<AerosolExtinction>& Aext,
 const boost::shared_ptr<AerosolProperty>& Aprop,
 const boost::shared_ptr<InitialGuessBuilder>& Ig)
{
  aext.push_back(Aext);
  aprop.push_back(Aprop);
  ig->add_builder(Ig);
}

//-----------------------------------------------------------------------
/// Determine latitude index to use.
//-----------------------------------------------------------------------

int MerraAerosol::file_index(const HdfFile& Merra_climatology, 
			     const std::string& Field_name,
			     DoubleWithUnit& V, bool Wrap_around)
{
  Array<double, 1> value_data = 
    Merra_climatology.read_field<double, 1>("/" + Field_name);
  double value = V.convert(Unit("deg")).value;
  // Special handling at the edges for longitude
  if(Wrap_around &&
     value < value_data(0) && 
     value >= (value_data(value_data.rows() - 1) - 360)) {
    if(fabs(value - value_data(0)) <
       fabs(value - (value_data(value_data.rows() - 1) - 360)))
      return 0;
    else
      return value_data.rows() - 1;
  }
  if(Wrap_around &&
     value >= value_data(value_data.rows() - 1) &&
     value < (value_data(0) + 360)) {
    if(fabs(value - (value_data(0) + 360)) <
       fabs(value - value_data(value_data.rows() - 1)))
      return 0;
    else
      return value_data.rows() - 1;
  }
  if(value < value_data(0) || value >= value_data(value_data.rows() - 1)) {
    Exception e;
    e << Field_name << " " << V << " is out or the range covered by the\n"
      << "Merra_climatology file " << Merra_climatology.file_name();
    throw e;
  }
  if(value == value_data(0))
    return 0;
  double* f = std::lower_bound(value_data.data(), 
			       value_data.data() + value_data.rows(), 
			       value);
  int i = f - value_data.data();
  if(((value - value_data(i - 1)) / (value_data(i) - value_data(i - 1))) < 0.5)
    return i - 1;
  else
    return i;
}

//-----------------------------------------------------------------------
/// Print to stream.
//-----------------------------------------------------------------------

void MerraAerosol::print(std::ostream& Os) const
{
  Os << "MerraAerosol: ("
     << (linear_aod ? "Linear, " : "Logarithmic, ")
     << (rh_aerosol ? "RH aerosol" : "Not RH aerosol")
     << ")\n"
     << "  Coefficient:\n"
     << "  Merra file name:            " << merra_fname << "\n"
     << "  Aerosol property file name: " << prop_fname << "\n"
     << "  Latitude:                   " << lat << "\n"
     << "  Longitude:                  " << lon << "\n";
}
