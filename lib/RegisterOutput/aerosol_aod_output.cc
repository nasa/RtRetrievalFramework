#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include "aerosol_aod_output.h"

using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(AerosolAodOutput, RegisterOutputBase)
.def(luabind::constructor<const boost::shared_ptr<Aerosol>&>())
.def(luabind::constructor<const boost::shared_ptr<Aerosol>&, bool>())
REGISTER_LUA_END()
#endif

/// The boundary for calculating retrieved_aerosol_aod_by_type_low
const double AerosolAodOutput::low_boundary = 800e2;

/// The boundary for calculating retrieved_aerosol_aod_by_type_high
const double AerosolAodOutput::high_boundary = 500e2;

AerosolAodOutput::AerosolAodOutput(const boost::shared_ptr<Aerosol>& A, 
		 bool Number_instead_of_name)
: number_instead_of_name(Number_instead_of_name)
{
  // Right now we only work with AerosolOptical. Not sure if this is
  // a *real* requirement, or just that our Aerosol base class should
  // have the additional functions AerosolOptical supplies. But for
  // now just fail miserably if we are using a different kind of
  // Aerosol. 
  a = boost::dynamic_pointer_cast<AerosolOptical>(A);
  if(!a)
    throw Exception("Currently only support AerosolOptical");
}

// See base class for description

void AerosolAodOutput::register_output(const boost::shared_ptr<Output>& out) const
{
  typedef boost::function<blitz::Array<double, 1> ()> ftype;
  double minv = std::numeric_limits<double>::min();
  double maxv = std::numeric_limits<double>::max();

  out->register_data_source("/Metadata/AerosolTypes", &AerosolOptical::aerosol_name_arr, a);

  typedef boost::function<double ()> func_type;
  std::vector<std::string> aerosol_names = a->aerosol_name();
  for(int aer_idx = 0; aer_idx < a->number_particle(); aer_idx++) {
    std::string aer_name; 
    if(number_instead_of_name) {
      aer_name = boost::lexical_cast<std::string>(aer_idx + 1);
    } else {
      aer_name = aerosol_names[aer_idx];
      boost::algorithm::to_lower(aer_name);
    }

    func_type func_all = 
      boost::bind(&AerosolOptical::aerosol_optical_depth, a, aer_idx, minv, maxv);
    out->register_data_source("/RetrievalResults/aerosol_" + aer_name + "_aod", func_all);

    func_type func_low = 
      boost::bind(&AerosolOptical::aerosol_optical_depth, a, aer_idx, low_boundary, maxv);
    out->register_data_source
      ("/RetrievalResults/aerosol_" + aer_name + "_aod_low", func_low);

    func_type func_mid = 
      boost::bind(&AerosolOptical::aerosol_optical_depth, a, aer_idx, high_boundary, low_boundary);
    out->register_data_source
      ("/RetrievalResults/aerosol_" + aer_name + "_aod_mid", func_mid);

    func_type func_high = 
      boost::bind(&AerosolOptical::aerosol_optical_depth, a, aer_idx, minv, high_boundary);
    out->register_data_source
      ("/RetrievalResults/aerosol_" + aer_name + "_aod_high", func_high);
  }

  func_type func_tot_all = 
    boost::bind(&AerosolOptical::aerosol_optical_depth_total, a, minv, maxv);
  out->register_data_source("/RetrievalResults/aerosol_total_aod", func_tot_all);

  func_type func_tot_low = 
    boost::bind(&AerosolOptical::aerosol_optical_depth_total, a, low_boundary, maxv);
  out->register_data_source
    ("/RetrievalResults/aerosol_total_aod_low", func_tot_low);

  func_type func_tot_mid = 
    boost::bind(&AerosolOptical::aerosol_optical_depth_total, a, high_boundary, low_boundary);
  out->register_data_source
    ("/RetrievalResults/aerosol_total_aod_mid", func_tot_mid);

  func_type func_tot_high = 
    boost::bind(&AerosolOptical::aerosol_optical_depth_total, a, minv, high_boundary);
  out->register_data_source
    ("/RetrievalResults/aerosol_total_aod_high", func_tot_high);

}

