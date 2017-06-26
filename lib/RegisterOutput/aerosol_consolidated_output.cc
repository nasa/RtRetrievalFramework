#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include "aerosol_consolidated_output.h"
#include "aerosol_extinction_imp_base.h"
#include "fill_value.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(AerosolConsolidatedOutput, RegisterOutputBase)
.def(luabind::constructor<const boost::shared_ptr<Aerosol>&, const std::vector<std::string>&>())
REGISTER_LUA_END()
#endif

class AerosolOutputHelper {
public:
    AerosolOutputHelper(const boost::shared_ptr<AerosolOptical>& Aerosol, const std::vector<std::string>& All_aer_names)
        : aerosol(Aerosol), all_aer_names(All_aer_names) { initialize(); }

    Array<double, 2> aerosol_aod_matrix() const {
        // Limits for low/high aerosol types
        double minv = std::numeric_limits<double>::min();
        double maxv = std::numeric_limits<double>::max();

        Array<double, 2> aer_aod_matrix_(all_aer_names.size(), 4);
        aer_aod_matrix_ = fill_value<double>(); 

        for(int out_aer_idx = 0; out_aer_idx < all_aer_names.size(); out_aer_idx++) {
            auto match_iter = std::find(ret_aer_names.begin(), ret_aer_names.end(), all_aer_names[out_aer_idx]);

            if(match_iter != ret_aer_names.end()) {
                int in_aer_idx = std::distance(ret_aer_names.begin(), match_iter);

                // total
                aer_aod_matrix_(out_aer_idx, 0) = aerosol->aerosol_optical_depth(in_aer_idx, minv, maxv); 

                // low
                aer_aod_matrix_(out_aer_idx, 1) = aerosol->aerosol_optical_depth(in_aer_idx, low_boundary, maxv);

                // mid
                aer_aod_matrix_(out_aer_idx, 2) = aerosol->aerosol_optical_depth(in_aer_idx, high_boundary, low_boundary);

                // low
                aer_aod_matrix_(out_aer_idx, 3) = aerosol->aerosol_optical_depth(in_aer_idx, minv, high_boundary);
            }
        }
        return aer_aod_matrix_; 
    }
    
    Array<double, 2> aerosol_param_matrix() const { 
        Array<double, 2> aer_param_matrix_(all_aer_names.size(), max_num_param);
        aer_param_matrix_ = fill_value<double>();

        for(int out_aer_idx = 0; out_aer_idx < all_aer_names.size(); out_aer_idx++) {
            auto match_iter = std::find(ret_aer_names.begin(), ret_aer_names.end(), all_aer_names[out_aer_idx]);

            if(match_iter != ret_aer_names.end()) {
                int in_aer_idx = std::distance(ret_aer_names.begin(), match_iter);

                auto aer_ext = boost::dynamic_pointer_cast<AerosolExtinctionImpBase>(aerosol->aerosol_extinction(in_aer_idx));
                if (!aer_ext)
                    throw Exception("Could not ast aerosol extinction into AerosolExctinctionImpBase");

                Range param_r = Range(0, aer_ext->aerosol_parameter().rows()-1);
                aer_param_matrix_(out_aer_idx, param_r) = aer_ext->aerosol_parameter();
            }
        }
        return aer_param_matrix_; 
    }

    Array<double, 2> aerosol_param_uncert_matrix() const { 
        Array<double, 2> aer_param_uncert_matrix_(all_aer_names.size(), max_num_param);
        aer_param_uncert_matrix_ = fill_value<double>();

        for(int out_aer_idx = 0; out_aer_idx < all_aer_names.size(); out_aer_idx++) {
            auto match_iter = std::find(ret_aer_names.begin(), ret_aer_names.end(), all_aer_names[out_aer_idx]);

            if(match_iter != ret_aer_names.end()) {
                int in_aer_idx = std::distance(ret_aer_names.begin(), match_iter);

                auto aer_ext = boost::dynamic_pointer_cast<AerosolExtinctionImpBase>(aerosol->aerosol_extinction(in_aer_idx));
                if (!aer_ext)
                    throw Exception("Could not ast aerosol extinction into AerosolExctinctionImpBase");

                Range param_r = Range(0, aer_ext->aerosol_parameter().rows()-1);
                aer_param_uncert_matrix_(out_aer_idx, param_r) = aer_ext->aerosol_parameter_uncertainty();
            }
        }
        return aer_param_uncert_matrix_; 
    }

    blitz::Array<std::string, 1> aerosol_model_types() const { return aer_model_types_; }
    blitz::Array<int, 1> aerosol_type_retrieved() const { return type_retrieved_; }

    blitz::Array<std::string, 1> all_aerosol_names() const { 
        // Convert to blitz array for output
        blitz::Array<std::string, 1> out_aer_names(all_aer_names.size());
        for(int i; i < all_aer_names.size(); i++)
            out_aer_names(i) = all_aer_names[i];
        return out_aer_names;
    }

private:

    void initialize() 
    {
        ret_aer_names = aerosol->aerosol_name();

        aer_model_types_.resize(all_aer_names.size());
        aer_model_types_ = "";

        type_retrieved_.resize(all_aer_names.size());
        type_retrieved_ = 0;

        // Find the maximum number of parameters and also fill in the model types array
        max_num_param = 0;
        for(int out_aer_idx = 0; out_aer_idx < all_aer_names.size(); out_aer_idx++) {
            auto match_iter = std::find(ret_aer_names.begin(), ret_aer_names.end(), all_aer_names[out_aer_idx]);

            if(match_iter != ret_aer_names.end()) {
                int in_aer_idx = std::distance(ret_aer_names.begin(), match_iter);

                auto aer_ext = boost::dynamic_pointer_cast<AerosolExtinctionImpBase>(aerosol->aerosol_extinction(in_aer_idx));
                if (!aer_ext)
                    throw Exception("Could not ast aerosol extinction into AerosolExctinctionImpBase");
                max_num_param = max(max_num_param, aer_ext->aerosol_parameter().rows());

                aer_model_types_(out_aer_idx) = aer_ext->model_short_name();
                type_retrieved_(out_aer_idx) = 1;
            }
        }
    }

private:
    /// The boundary for calculating retrieved_aerosol_aod_by_type_low
    const double low_boundary = 800e2;

    /// The boundary for calculating retrieved_aerosol_aod_by_type_high
    const double high_boundary = 500e2;

    boost::shared_ptr<AerosolOptical> aerosol;
    std::vector<std::string> all_aer_names;
    std::vector<std::string> ret_aer_names;
    int max_num_param;

    Array<std::string, 1> aer_model_types_;
    Array<int, 1> type_retrieved_;
  
};

AerosolConsolidatedOutput::AerosolConsolidatedOutput(const boost::shared_ptr<Aerosol>& Aerosol, const std::vector<std::string>& All_aer_names)
: all_aer_names(All_aer_names)
{
    // Right now we only work with AerosolOptical. Not sure if this is
    // a *real* requirement, or just that our Aerosol base class should
    // have the additional functions AerosolOptical supplies. But for
    // now just fail miserably if we are using a different kind of
    // Aerosol. 
    aerosol = boost::dynamic_pointer_cast<AerosolOptical>(Aerosol);
    if(!aerosol)
        throw Exception("Currently only support AerosolOptical");
}



// See base class for description
void AerosolConsolidatedOutput::register_output_apriori(const boost::shared_ptr<Output>& out) const
{
    // Freeze the state
    boost::shared_ptr<AerosolOptical> afreeze = 
        boost::dynamic_pointer_cast<AerosolOptical>(aerosol->clone());
    boost::shared_ptr<AerosolOutputHelper> h(new AerosolOutputHelper(afreeze, all_aer_names));
    
    out->register_data_source("/RetrievalResults/aerosol_param_apriori", &AerosolOutputHelper::aerosol_param_matrix, h);
}

void AerosolConsolidatedOutput::register_output(const boost::shared_ptr<Output>& out) const
{
    boost::shared_ptr<AerosolOutputHelper> h(new AerosolOutputHelper(aerosol, all_aer_names));

    out->register_data_source("/Metadata/AllAerosolTypes", &AerosolOutputHelper::all_aerosol_names, h);

    out->register_data_source("/RetrievalResults/aerosol_model", &AerosolOutputHelper::aerosol_model_types, h);
    out->register_data_source("/RetrievalResults/aerosol_type_retrieved", &AerosolOutputHelper::aerosol_type_retrieved, h);

    out->register_data_source("/RetrievalResults/aerosol_aod", &AerosolOutputHelper::aerosol_aod_matrix, h);
    out->register_data_source("/RetrievalResults/aerosol_param", &AerosolOutputHelper::aerosol_param_matrix, h);
    out->register_data_source("/RetrievalResults/aerosol_param_uncert", &AerosolOutputHelper::aerosol_param_uncert_matrix, h);

    // Also include total optical depth. This is a duplicate of what
    // we had in AerosolAodOutput. Hopefully this doesn't cause any
    // problems, I don't believe we plan on having this both on at the
    // same time. But if we do, we might need to add some simple logic
    // here (e.g., add a flag to AerosolConsolidatedOutput to turn
    // this on or off here).
    double minv = std::numeric_limits<double>::min();
    double maxv = std::numeric_limits<double>::max();
    boost::function<double ()> func_tot_all = 
      boost::bind(&AerosolOptical::aerosol_optical_depth_total, aerosol,
		  minv, maxv);
    out->register_data_source("/RetrievalResults/aerosol_total_aod", func_tot_all);
}

