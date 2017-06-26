#include "level_1b_fts.h"
#include "fp_exception.h"
// This won't be needed once version 3 becomes the default for boost filesystem
#define BOOST_FILESYSTEM_VERSION 3
#include <boost/filesystem.hpp>

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"

REGISTER_LUA_DERIVED_CLASS(Level1bFts, Level1b)
.def(luabind::constructor<const std::string&,
     const std::vector<std::string>&,
     const SpectralBound&>())
.def(luabind::constructor<const HdfFile&,
                          const std::vector<std::string>&,
                          const std::string& >())
.def(luabind::constructor<const HdfFile&,
                          const std::vector<std::string>& >())
.def("run_log", (const FtsRunLogRecord& (Level1bFts::*)(int) const) &Level1bFts::run_log)
.def("run_log", (const std::vector<FtsRunLogRecord>& (Level1bFts::*)() const) &Level1bFts::run_log)
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------

extern "C" {
  void jetspe(const char* specname_c,const int *speclen, const double *opd,
              const double *graw,
              const int *ifirst, const int *ilast, const int *possp,
              const int *bytepw,
              const double *nus, const double *nue,
              const int *apo_m, const int *interp, const float *foff,
              const float *res, float *slit, const int *mii,
              const int *nscycle, float *yobs, const int *mmp, int *nmp,
              double *nustrt, double *delwav, int *status);
  double riair(const double *w, const double *t,
               const double *p, const double *h);
}

//-----------------------------------------------------------------------
/// Constructor. This takes the run log file and the set of spectra to
/// use. We also give the spectral range. Spectral_range should be a
/// spectra.size() x 2 array, with the first column giving the start
/// of the frequency range and the second the end of it.
//-----------------------------------------------------------------------

Level1bFts::Level1bFts(const std::string& Runlog, 
                       const std::vector<std::string>& Spectra,
                       const ArrayWithUnit<double, 2>& Spectral_range)
{
  initialize(Runlog, Spectra, Spectral_range);
}

//-----------------------------------------------------------------------
/// Creates a L1B FTS object using spectral window file instead of 
/// an array of window ranges
//-----------------------------------------------------------------------

Level1bFts::Level1bFts(const std::string& Runlog,
                       const std::vector<std::string>& Spectra,
                       const SpectralBound& Spec_bound)
{
  ArrayWithUnit<double, 2> spec_range;
  spec_range.value.resize(Spec_bound.number_spectrometer(),2);
  spec_range.units = Spec_bound.lower_bound(0).units;
  for(int spec_idx = 0; spec_idx < Spec_bound.number_spectrometer(); 
      spec_idx++) {
    spec_range.value(spec_idx, 0) = 
      Spec_bound.lower_bound(spec_idx, spec_range.units).value;
    spec_range.value(spec_idx, 1) = 
      Spec_bound.upper_bound(spec_idx, spec_range.units).value;
  }
  initialize(Runlog, Spectra, spec_range);
}

//-----------------------------------------------------------------------
/// Initializes the object for the constructor
//-----------------------------------------------------------------------

void Level1bFts::initialize(const std::string& Runlogfile, 
                            const std::vector<std::string>& spectra,
                            const ArrayWithUnit<double, 2>& Spectral_range)
{

  if((int) spectra.size() != Spectral_range.value.rows())
    throw Exception("Spectra and Spectral_range need to be the same size");
  if(Spectral_range.value.cols() != 2)
    throw Exception("Spectral_range needs to have two columns");

  FtsRunLog ftr(Runlogfile);

  int num_spec = (int) Spectral_range.value.rows();
  freq_start.resize(num_spec);
  freq_end.resize(num_spec);
  freq_spacing.resize(num_spec);
  
  freq_spacing = 0;
  freq_start   = 0;
  freq_end     = 0;

  for(int specnum = 0; specnum < num_spec; ++specnum) {
    bool non_empty_window = Spectral_range.value(specnum, 1) - Spectral_range.value(specnum, 0) >  0.0;

    Array<double, 1> radiance_arr;
    Array<double, 1> uncert_arr;
    FtsRunLogRecord frlr;

    // Only process spectrometers that have non-empty ranges
    if(non_empty_window) {
      boost::filesystem::path path(spectra[specnum]);
      std::string ext = boost::filesystem::extension(spectra[specnum]);
      string basespectrum = boost::filesystem::basename (path);
      frlr = ftr.read(basespectrum + ext);

      const char* spectra_c = spectra[specnum].c_str();
      int speclen = spectra[specnum].size();
      double opd=frlr.optical_path_difference;
      double spec_spacing = frlr.spacing_raw_spectrum;
      
      //adjust spec_spacing
      double center_f= 
        (Spectral_range.value(specnum, 0) + Spectral_range.value(specnum, 1)) / 2;

      spec_spacing *=  riair(&frlr.laser_frequency,
                             &frlr.inside_temperature,
                             &frlr.inside_pressure, &frlr.inside_humidity)/
        
        riair(&center_f,&frlr.inside_temperature,
              &frlr.inside_pressure,&frlr.inside_humidity);
      int ifirst = frlr.index_first;
      int ilast  = frlr.index_last;
      int possp=frlr.length_attached_header;
      int bytepw= frlr.byte_per_word;
      int interp=1;
      float  foff=0;
      float res=0; 
      int mii=103847;
      float slit[mii];
      
      // Not really sure what ilscycle is suppose to be, but the
      // documentation in the old Fortran code seem to require this to
      // be "25". For now, just hardcode this
      // int nscycle= ilscycle[specnum];
      const int nscycle= 25;

      // This definition matches the one used internally by jetspe for calculating the space
      // for reading in data
      int mmp = (Spectral_range.value(specnum, 1) - Spectral_range.value(specnum, 0)) / spec_spacing
        + int(nscycle/(opd*spec_spacing)) + 2;

      Array<float, 1> rad_out(mmp);
      int nmp = 0;
      double nustrt;
      double delwav;
      int status;
      
      int apodization_method = 1;
      double nus = Spectral_range.value(specnum, 0);
      double nue = Spectral_range.value(specnum, 1);

      jetspe(spectra_c, &speclen, &opd, &spec_spacing, &ifirst, &ilast, &possp, 
             &bytepw, &nus, &nue, &apodization_method, &interp, &foff, &res, slit, &mii, 
             &nscycle, rad_out.dataFirst(), &mmp, &nmp, &nustrt, &delwav, &status);
      
      if(status != 0)
        throw Exception("Call to jetspe failed");
      
      Array<double, 1> rad_mem(cast<double>(rad_out(Range(0, nmp-1))));
      radiance_arr.reference(rad_mem);
      uncert_arr.resize(radiance_arr.shape());
      uncert_arr = max(radiance_arr) / frlr.snr;

      freq_spacing(specnum) = delwav;
      freq_start(specnum)   = nustrt;
      freq_end(specnum)     = nustrt + (nmp - 1) * delwav;
    }

    // Add something to the vector even if its empty...
    SpectralRange spec_range(radiance_arr, Unit("W / m^2 / sr / cm^-1"), uncert_arr);
    radiances_.push_back(spec_range);

    run_log_.push_back(frlr);

  } // end of for loop over spectometer

}

//-----------------------------------------------------------------------
/// Creates a L1B FTS object from an HDF file in the format of a L2
/// output file
//-----------------------------------------------------------------------

Level1bFts::Level1bFts(const HdfFile& Hfile,
                       const std::vector<std::string>& Band_names,
                       const std::string& Radiance_dataset)
{
  // Read runlog from HDF file, assume group is FtsRunLog
  // since other parts of this reader also assume L2 HDF group names
  FtsRunLog rl_hdf(Hfile, "FtsRunLog", Band_names);

  int num_spec = (int) Band_names.size(); 
  freq_start.resize(num_spec);
  freq_end.resize(num_spec);
  freq_spacing.resize(num_spec);
  
  freq_spacing = 0;
  freq_start   = 0;
  freq_end     = 0;
  
  Range all = Range::all();
  Array<int, 1> num_colors = Hfile.read_field<int, 2>("/SpectralParameters/num_colors_per_band")(0, all); 
  Array<double, 1> radiance_all = Hfile.read_field<double, 2>(Radiance_dataset)(0, all);

  int rad_idx = 0; // Index of where in all radiance field we are
  for(int specnum = 0; specnum < num_spec; ++specnum) {
    bool non_empty_window = num_colors(specnum) >  0;

    Array<double, 1> radiance_arr;
    Array<double, 1> uncert_arr;
    FtsRunLogRecord frlr;

    // Only process spectrometers that have non-empty ranges
    if(non_empty_window) {
      int n_band_colors = num_colors(specnum);
      frlr = rl_hdf.read(Band_names[specnum]);

      // Copy radiance out of HDF file
      radiance_arr.reference( radiance_all(Range(rad_idx, rad_idx+n_band_colors-1)) );
     
      // Calculate uncertainty from radiance and runlog record
      uncert_arr.resize(radiance_arr.shape());
      uncert_arr = max(radiance_arr) / frlr.snr;

      // Figure out frequency spacing, offset, ending
      double disp_offset = Hfile.read_field<double, 1>("/RetrievalResults/dispersion_offset_" + Band_names[specnum])(0);
      double disp_spacing = Hfile.read_field<double, 1>("/RetrievalResults/dispersion_spacing_" + Band_names[specnum])(0);

      freq_start(specnum) = disp_offset;
      freq_end(specnum) = disp_offset + (n_band_colors-1)*disp_spacing; 
      freq_spacing(specnum) = disp_spacing;

      rad_idx += n_band_colors;
    }

    // Add something to the vector even if its empty...
    SpectralRange spec_range(radiance_arr, Unit("W / m^2 / sr / cm^-1"), uncert_arr);
    radiances_.push_back(spec_range);

    run_log_.push_back(frlr);
  } 
}
