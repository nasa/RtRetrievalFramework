#include "output_hdf.h"
#include "fp_exception.h"
#include <boost/regex.hpp>
#include <cstdio>

using namespace FullPhysics;
using namespace blitz;

//-----------------------------------------------------------------------
/// Constructor. This takes the file name to write.
//-----------------------------------------------------------------------

OutputHdf::OutputHdf(const std::string& Fname, int Num_level, 
		     int Statevector_size, int Num_aerosol, 
		     int Num_band) 
: num_level(Num_level), statevector_size(Statevector_size),
  num_aerosol(Num_aerosol), num_band(Num_band)
{
  h.reset(new HdfFileGenerating(Fname));
  initialize();
}

//-----------------------------------------------------------------------
/// Constructor. This takes the file to write to.
//-----------------------------------------------------------------------

OutputHdf::OutputHdf(const boost::shared_ptr<HdfFileGenerating>& H, 
		     int Num_level, 
		     int Statevector_size, int Num_aerosol, 
		     int Num_band) 
: num_level(Num_level), statevector_size(Statevector_size),
  num_aerosol(Num_aerosol), num_band(Num_band), h(H)
{
  initialize();
}

void OutputHdf::initialize()
{
  range_min_check(num_level, 1);
  range_min_check(num_aerosol, 1);

  // See discussion in output_hdf.h for explanation of this check.

  if(num_level == statevector_size ||
     num_aerosol == statevector_size ||
     num_aerosol == num_level) {
    Exception e;
    e << "The current implementation assumes Num_level, Statevector_size\n"
      << "and Num_aerosol are all different. This is entirely to simplify\n"
      << "the implementation, see the discussion in output_hdf.h. If\n"
      << "needed, we could relax this constraint by changing the code.\n"
      << "The values of Num_level = " << num_level << ", Statevector_size = "
      << statevector_size << ",\n"
      << "and Num_aerosol = " << num_aerosol << " violates this assumption.";
    throw e;
  }
}

// See base class for descriptions of all these functions.

void OutputHdf::start_write()
{
  write_dimension_metadata();
  write_shape_metadata();
}

void OutputHdf::end_because_of_error()
{
  h->abandon();
}

//-----------------------------------------------------------------------
/// Add the dimension metadata to the file.
///
/// Right now, there are only a handful of dimensions. We go ahead and
/// hard code this, since this is the easiest implementation. If this
/// list grows in complexity, we may want to come up with a more
/// flexible (but also more complicated) implementation.
//-----------------------------------------------------------------------

void OutputHdf::write_dimension_metadata()
{
  h->hdf_file().dimension_metadata("Retrieval", 
       "Number retrievals reported.", 1);
  h->hdf_file().dimension_metadata("StateVectorElement", 
       "Retrieved state vector elements.", statevector_size);
  h->hdf_file().dimension_metadata("Level",
       "Atmospheric retrieval levels.", num_level);
  h->hdf_file().dimension_metadata("Aerosol",
       "Retrieved aerosol type.", num_aerosol);
  h->hdf_file().dimension_metadata("Band",
       "Spectral band.", num_band);
}

//-----------------------------------------------------------------------
/// Add the shape metadata to the file.
///
/// Right now, there are only a handful of shapes. We go ahead and
/// hard code this, since this is the easiest implementation. If this
/// list grows in complexity, we may want to come up with a more
/// flexible (but also more complicated) implementation.
//-----------------------------------------------------------------------

void OutputHdf::write_shape_metadata()
{
  h->hdf_file().shape_metadata("Retrieval_Array", "Retrieval");
  h->hdf_file().shape_metadata(
     "Retrieval_StateVectorElement_StateVectorElement_Array",
     "Retrieval", "StateVectorElement", "StateVectorElement");
  h->hdf_file().shape_metadata("Retrieval_Level_Array", "Retrieval", "Level");
  h->hdf_file().shape_metadata("Retrieval_StateVectorElement_Array",
			   "Retrieval", "StateVectorElement");
  h->hdf_file().shape_metadata("Retrieval_Aerosol_Array", "Retrieval", "Aerosol");
  h->hdf_file().shape_metadata("Retrieval_Band_Array", "Retrieval", "Band");
}

template<class T> void OutputHdf::write_data_t(
    const std::string& Dataset_name, T Val)
{
  h->hdf_file().write_field(Dataset_name, Val);
  // Metadata fields are scalars, all others are indexed by Retrieval.
  if(Dataset_name.find("/Metadata/") == std::string::npos)
    h->hdf_file().write_attribute(Dataset_name + "/Shape", "Retrieval_Array");
  else
    h->hdf_file().write_attribute(Dataset_name + "/Shape", "Scalar");
}

template<class T> void OutputHdf::write_data_t(
    const std::string& Dataset_name, const blitz::Array<T, 1>& Val)
{
  using namespace blitz;
  // We add an extra retrieval dimension, even though it has only 1
  // value. This makes this look like the aggregated product SDOS
  // produces. 
  Array<T, 2> val2(1, Val.rows());
  val2(0, Range::all()) = Val;
  h->hdf_file().write_field(Dataset_name, val2);

  // Check size, and use this to fill in Shape metadata.  We silently
  // leave off Shape metadata if don't recognize the size. See discussion
  // of this design decision in output_hdf.h
  
  if(Val.rows() == statevector_size)
    h->hdf_file().write_attribute(Dataset_name + "/Shape", 
			      "Retrieval_StateVectorElement_Array");
  else if(Val.rows() == num_level)
    h->hdf_file().write_attribute(Dataset_name + "/Shape", 
			      "Retrieval_Level_Array");
  // We can't depend on num_aerosol and num_band being different, so
  // we hardcode the few fields we have that depend on band.
  else if(boost::regex_match(Dataset_name, 
			     boost::regex(".*/num_colors_per_band")) ||
	  boost::regex_match(Dataset_name, boost::regex(".*/Absco.*Scale")))
    h->hdf_file().write_attribute(Dataset_name + "/Shape", 
				  "Retrieval_Band_Array");
  else if(Val.rows() == num_aerosol)
    h->hdf_file().write_attribute(Dataset_name + "/Shape", 
			      "Retrieval_Aerosol_Array");
}

template<class T> void OutputHdf::write_data_t(
    const std::string& Dataset_name, const blitz::Array<T, 2>& Val)
{
  using namespace blitz;
  // We add an extra retrieval dimension, even though it has only 1
  // value. This makes this look like the aggregated product SDOS
  // produces. 
  Array<T, 3> val2(1, Val.rows(), Val.cols());
  val2(0, Range::all(), Range::all()) = Val;
  h->hdf_file().write_field(Dataset_name, val2);

  // Check size, and use this to fill in Shape metadata.  We silently
  // leave off Shape metadata if don't recognize the size. See discussion
  // of this design decision in output_hdf.h
  
  if(Val.rows() == statevector_size &&
     Val.cols() == statevector_size)
    h->hdf_file().write_attribute(Dataset_name + "/Shape", 
	      "Retrieval_StateVectorElement_StateVectorElement_Array");
}

template<class T> void OutputHdf::write_data_t(
    const std::string& Dataset_name, const blitz::Array<T, 3>& Val)
{
  using namespace blitz;
  // We add an extra retrieval dimension, even though it has only 1
  // value. This makes this look like the aggregated product SDOS
  // produces. 
  Array<T, 4> 
    val2(1, Val.rows(), Val.cols(),
	 Val.extent(thirdDim));
  val2(0, Range::all(), Range::all(), Range::all()) = Val;
  h->hdf_file().write_field(Dataset_name, val2);

  // Don't have any hardcoded shapes of rank 3, so silently leave
  // metadata off. See discussion in output_hdf.h for this design
  // decision. 
}

// Instantiation of the templates for the various types.

template void OutputHdf::write_data_t(const std::string& Dataset_name, 
				      int Val);
template void OutputHdf::write_data_t(const std::string& Dataset_name, 
				      int64_t Val);
template void OutputHdf::write_data_t(const std::string& Dataset_name, 
				      double Val);
template void OutputHdf::write_data_t(const std::string& Dataset_name, 
				      const std::string& Val);
template void OutputHdf::write_data_t(const std::string& Dataset_name, 
				      const char* Val);

template void OutputHdf::write_data_t(const std::string& Dataset_name, 
				      const blitz::Array<int, 1>& Val);
template void OutputHdf::write_data_t(const std::string& Dataset_name, 
				      const blitz::Array<std::string, 1>& Val);
template void OutputHdf::write_data_t(const std::string& Dataset_name, 
				      const blitz::Array<const char*, 1>& Val);
template void OutputHdf::write_data_t(const std::string& Dataset_name, 
				      const blitz::Array<double, 1>& Val);

template void OutputHdf::write_data_t(const std::string& Dataset_name, 
				      const blitz::Array<int, 2>& Val);
template void OutputHdf::write_data_t(const std::string& Dataset_name, 
				      const blitz::Array<std::string, 2>& Val);
template void OutputHdf::write_data_t(const std::string& Dataset_name, 
				      const blitz::Array<const char*, 2>& Val);
template void OutputHdf::write_data_t(const std::string& Dataset_name, 
				      const blitz::Array<double, 2>& Val);

template void OutputHdf::write_data_t(const std::string& Dataset_name, 
				      const blitz::Array<int, 3>& Val);
template void OutputHdf::write_data_t(const std::string& Dataset_name, 
				      const blitz::Array<std::string, 3>& Val);
template void OutputHdf::write_data_t(const std::string& Dataset_name, 
				      const blitz::Array<const char*, 3>& Val);
template void OutputHdf::write_data_t(const std::string& Dataset_name, 
				      const blitz::Array<double, 3>& Val);
