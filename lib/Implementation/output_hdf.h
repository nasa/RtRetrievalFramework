#ifndef OUTPUT_HDF_H
#define OUTPUT_HDF_H
#include "output.h"
#include "hdf_file_generating.h"

namespace FullPhysics {
/****************************************************************//**
  This write the output of the Level 2 Full physics. This particular
  implementation writes out a HDF5 file.

  Note that we make a few assumptions to simplify the
  implementation. We could relax any of these constraints, we would
  just need to modify the implementation.

  1. We assume that there are only a few hardcoded Shapes and
     Dimensions. The current implementation just uses a fixed
     hardcoded set. We could create a more general, flexible (but also
     more complicated) implementation if needed in the future.

  2. We determine the shape metadata of any particular field by
     looking at the size. If we don't recognize the size, we simply
     silently leave off the shape metadata information. This allows
     new fields to be added without needing to necessarily update the
     shape information (e.g., a new diagnostic field). A consequence
     of this design decision is that actual mistakes (e.g., wrong
     number of pressures level reported) won't be caught, and new
     field may be missing shape information until we change this
     class. I think this is the right trade, but we may need to
     reevaluate this at some point in the future.

  3. We assume that the different dimensions that can make the same
     rank shape are different sizes (e.g., Retrieval_Level_Array and
     Retrieval_StateVectorElement_Array). It is hard to think of a
     case where this wouldn't be true, but this assumption is entirely
     to make the implementation easier. There is no intrinsic reason
     why we couldn't support these being the same. We check and
     trigger an error if these are the same, but we can change the
     code to support this if needed.

  An alternative implementation would be to have a table mapping each
  field to the shape it is. This isn't hugely complicated, but would
  require us to keep a table of all possible fields. For now, we will
  do the simpler implementation.
*******************************************************************/

class OutputHdf : public OutputTemplate<OutputHdf> {
public:
  OutputHdf(const std::string& Fname, int Num_level, int Statevector_size, 
	    int Num_aerosol, int Number_band);
  OutputHdf(const boost::shared_ptr<HdfFileGenerating>& H, int Num_level, 
	    int Statevector_size, int Num_aerosol, int Number_band);
  virtual ~OutputHdf() {}
  virtual void print(std::ostream& Os) const {Os << "OutputHdf";}
protected:
  virtual void start_write();
  virtual void end_because_of_error();
  template<class T> void write_data_t(const std::string& Dataset_name,
				      T Val);
  template<class T> void write_data_t(const std::string& Dataset_name,
				      const blitz::Array<T, 1>& Val);
  template<class T> void write_data_t(const std::string& Dataset_name,
				      const blitz::Array<T, 2>& Val);
  template<class T> void write_data_t(const std::string& Dataset_name,
				      const blitz::Array<T, 3>& Val);
  friend class OutputTemplate<OutputHdf>; 
private:
  int num_level, statevector_size, num_aerosol, num_band;
  boost::shared_ptr<HdfFileGenerating> h;
  void write_dimension_metadata();
  void write_shape_metadata();
  void initialize();
};
}
#endif
