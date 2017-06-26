// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "output_hdf.h"
%}
%base_import(output)
%import "hdf_file_generating.i"
%fp_shared_ptr(FullPhysics::OutputHdf);

namespace FullPhysics {
class OutputHdf : public Output {
public:
  OutputHdf(const std::string& Fname, int Num_level, int Statevector_size, 
	    int Num_aerosol, int Num_band);
  OutputHdf(const boost::shared_ptr<FullPhysics::HdfFileGenerating>& H, int Num_level, 
	    int Statevector_size, int Num_aerosol, int Num_band);
};
}
