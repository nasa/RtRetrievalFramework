#ifndef OUTPUT_HDF_ITERATION_H
#define OUTPUT_HDF_ITERATION_H
#include "output.h"
#include "hdf_file_generating.h"
#include "connor_solver.h"

namespace FullPhysics {
// Helper class for OutputHdfIteration
// Don't have Doxygen document this class
/// @cond
class OutputHdfIterationHelperBase {
public:
  virtual ~OutputHdfIterationHelperBase() {}
  virtual void write() = 0;
};

/// @endcond

/****************************************************************//**
  This write the output of the Level 2 Full physics. This particular
  implementation writes out a HDF5 file. 

  This class handles writing iteration output. Note that this *only*
  handles the iteration output, normally you would combine this with
  OutputHdf to get both the iteration output and final results.
*******************************************************************/

class OutputHdfIteration : public OutputTemplate<OutputHdfIteration>,
		           public Observer<ConnorSolver> {
public:
  OutputHdfIteration(const boost::shared_ptr<HdfFileGenerating>& H,
		     const std::string& Iteration_group = "Iteration")
    :h(H), iteration_group(Iteration_group) {}
  virtual ~OutputHdfIteration() { close(); }
  void close();
  virtual void print(std::ostream& Os) const {Os << "OutputHdfIteration";}
  virtual void notify_update(const ConnorSolver& Solver)
  { write(); }
protected:
  virtual void end_because_of_error() { h->abandon(); }
  template<class T> void write_data_t(const std::string& Dataset_name,
				      T Val);
  template<class T, int D> void write_data_t(const std::string& Dataset_name,
					     const blitz::Array<T, D>& Val);
  friend class OutputTemplate<OutputHdfIteration>; 
private:
  boost::shared_ptr<HdfFileGenerating> h;
  std::map<std::string, boost::shared_ptr<OutputHdfIterationHelperBase> >
  data;
  std::string iteration_group;
};
}
#endif
