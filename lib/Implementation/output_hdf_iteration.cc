#include "output_hdf_iteration.h"

using namespace FullPhysics;
using namespace blitz;

// Helper class for OutputHdfIteration
// Don't have Doxygen document this class
/// @cond
namespace FullPhysics {
template<class T, int D> class OutputHdfIterationHelper:
    public OutputHdfIterationHelperBase {
public:
  OutputHdfIterationHelper(const boost::shared_ptr<HdfFileGenerating> H,
                           const std::string& Data_name)
    : h(H), d(Data_name) {}
  // This tacks dt to the end of data.
  void add_iteration(const Array<T, D-1>& dt)
  {
    TinyVector<int, D+1> sz;
    bool incoming_smaller = false;
    if(data.cols() == 0) {
      for(int i = 0; i < D - 1; ++i)
        sz[i + 2] = dt.shape()[i];
      sz[0] = 1;
      sz[1] = 1;
    } else {
      for(int i = 0; i < D - 1; ++i) {
        // Grow data array if incoming data is bigger in 
        // a dimension
        if(dt.shape()[i] > data.shape()[i + 2]) {
          sz[i+2] = dt.shape()[i];
        } else if(dt.shape()[i] < data.shape()[i + 2]) {
          incoming_smaller = true;
          sz[i+2] = data.shape()[i+2];
        } else {
          sz[i+2] = data.shape()[i+2];
        }
      }
      sz[0] = 1;
      sz[1] = data.cols() + 1;
    }
    data.resizeAndPreserve(sz);

    if (incoming_smaller) {
      // Match data portion that is only as
      // big as incoming data
      TinyVector<int, D-1> lower(0);
      RectDomain<D-1> dt_subdomain(lower, dt.shape());

      // This subdomain use WILL die for string types
      // but hopefully unlikely they chagne size during execution lifetime
      subset(data) = 0; // set all entries to 0
      subset(data)(dt_subdomain) = dt;
    } else {
      subset(data) = dt;
    }
  }
  virtual void write() { h->hdf_file().write_field(d, data); }
private:
  Array<T, 1> subset(const Array<T, 3>& V) 
  {return V(0, V.cols() - 1, Range::all());}
  Array<T, 2> subset(const Array<T, 4>& V) 
  {return V(0, V.cols() - 1, Range::all(), Range::all());}
  Array<T, 3> subset(const Array<T, 5>& V) 
  {return V(0, V.cols() - 1, Range::all(), Range::all(), Range::all());}
  boost::shared_ptr<HdfFileGenerating> h;
  std::string d;
  Array<T, D+1> data;
};

template<class T> class OutputHdfIterationHelper<T, 1> : 
    public OutputHdfIterationHelperBase {
public:
  OutputHdfIterationHelper(const boost::shared_ptr<HdfFileGenerating> H,
			   const std::string& Data_name)
    : h(H), d(Data_name) {}
  void add_iteration(const T& dt)
  {
    data.resizeAndPreserve(1, data.cols() + 1);
    data(0, data.cols() - 1) = dt;
  }
  virtual void write() { h->hdf_file().write_field(d, data); }
private:
  boost::shared_ptr<HdfFileGenerating> h;
  std::string d;
  Array<T, 2> data;
};
}
/// @endcond

void OutputHdfIteration::close()
{
  typedef std::map<std::string, 
		   boost::shared_ptr<OutputHdfIterationHelperBase> >::value_type
    vt;
  BOOST_FOREACH(vt& p, data)
    p.second->write();
  h.reset();
}

template<class T> void OutputHdfIteration::write_data_t
(const std::string& Dataset_name,
 T Val)
{ 
  if(data.count(Dataset_name) == 0)
    data[Dataset_name] = boost::shared_ptr<OutputHdfIterationHelperBase>
      (new OutputHdfIterationHelper<T, 1>
       (h, iteration_group + "/" + Dataset_name));
  boost::shared_ptr<OutputHdfIterationHelper<T, 1> >
    hlp = boost::dynamic_pointer_cast<OutputHdfIterationHelper<T, 1> >
    (data[Dataset_name]);
  hlp->add_iteration(Val);
}

template<class T, int D> void OutputHdfIteration::write_data_t
(const std::string& Dataset_name,
 const blitz::Array<T, D>& Val)
{ 
  if(data.count(Dataset_name) == 0)
    data[Dataset_name] = boost::shared_ptr<OutputHdfIterationHelperBase>
      (new OutputHdfIterationHelper<T, D + 1>
       (h, iteration_group + "/" + Dataset_name));
  boost::shared_ptr<OutputHdfIterationHelper<T, D + 1> >
    hlp = boost::dynamic_pointer_cast<OutputHdfIterationHelper<T, D + 1> >
    (data[Dataset_name]);
  hlp->add_iteration(Val);
}

// Instantiation of the templates for the various types.

template void OutputHdfIteration::write_data_t(const std::string& Dataset_name, 
				      int Val);
template void OutputHdfIteration::write_data_t(const std::string& Dataset_name, 
				      int64_t Val);
template void OutputHdfIteration::write_data_t(const std::string& Dataset_name, 
				      double Val);
template void OutputHdfIteration::write_data_t(const std::string& Dataset_name, 
				      std::string Val);
template void OutputHdfIteration::write_data_t(const std::string& Dataset_name, 
				      const char* Val);

template void OutputHdfIteration::write_data_t(const std::string& Dataset_name, 
				      const blitz::Array<int, 1>& Val);
template void OutputHdfIteration::write_data_t(const std::string& Dataset_name, 
				      const blitz::Array<std::string, 1>& Val);
template void OutputHdfIteration::write_data_t(const std::string& Dataset_name, 
				      const blitz::Array<const char*, 1>& Val);
template void OutputHdfIteration::write_data_t(const std::string& Dataset_name, 
				      const blitz::Array<double, 1>& Val);

template void OutputHdfIteration::write_data_t(const std::string& Dataset_name, 
				      const blitz::Array<int, 2>& Val);
template void OutputHdfIteration::write_data_t(const std::string& Dataset_name, 
				      const blitz::Array<std::string, 2>& Val);
template void OutputHdfIteration::write_data_t(const std::string& Dataset_name, 
				      const blitz::Array<const char*, 2>& Val);
template void OutputHdfIteration::write_data_t(const std::string& Dataset_name, 
				      const blitz::Array<double, 2>& Val);

template void OutputHdfIteration::write_data_t(const std::string& Dataset_name, 
				      const blitz::Array<int, 3>& Val);
template void OutputHdfIteration::write_data_t(const std::string& Dataset_name, 
				      const blitz::Array<std::string, 3>& Val);
template void OutputHdfIteration::write_data_t(const std::string& Dataset_name, 
				      const blitz::Array<const char*, 3>& Val);
template void OutputHdfIteration::write_data_t(const std::string& Dataset_name, 
				      const blitz::Array<double, 3>& Val);

