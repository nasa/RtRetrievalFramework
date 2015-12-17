#include "max_a_posteriori_output.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(MaxAPosterioriOutput, RegisterOutputBase)
.def(luabind::constructor<const boost::shared_ptr<MaxAPosteriori>& >())
.def(luabind::constructor<const boost::shared_ptr<MaxAPosteriori>&, bool>())
REGISTER_LUA_END()
#endif

void MaxAPosterioriOutput::register_output(const boost::shared_ptr<Output>& out) const
{
  { boost::function<blitz::Array<double,2> ()> f = boost::bind(&MaxAPosterioriOutput::a_posteriori_covariance, this);
  out->register_data_source("/RetrievalResults/aposteriori_covariance_matrix", f); }

  { boost::function<blitz::Array<double,2> ()> f = boost::bind(&MaxAPosterioriOutput::a_priori_cov, this);
  out->register_data_source("/RetrievalResults/apriori_covariance_matrix", f); }

  { boost::function<blitz::Array<double,2> ()> f = boost::bind(&MaxAPosterioriOutput::averaging_kernel, this);
  out->register_data_source("/RetrievalResults/averaging_kernel_matrix", f); }

  { boost::function<int ()> f = boost::bind(&MaxAPosterioriOutput::parameter_size, this);
  out->register_data_source("/RetrievalResults/num_state_vector_elements", f); }


  { boost::function<blitz::Array<double,1> ()> f = boost::bind(&MaxAPosterioriOutput::a_priori_params, this);
  out->register_data_source("/RetrievedStateVector/state_vector_apriori", f); }

  { boost::function<blitz::Array<double,1> ()> f = boost::bind(&MaxAPosterioriOutput::param_a_priori_uncertainty, this);
  out->register_data_source("/RetrievedStateVector/state_vector_apriori_uncert", f); }

  //  The following two may need to go somewhere else.
  //  The two outputs must have zeros for the "unused"
  //  parameters.
  { boost::function<blitz::Array<double,1> ()> f = boost::bind(&MaxAPosterioriOutput::parameters, this);
  out->register_data_source("/RetrievedStateVector/state_vector_result", f); }

  { boost::function<blitz::Array<double,1> ()> f = boost::bind(&MaxAPosterioriOutput::param_a_posteriori_uncertainty, this);
  out->register_data_source("/RetrievedStateVector/state_vector_aposteriori_uncert", f); }

  if(write_jacobian)
    { boost::function<blitz::Array<double,2> ()> f = boost::bind(&MaxAPosterioriOutput::jacobian, this);
    out->register_data_source("/RetrievalResults/jacobian", f); }
}
