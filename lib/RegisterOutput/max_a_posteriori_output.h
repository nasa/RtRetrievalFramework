#ifndef MAX_A_POSTERIORI_OUTPUT_H
#define MAX_A_POSTERIORI_OUTPUT_H
#include "register_output_base.h"
#include "max_a_posteriori.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the MaxAPosteriori class that should be
  written as output.

  See the discussion in RegisterOutputBase why this isn't just part of
  the ConnnorSolver class.
*******************************************************************/
class MaxAPosterioriOutput : public RegisterOutputBase {
public:
//-----------------------------------------------------------------------
/// You can optionally turn on the writing of the Jacobian. We don't
/// do this by default because the files get large if we do, and most
/// of the time we don't need the Jacobian values.
//-----------------------------------------------------------------------

  MaxAPosterioriOutput(const boost::shared_ptr<MaxAPosteriori>& MAP, 
                       bool Write_jacobian = false)
    : map(MAP), write_jacobian(Write_jacobian) {}
  virtual ~MaxAPosterioriOutput() {}
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
private:
  blitz::Array<double,2> a_posteriori_covariance() const
  { return map->a_posteriori_covariance(); }

  blitz::Array<double,2> a_priori_cov() const
  { return map->a_priori_cov(); }

  blitz::Array<double,2> averaging_kernel() const
  { return map->averaging_kernel(); }

  int parameter_size() const
  { return map->parameter_size(); }

  blitz::Array<double,1> a_priori_params() const
  { return map->a_priori_params(); }

  blitz::Array<double,1> param_a_priori_uncertainty() const
  { return map->param_a_priori_uncertainty(); }

  blitz::Array<double,1> parameters() const
  { return map->parameters(); }

  blitz::Array<double,1> param_a_posteriori_uncertainty() const
  { return map->param_a_posteriori_uncertainty(); }

  blitz::Array<double,2> jacobian() const
  { return map->jacobian(); }

  boost::shared_ptr<MaxAPosteriori> map;
  bool write_jacobian;
};
}
#endif
