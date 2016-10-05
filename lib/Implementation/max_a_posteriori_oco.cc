#include <max_a_posteriori_oco.h>


using namespace FullPhysics;
using namespace blitz;



#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(MaxAPosterioriOCO, MaxAPosteriori)
.def(luabind::constructor< const boost::shared_ptr<ForwardModel>&,
                           const Array<double, 1>,
                           const Array<double, 2> >())
REGISTER_LUA_END()
#endif



MaxAPosterioriOCO::MaxAPosterioriOCO(const boost::shared_ptr<ForwardModel>& fm,
                                     const Array<double, 1> a_priori_params,
                                     const Array<double, 2> a_priori_cov)
  : ModelMeasure(
      fm->measured_radiance_all().spectral_range().data(),
      Array<double, 1>(sqr(fm->measured_radiance_all().spectral_range().uncertainty()))), 
    MaxAPosteriori(a_priori_params, a_priori_cov),
    ModelMeasureOCO(fm) 
{
  if(Xa.rows() != fm->state_vector()->observer_claimed_size())
    throw Exception("A priori state vector size and state vector size expected by the model are not equal. :( ");
}


//  TEMPORARY
//
// Should go away after we end support for 
// fixed pressure level grid. 
// Not implemented efficiently.
#include <linear_algebra.h>
void MaxAPosterioriOCO::vanishing_params_update()
{
  ModelMeasureOCO::vanishing_params_update();
  Array<bool, 1> used(FM->state_vector()->used_flag());
  if(used.rows() != Sa_chol.rows())
    throw Exception("Size of a-priori cov. matrix and the number of the elements of used-flag inconsistent! :( ");
  if(!all(used)) {
//    throw Exception("Handling vanishing parameters is not supported yet! :( ");
    Sa_chol.reference(cholesky_decomposition(Sa));
    for(int i=0; i<used.rows(); i++)
      if(!used(i))
        Sa_chol(i,Range::all()) = 0.0;
    Sa_chol_inv.reference(generalized_inverse(Sa_chol,1e-20));
    // Theoretically the selected columns of Sa_chol_inv should
    // be zero; however, they may have very small nonzero values.
    for(int i=0; i<used.rows(); i++)
      if(!used(i))
        Sa_chol_inv(Range::all(),i) = 0.0;
  }
}
