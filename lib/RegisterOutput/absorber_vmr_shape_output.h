#ifndef ABSORBER_VMR_SHAPE_OUTPUT_H
#define ABSORBER_VMR_SHAPE_OUTPUT_H
#include "register_output_base.h"
#include "absorber_vmr_shape.h"
#include "state_vector.h"

namespace FullPhysics {
    /****************************************************************//**
      This registers the portions of the AbsorberVmrShape class that
      should be written as output.

      See the discussion in RegisterOutputBase why this isn't just part of
      the AbsorberVmrShape class.
    *******************************************************************/
    class AbsorberVmrShapeOutput : public RegisterOutputBase {
    public:
        AbsorberVmrShapeOutput(const boost::shared_ptr<AbsorberVmrShape>& A)
            : a(A) {}
        virtual ~AbsorberVmrShapeOutput() {}
        virtual void register_output(const boost::shared_ptr<Output>& out) const;
        virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
    private:
        boost::shared_ptr<AbsorberVmrShape> a;
    };
}
#endif
