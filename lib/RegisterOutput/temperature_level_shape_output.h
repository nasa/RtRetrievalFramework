#ifndef TEMPERATURE_LEVEL_SHAPE_OUTPUT_H
#define TEMPERATURE_LEVEL_SHAPE_OUTPUT_H
#include "register_output_base.h"
#include "temperature_level_shape.h"

namespace FullPhysics {
    /****************************************************************//**
      This registers the portions of the TemperatureLevelShape class that
      should be written as output.

      See the discussion in RegisterOutputBase why this isn't just part of
      the TemperatureLevelShape class.
    *******************************************************************/
    class TemperatureLevelShapeOutput : public RegisterOutputBase {
    public:
        TemperatureLevelShapeOutput(const boost::shared_ptr<TemperatureLevelShape>& T)
            : t(T) {}
        virtual ~TemperatureLevelShapeOutput() {}
        virtual void register_output(const boost::shared_ptr<Output>& out) const;
        virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
    private:
        boost::shared_ptr<TemperatureLevelShape> t;
    };
}
#endif
