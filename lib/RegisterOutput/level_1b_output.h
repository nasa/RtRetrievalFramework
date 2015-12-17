#ifndef LEVEL_1B_OUTPUT_H
#define LEVEL_1B_OUTPUT_H
#include "register_output_base.h"
#include "level_1b.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the Level1b class that should be
  written as output.

  See the discussion in RegisterOutputBase why this isn't just part of
  the Level_1b class.
*******************************************************************/
class Level1bOutput : public RegisterOutputBase {
public:
  Level1bOutput(const boost::shared_ptr<Level1b>& F) : f(F) {}
  virtual ~Level1bOutput() {}
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
private:
  boost::shared_ptr<Level1b> f;
};
}
#endif
