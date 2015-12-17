#ifndef FE_DISABLE_EXCEPTION_H
#define FE_DISABLE_EXCEPTION_H
#include <fenv.h>

namespace FullPhysics {
/****************************************************************//**
  To detect things like divide by zero, we may turn on floating point
  exceptions.  

  However, there may be code that legitimately cause floating point
  exceptions. This might be due to copy of garbage value which by
  chance cause an over or under flow (Lua does this). These garbage
  values aren't actually a problem unless they are used somewhere.

  To allow code like this to execute, this class turns off floating
  point exceptions and then turn them back on in a destructor. So you
  can create an object of this type in any function where we don't
  want to trigger floating point exceptions.
*******************************************************************/

class FeDisableException {
public:
  FeDisableException()
  {
    // Mac doesn't have this function, even though it is a C99
    // function. We check for this during configuration.
#ifdef HAVE_FEENABLEEXCEPT
    fegetexceptflag(&initial_flag, FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
    initial_trap = fedisableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif    
  }
  ~FeDisableException() 
  {
#ifdef HAVE_FEENABLEEXCEPT
    feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
    fesetexceptflag(&initial_flag, FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
    feenableexcept(initial_trap);
#endif
  }
private:
  fexcept_t initial_flag;
  int initial_trap;
};

}
#endif
