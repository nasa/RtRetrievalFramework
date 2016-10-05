#include "initial_guess.h"

using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(InitialGuess)
.def("initial_guess", &InitialGuess::initial_guess)
.def("apriori", &InitialGuess::apriori)
.def("apriori_covariance", &InitialGuess::apriori_covariance)
REGISTER_LUA_END()
#endif
