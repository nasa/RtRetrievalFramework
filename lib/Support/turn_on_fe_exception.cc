#include "turn_on_fe_exception.h"
#include "fp_exception.h"
#include <signal.h>
#include <fenv.h>
using namespace FullPhysics;

// This is a handler in case we get a floating point exception.
void floating_point_exception2(int signum, siginfo_t *sip, void *scp)
{
  Exception e;
  e << "Exception thrown by Full Physics code:\n"
            << "   Floating point exception\n"
	    << "   Address exception occurred: " << sip->si_addr << "\n"
	    << "   Exception type:\n" << "      ";
  int fe_code = sip->si_code;
  switch(fe_code) {
  case FPE_INTDIV:
    e << "Integer divide by zero\n";
    break;
  case FPE_INTOVF:
    e << "Integer overflow\n";
    break;
  case FPE_FLTDIV:
    e << "Floating point divide by zero\n";
    break;
  case FPE_FLTOVF:
    e << "Floating point overflow\n";
    break;
  case FPE_FLTUND:
    e << "Floating point underflow\n";
    break;
  case FPE_FLTRES:
    e << "Floating point inexact result\n";
    break;
  case FPE_FLTINV:
    e << "Floating point invalid operation\n";
    break;
  case FPE_FLTSUB:
    e << "Subscript out of range\n";
    break;
  default:
    e << "Unknown floating point exception.\n";
  }
  throw e;
}

//-----------------------------------------------------------------------
/// It can be useful during unit test/debugging/python usage to be
/// able to have any floating point exception throw an error. This
/// turns that on.
///
/// Note that the top level executable l2_fp already handles this, you
/// don't want to call this function from there.
//-----------------------------------------------------------------------

void FullPhysics::turn_on_fe_exception()
{
  feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
  struct sigaction act;
  act.sa_sigaction = floating_point_exception2;
  sigemptyset(&act.sa_mask);
  act.sa_flags = SA_SIGINFO;    
  sigaction(SIGFPE, &act, (struct sigaction *)0);
}
