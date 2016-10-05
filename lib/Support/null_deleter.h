#ifndef NULL_DELETER_H
#define NULL_DELETER_H
#include <boost/shared_ptr.hpp>

namespace FullPhysics {
/****************************************************************//**
  Class for use with boost::shared_ptr when we don't want data to
  be deleted.
*******************************************************************/

class null_deleter
{
public:
    void operator()(void const *) const
    {
    }
};

//-----------------------------------------------------------------------
/// Helper routine to get a shared_ptr from a reference. The reference
/// is not deleted when the shared_ptr is.
//-----------------------------------------------------------------------

template<class T> boost::shared_ptr<T> to_ptr(T& t)
{ return boost::shared_ptr<T>(&t, null_deleter());}
}
#endif
