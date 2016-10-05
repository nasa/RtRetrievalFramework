//--------------------------------------------------------------
// Add support for boost rational.
//--------------------------------------------------------------

%{
#include <boost/rational.hpp>
%}

namespace boost {
  template<class I>  class rational {
  public:
    rational(I n);
    rational(I n, I d);
  };
}

%template(rational_int) boost::rational<int>;
