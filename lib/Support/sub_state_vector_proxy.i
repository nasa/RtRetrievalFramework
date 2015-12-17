// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "sub_state_vector_proxy.h"
%}

%import "state_vector.i"

%fp_shared_ptr(FullPhysics::SubStateVectorProxy);

namespace FullPhysics {
class SubStateVectorProxy : public SubStateVectorObserver {
    public:
        virtual ~SubStateVectorProxy();
        virtual void update_sub_state(const ArrayAd<double, 1>& Sv_sub, const blitz::Array<double, 2>& Cov_sub);
        virtual void mark_used_sub(blitz::Array<bool, 1>& Used) const;
        virtual void state_vector_name_sub(blitz::Array<std::string, 1>& Sv_name) const;
        virtual void notify_add(StateVector& Sv);
        virtual void notify_remove(StateVector& Sv);
        virtual void print(std::ostream& Os) const = 0;
    protected:
        SubStateVectorProxy();
        void initialize(const std::vector<boost::shared_ptr<SubStateVectorObserver> >& Proxied);
};
}
