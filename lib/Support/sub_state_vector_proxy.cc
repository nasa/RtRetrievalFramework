#include "sub_state_vector_proxy.h"
#include "ostream_pad.h"

using namespace FullPhysics;
using namespace blitz;

//-----------------------------------------------------------------------
/// Registers the classes that will be proxied. Must be called from
/// subclasses constructor or else state vector will not be set up
/// correctly.
//-----------------------------------------------------------------------

void SubStateVectorProxy::initialize(const std::vector<boost::shared_ptr<SubStateVectorObserver> >& Proxied)
{
    proxied_observers = Proxied;

    int plen = 0;
    BOOST_FOREACH(boost::shared_ptr<SubStateVectorObserver> curr_obs, proxied_observers) {
        plen += curr_obs->sub_vector_size();
    }
    state_vector_observer_initialize(plen);
}

//-----------------------------------------------------------------------
/// Extracts the relevant portions of the passed arrays and passes
/// them to the proxied objects
//-----------------------------------------------------------------------

void SubStateVectorProxy::update_sub_state(const ArrayAd<double, 1>& Sv_sub, const blitz::Array<double, 2>& Cov_sub)
{
    int offset = 0;
    BOOST_FOREACH(boost::shared_ptr<SubStateVectorObserver> curr_obs, proxied_observers) {
        Range prox_range(offset, offset + curr_obs->sub_vector_size() - 1);
        curr_obs->update_sub_state(Sv_sub(prox_range), Cov_sub(prox_range, prox_range) );
        offset += curr_obs->sub_vector_size();
    }
}

//-----------------------------------------------------------------------
/// Extracts the relevant portions of the passed arrays and passes
/// them to the proxied objects
//-----------------------------------------------------------------------

void SubStateVectorProxy::mark_used_sub(blitz::Array<bool, 1>& Used) const
{
    int offset = 0;
    BOOST_FOREACH(boost::shared_ptr<SubStateVectorObserver> curr_obs, proxied_observers) {
        Range prox_range(offset, offset + curr_obs->sub_vector_size() - 1);
        Array<bool, 1> prox_used(Used(prox_range));
        curr_obs->mark_used_sub(prox_used);
        offset += curr_obs->sub_vector_size();
    }
}

//-----------------------------------------------------------------------
/// Extracts the relevant portions of the passed arrays and passes
/// them to the proxied objects
//-----------------------------------------------------------------------


void SubStateVectorProxy::state_vector_name_sub(blitz::Array<std::string, 1>& Sv_name) const
{
    int offset = 0;
    BOOST_FOREACH(boost::shared_ptr<SubStateVectorObserver> curr_obs, proxied_observers) {
        Range prox_range(offset, offset + curr_obs->sub_vector_size() - 1);
        blitz::Array<std::string, 1> prox_name(Sv_name(prox_range));
        curr_obs->state_vector_name_sub(prox_name);
        offset += curr_obs->sub_vector_size();
    }
}

//-----------------------------------------------------------------------
/// Output the print results from the proxied classes
//-----------------------------------------------------------------------

void SubStateVectorProxy::print(std::ostream& Os) const
{
    OstreamPad opad(Os, "    ");
    Os << "SubStateVectorProxy";
    BOOST_FOREACH(boost::shared_ptr<SubStateVectorObserver> curr_obs, proxied_observers) {
        curr_obs->print(opad);
    }
    opad.strict_sync();
}
