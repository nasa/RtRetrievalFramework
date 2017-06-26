#include "absorber_absco.h"
#include "fp_exception.h"
#include "linear_algebra.h"
#include "old_constant.h"
#include "ostream_pad.h"
#include "absco.h"
#include <boost/lexical_cast.hpp>
#include <boost/bind.hpp>
using namespace FullPhysics;
using namespace blitz;
inline double sqr(double x) {return x * x;}

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(AbsorberAbsco, Absorber)
.def(luabind::constructor<const std::vector<boost::shared_ptr<AbsorberVmr> >,
			  const boost::shared_ptr<Pressure>&,	   
			  const boost::shared_ptr<Temperature>&,
			  const std::vector<boost::shared_ptr<Altitude> >&,
                          const std::vector<boost::shared_ptr<GasAbsorption> >&,
                          const boost::shared_ptr<Constant>&>())
.def(luabind::constructor<const std::vector<boost::shared_ptr<AbsorberVmr> >,
			  const boost::shared_ptr<Pressure>&,	   
			  const boost::shared_ptr<Temperature>&,
			  const std::vector<boost::shared_ptr<Altitude> >&,
			  const std::vector<boost::shared_ptr<GasAbsorption> >&,
                          const boost::shared_ptr<Constant>&,
			  int>())
REGISTER_LUA_END()
#endif

//----------------------------------------------------------------
/// This is the portion of the optical depth calculation integrand
/// that is independent on the wave number. We separate this out
/// because we can calculate this portion once and cache it.
/// This only changes if the Pressure, Temperature, Altitude, or VMR has
/// changed (e.g., the state vector is updated).
//----------------------------------------------------------------

AutoDerivativeWithUnit<double> AbsorberAbsco::integrand_independent_wn
(int Spec_index, int Species_index, const DoubleWithUnit& P)
const
{
  return vmr_func(Species_index, P) * c->avogadro_constant() /
    (gravity_func(Spec_index, P) * 
     (1 + h2o_vmr_func(P) * c->molar_weight_water() / 
      c->molar_weight_dry_air()) * c->molar_weight_dry_air());
}

//----------------------------------------------------------------
/// Integrand used in the absorption calculation. This is in Pa^-1.
/// wn should be in cm^-1, and P in Pascals.
//----------------------------------------------------------------

double AbsorberAbsco::integrand
(double wn, double P, int Spec_index, int Species_index) const
{
  range_check(Species_index, 0, number_species());
  DoubleWithUnit punit(P, "Pa");
  double v1 = 
    integrand_independent_wn(Spec_index, Species_index, 
			     punit).convert(Unit("cm^-2 / Pa")).value.value();
  AutoDerivativeWithUnit<double> tder = temp_func(punit);
  DoubleWithUnit t(tder.value.value(), tder.units);
  AutoDerivativeWithUnit<double> bvmrder = h2o_vmr_func(punit);
  DoubleWithUnit bvmr(bvmrder.value.value(), bvmrder.units);
  if(gas_absorption[Species_index]->have_data(wn)) {
    DoubleWithUnit cross_sec = gas_absorption[Species_index]->absorption_cross_section(wn, punit, t, bvmr);
    return v1 * cross_sec.value;
  } else
    return 0.0;
}

//----------------------------------------------------------------
/// Find the optical depth for a layer by directly integrating 
/// the integrand function. Note that this is slower than the 
/// optical_depth_each_layer function, and also doesn't calculate
/// the Jacobians. But it is more accurate, and can be used as 
/// a diagnostic tool to investigate the accuracy of 
/// optical_depth_each_layer
//----------------------------------------------------------------

double AbsorberAbsco::optical_depth_each_layer_direct_integrate(
double wn, int Spec_index, int Species_index, int Layer_index,
double eps_abs, double eps_rel) const
{
  range_check(Spec_index, 0, number_spectrometer());
  range_check(Species_index, 0, number_species());
  range_check(Layer_index, 0, press->number_layer());
  if(!gas_absorption[Species_index]->have_data(wn))
    return 0.0;
  boost::function<double (double)> f;
  f = boost::bind(&AbsorberAbsco::integrand, this, wn, _1, Spec_index, 
		  Species_index);
  ArrayAdWithUnit<double, 1> pgrid(press->pressure_grid().convert(units::Pa));
  std::vector<double> bp;
  Array<double, 1> imp = 
    temp->important_pressure_level(). convert(units::Pa).value;
  for(int i = 0; i < imp.rows(); ++i)
    bp.push_back(imp(i));
  return intg.integrate(f, pgrid.value.value()(Layer_index),
			pgrid.value.value()(Layer_index + 1),
			bp,
			eps_abs, eps_rel);
}

//----------------------------------------------------------------
/// Version of optical_depth_each_layer_direct_integrate that returns
/// an array of species/layer like optical_depth_each_layer does.
//----------------------------------------------------------------

blitz::Array<double, 2> 
AbsorberAbsco::optical_depth_each_layer_direct_integrate
(double wn, int Spec_index, double eps_abs, double eps_rel) const
{
  blitz::Array<double, 2> res(press->number_layer(), number_species());
  for(int i = 0; i < res.rows(); ++i)
    for(int j = 0; j < res.cols(); ++j)
      res(i, j) = optical_depth_each_layer_direct_integrate
	(wn, Spec_index, j, i, eps_abs, eps_rel);
  return res;
}

//----------------------------------------------------------------
/// This fills the sublayer structure. This populates pgrid, psub, 
/// weight, and layer range. 
///
/// We set up for doing a Simpsons rule integration, and choose the
/// weights accordingly. There is no assumption that the spacing
/// between points is constant.
///
/// For a layer, we divide the area up into nsub / 2 regions. For 
/// each region, we have 3 equally spaced points that we apply
/// Simpsons rule to. This is:
/// integral from a to b f(x) dx = (b - a) / 6 * (f(a) + 4f((a+b)/2)
/// + f(b))
///
/// Note that f(b) is then the first point used in the next region, 
/// so it gets used twice. We account for this when creating the
/// sublayers, and the weights to apply.
///
/// In addition of the nsub / 2 points, we also add any points marked
/// as "important" by the temperature object. For example, if the
/// temperate is TemperatureMet this is the 91 ECMWF levels that 
/// the temperature is reported on.
//----------------------------------------------------------------

void AbsorberAbsco::create_sublayer() const
{
  firstIndex i1;

  // ****** IMPORTANT **********
  // Note that there is an implicit assumption in tau_gas_der that the end
  // point of the range are at the edge of the integral, so for a
  // layer going from p1 to p2 the edges are p1 and p2. This is the
  // case for Simpson's rule. But if you change this, you may need to
  // modify tau_gas_der


  pgrid.reference(press->pressure_grid().convert(units::Pa));

  // Note that although the pressure grid itself depends on the state vector,
  // the grid psub we create here does *not*. This is because this just the
  // points we happened to pick for doing the integral on. The
  // dependency on pressure actually comes in when we calculate the
  // definite integral integral p1 to p2 (integrand) dp, where we have
  // the derivative with respect to p1 and p2. 

  int npoint = nsub / 2;	// We'll double this when we add
				// midpoints for Simpsons rule.
  Array<double, 1> sublay_fac(npoint);
  sublay_fac = (i1 + 1.0) / (npoint);

  // "Important" pressure levels that should be included in the
  // integral.
  Array<double, 1> pimp_arr = 
    temp->important_pressure_level().convert(units::Pa).value;
  // Might not be sorted, so stick into a list that we will sort.
  std::list<double> pimp;
  for(int i = 0; i < pimp_arr.rows(); ++i)
    pimp.push_back(pimp_arr(i));
  pimp.sort();

  // We set up the pressure grid for the Simpsons rule. For integrating
  // from p1 to p2, we have the points p1, p1 + spacing, ... p2. Note
  // that we reuse p2 in the next interval, because to integrate from
  // p2 to p3 we have p2, p2 + spacing, ... p3. This is taken care of
  // by having overlapping layer_range.
  layer_range.clear();
  std::vector<double> psub_vec;
  psub_vec.push_back(pgrid.value.value()(0));
  // Drop all "important" pressure point that come before the first
  // pressure level.
  while(pimp.size() > 0 && pimp.front() <= psub_vec.front())
    pimp.pop_front();

  // The logic here is a little complicated. In words we do the
  // following:
  // 1. Go through the nsub / 2 equally spaced points between p1 and
  //    p2.
  // 2. When we get a point, first check if there is any "important"
  //    points that come before it. If there are, then add the important
  //    point + the midpoint used by simpsons rule. So 2 points get
  //    added, not just one. The midpoint is the f((a+b)/2) point
  //    mentioned above in Simpsons rule.
  // 3. After dealing with the important points, if any, add in our
  //    equally spaced point, along with the midpoint (so again 2
  //    points).
  // At the same time we are doing this, we keep track of what points
  // fall in a particular layer. We create a Range object that covers
  // this, and save it in layer_range.

  int start = 0;
  for(int i = 0; i < press->number_layer(); ++i) {
    for(int j = 0; j < sublay_fac.rows(); ++j) {
      double pv = pgrid.value.value()(i) * (1 - sublay_fac(j)) +
	pgrid.value.value()(i + 1) * sublay_fac(j);
      while(pimp.size() > 0 && pimp.front() < pv) {
	double midpoint = (psub_vec.back() + pimp.front()) / 2.0;
	psub_vec.push_back(midpoint);
	psub_vec.push_back(pimp.front());
	pimp.pop_front();
      }
      // Unlikely, but possible case where a important pressure is
      // exactly the same as pv. Drop the point if this is the case.
      while(pimp.size() > 0 && pimp.front() == pv)
	pimp.pop_front();
      // Add midpoint, and then the point
      double midpoint = (psub_vec.back() + pv) / 2.0;
      psub_vec.push_back(midpoint);
      psub_vec.push_back(pv);
    }
    int end = ((int) psub_vec.size()) - 1;
    layer_range.push_back(Range(start, end));
    // Start next layer with the last pressure in the current layer.
    start = end;
  }
  // Move points stored in a std::vector into a Array. This is the
  // same data, just in a different structure.
  psub.units = pgrid.units;
  psub.value.resize((int) psub_vec.size());
  for(int i = 0; i < psub.rows(); ++i)
    psub.value(i) = psub_vec[i];

  // Now go through and figure out the weights for all our
  // points. This is just the weights for Simpson's rule, but with the
  // complication that the spacing between regions we apply Simpson's
  // rule may change.
  weight.resize(layer_range.size());
  for(int i = 0; i < (int) weight.size(); ++i) {
    Array<double, 1> play(psub.value(layer_range[i]));
    weight[i].value.resize(play.shape());
    weight[i].units = units::Pa;
    weight[i].value = 0;
    for(int j = 0; j < play.rows() - 1; j += 2) {
      double delta = play(j + 2) - play(j);
      weight[i].value(j) += delta / 6.0;
      weight[i].value(j + 1) += 4 * delta / 6.0;
      weight[i].value(j + 2) += delta / 6.0;
    }
  }

  // For each layer, store the nonzero columns of p1 and p2.
  pressure_nonzero_column.clear();
  p1_grad.clear();
  p2_grad.clear();
  for(int i = 0; i < number_layer(); ++i) {
    Array<double, 1> p1g = pgrid.value(i + 1).gradient();
    Array<double, 1> p2g = pgrid.value(i).gradient();
    std::vector<int> nzero_col;
    for(int j = 0; j < p1g.rows(); ++j)
      if(abs(p1g(j)) > 1e-20  ||
	 abs(p2g(j)) > 1e-20)
	nzero_col.push_back(j);
    Array<double, 1> p1g_sparse((int) nzero_col.size());
    Array<double, 1> p2g_sparse(p1g_sparse.shape());
    int k = 0;
    BOOST_FOREACH(int j, nzero_col) {
      p1g_sparse(k) = p1g(j);
      p2g_sparse(k) = p2g(j);
      ++k;
    }
    pressure_nonzero_column.push_back(nzero_col);
    p1_grad.push_back(p1g_sparse);
    p2_grad.push_back(p2g_sparse);
  }
}

//----------------------------------------------------------------
/// This calculates all the tau_gas that is independent of wn. We
/// cache this so we don't recalculate this for every high resolution
/// grid point, instead we only calculate when the Pressure,
/// Temperature, Altitude, or VMR has changed.
//----------------------------------------------------------------

void AbsorberAbsco::fill_tau_gas_cache() const
{
  if(!cache_tau_gas_stale)
    return;
  
  firstIndex i1; secondIndex i2; thirdIndex i3;
  Range ra(Range::all());

  //----------------------------------------------------------------
  // Fill in pressure, temperature, gravity and mixing ratio for
  // the sublayers.
  //----------------------------------------------------------------

  create_sublayer();
  Array<AutoDerivative<double>, 3> 
    integrand_independent_wn_sub_t(number_spectrometer(), number_species(), 
			       psub.rows());
  for(int i = 0; i < integrand_independent_wn_sub_t.rows(); ++i)
    for(int j = 0; j < integrand_independent_wn_sub_t.cols(); ++j)
      for(int k = 0; k < integrand_independent_wn_sub_t.depth(); ++k)
	integrand_independent_wn_sub_t(i, j, k) = 
	  integrand_independent_wn(i, j, psub(k)).
	  convert(Unit("cm^-2 / Pa")).value;
  integrand_independent_wn_sub.units = Unit("cm^-2 / Pa");
  integrand_independent_wn_sub.value.
    reference(ArrayAd<double, 3>(integrand_independent_wn_sub_t));

  // The integrand tends to be fairly sparse. Set up the nonzero
  // columns and compress jacobian for each layer.
  integrand_nonzero_column.clear();
  integrand_jac.clear();
  for(int i = 0; i < number_layer(); ++i) {
    Range r = layer_range[i];
    ArrayAd<double, 3> piv(integrand_independent_wn_sub.value(ra, ra,r));
    Array<double, 4> jac = piv.jacobian();
    std::vector<int> nz;
    for(int j = 0; j < piv.number_variable(); ++j) {
      double mv = max(abs(jac(ra, ra, ra, j)));
      if(mv > 1e-20)
	nz.push_back(j);
    }
    int k = 0;
    Array<double, 4> ijac(jac.rows(), jac.cols(), jac.depth(), (int) nz.size());
    BOOST_FOREACH(int m, nz) {
      ijac(ra, ra, ra, k) = jac(ra,ra,ra,m);
      ++k;
    }
    integrand_nonzero_column.push_back(nz);
    integrand_jac.push_back(ijac);
  }

  // Set up Absco interpolator for our sublayers
  Array<AutoDerivative<double>, 1> tsub_t(psub.rows());
  Array<AutoDerivative<double>, 1> h2o_vmr_t(psub.rows());
  for(int i = 0; i < psub.rows(); ++i) {
    tsub_t(i) = temp_func(psub(i)).convert(units::K).value;
    h2o_vmr_t(i) = h2o_vmr_func(psub(i)).value;
  }
  ArrayAdWithUnit<double, 1> tsub(ArrayAd<double,1>(tsub_t), units::K);
  ArrayAdWithUnit<double, 1> h2o_vmr(ArrayAd<double,1>(h2o_vmr_t), 
				     units::dimensionless);
  // Absco derivatives tend to be fairly sparse. They depend only on
  // tsub and h2o_vmr, so we can go through and see which columns are
  // nonzero.
  absco_nonzero_column.clear();
  for(int i = 0; i < std::max(tsub.value.number_variable(),
			      h2o_vmr.value.number_variable());
      ++i) {
    double mv1;
    double mv2;
    if(tsub.value.is_constant())
      mv1 = 0;
    else
      mv1 = max(abs(tsub.value.jacobian()(ra, i)));
    if(h2o_vmr.value.is_constant())
      mv2 = 0;
    else
      mv2 = max(abs(h2o_vmr.value.jacobian()(ra, i)));
    if(mv1 > 1e-20 || mv2 > 1e-20)
      absco_nonzero_column.push_back(i);
  }
  Array<double, 2> tsub_sjac(tsub.value.rows(), 
			     (int) absco_nonzero_column.size());
  Array<double, 2> h2o_vmr_sjac(tsub_sjac.shape());
  int k = 0;
  tsub_sjac = 0;
  h2o_vmr_sjac = 0;
  BOOST_FOREACH(int m, absco_nonzero_column) {
    if(!tsub.value.is_constant())
      tsub_sjac(ra, k) = tsub.value.jacobian()(ra, m);
    if(!h2o_vmr.value.is_constant())
      h2o_vmr_sjac(ra, k) = h2o_vmr.value.jacobian()(ra, m);
    ++k;
  }
  
  ArrayAdWithUnit<double, 1> 
    tsub_sparse(ArrayAd<double,1>(tsub.value.value(), tsub_sjac), units::K);
  ArrayAdWithUnit<double, 1> 
    h2o_vmr_sparse(ArrayAd<double,1>(h2o_vmr.value.value(), h2o_vmr_sjac),
		   units::dimensionless);

  absco_interp.resize(gas_absorption.size());
  for(int i = 0; i < (int) absco_interp.size(); ++i) {
    boost::shared_ptr<Absco> a = 
      boost::dynamic_pointer_cast<Absco>(gas_absorption[i]);
    if(!a)
      throw Exception("AbsorberAbsco requires Absco, not just any GasAbsorption object");
    absco_interp[i].reset(new AbscoInterpolator(a, psub, tsub_sparse, 
						h2o_vmr_sparse));
  }

  taug.resize((int) layer_range.size(), number_species(), 
	      std::max(integrand_independent_wn_sub.number_variable(),
		       tsub.number_variable()));
  taug.jacobian() = 0;
  cache_tau_gas_stale = false;
}

//----------------------------------------------------------------
/// Calculate tau_gas. Much of the work is cached, because it is
/// independent of wn and spec_index, and only changes when the
/// Pressure, Temperature, or Altitude changes. This version skips the
/// calculation of derivatives
//----------------------------------------------------------------

Array<double, 2> 
AbsorberAbsco::tau_gas_nder(double wn, int spec_index) const
{
  firstIndex i1; secondIndex i2; thirdIndex i3;
  Range ra(Range::all());
  fill_tau_gas_cache();		// Calculate part independent of wn if
				// needed.

  //----------------------------------------------------------------
  // Calculate tau for each sublayer, and then combine to generate
  // taug. 
  //
  // Note that this calculation is a bottle neck, since it gets called
  // for every high resolution wavelength. To speed the calculation up,
  // we have explicit calculation of the Jacobian, rather than
  // automatically doing this with AutoDerivative. This code is less
  // clear, but faster.
  // 
  // You can time this by using the unit test
  // atmosphere_oco/optical_depth_timing (or of course by just
  // profiling a l2_fp run).
  //----------------------------------------------------------------

  for(int i = 0; i < taug.cols(); ++i)
    if(gas_absorption[i]->have_data(wn)) {
      Array<double, 1> abcsub = 
	absco_interp[i]->absorption_cross_section_noderiv(wn);
      for(int j = 0; j < taug.rows(); ++j) {
        Range r = layer_range[j];
        Array<double, 1> av(abcsub(r));
        Array<double, 1> piv(integrand_independent_wn_sub.value.value()(spec_index,i,r));
        taug.value()(j, i) = sum(av * piv * weight[j].value);
      }
    } else
      for(int j =0; j < taug.rows(); ++j) {
        taug.value()(j, i) = 0;
      }
  return taug.value();
}

//----------------------------------------------------------------
/// Calculate tau_gas. Much of the work is cached, because it is
/// independent of wn and spec_index, and only changes when the
/// Pressure, Temperature, or Altitude changes. This version includes the
/// calculation of derivatives.
//----------------------------------------------------------------

ArrayAd<double, 2> 
AbsorberAbsco::tau_gas_der(double wn, int spec_index) const
{
  firstIndex i1; secondIndex i2; thirdIndex i3;
  Range ra(Range::all());
  fill_tau_gas_cache();		// Calculate part independent of wn if
				// needed.

  //----------------------------------------------------------------
  // Calculate tau for each sublayer, and then combine to generate
  // taug. 
  //
  // Note that this calculation is a bottle neck, since it gets called
  // for every high resolution wavelength. To speed the calculation up,
  // we have explicit calculation of the Jacobian, rather than
  // automatically doing this with AutoDerivative. This code is less
  // clear, but faster.
  // 
  // You can time this by using the unit test
  // atmosphere_oco/optical_depth_timing (or of course by just
  // profiling a l2_fp run).
  //----------------------------------------------------------------

  // Temporary, if this works we can cache this
  ArrayAdWithUnit<double, 1> pgrid(press->pressure_grid().convert(Unit("Pa")));

  for(int i = 0; i < taug.cols(); ++i)
    if(gas_absorption[i]->have_data(wn)) {
      ArrayAd<double, 1> abcsub = 
	absco_interp[i]->absorption_cross_section_deriv(wn);
      for(int j = 0; j < taug.rows(); ++j) {
        Range r = layer_range[j];
        ArrayAd<double, 1> av(abcsub(r));
        Array<double, 1> piv(integrand_independent_wn_sub.value(spec_index,i,r).value());
	Array<double, 1> jacv(taug.jacobian()(j, i, Range::all()));
	taug.value()(j, i) = sum(av.value() * piv * weight[j].value);

	// Zero out taug so we can add in the jacobian. Note that the 
        // rest of the Jacobian is already zero from when we sized it.
	BOOST_FOREACH(int m, pressure_nonzero_column[j])
	  jacv(m) = 0;
	
	BOOST_FOREACH(int m, absco_nonzero_column)
	  jacv(m) = 0;

	// Can skip integrand_nonzero_column because we do that first,
	// and just assign to it.

	// BOOST_FOREACH(int m, integrand_nonzero_column[j])
	//   jacv(m) = 0;

	int k = 0;
	BOOST_FOREACH(int m, integrand_nonzero_column[j]) {
	  Array<double, 1> intjac(integrand_jac[j](spec_index, i, ra, k));
	  jacv(m) = sum(av.value() * intjac * weight[j].value);
	  ++k;
	}

	k = 0;
	BOOST_FOREACH(int m, absco_nonzero_column) {
	  Array<double, 1> avjac(av.jacobian()(ra, k));
	  jacv(m) += sum(avjac * piv * weight[j].value);
	  ++k;
	}

	// If p1 and p2 are the pressures at the edge of the layer,
	// then we want to calculate (d tau / dp1) * (dp1 / dstate) +  
        // d (tau / dp2) * dp2 / dstate.  The fundamental theorem of
	// calculus then gives us integrand(p1) * dp 1 /state -
	// integrand(p2) * dp 1 / state

	// Note assumption here that ends of range are for p1 and
	// p2. This is true for Simpsons rule, but if you change
	// create_sublayer this may need to be modified.
	int i1 = av.rows() - 1;
	int i2 = 0;
	k = 0;
	BOOST_FOREACH(int m, pressure_nonzero_column[j]) {
	  jacv(m) += av.value()(i1) * piv(i1) * p1_grad[j](k) -
	    av.value()(i2) * piv(i2) * p2_grad[j](k);
	  ++k;
	}
      }
    } else
      for(int j =0; j < taug.rows(); ++j) {
        taug.value()(j, i) = 0;
        taug.jacobian()(j, i, Range::all()) = 0;
      }
  return taug;
}

//-----------------------------------------------------------------------
/// Create an absorber. 
//-----------------------------------------------------------------------

AbsorberAbsco::AbsorberAbsco
(const std::vector<boost::shared_ptr<AbsorberVmr> > Vmr,
 const boost::shared_ptr<Pressure>& Press,
 const boost::shared_ptr<Temperature>& Temp,
 const std::vector<boost::shared_ptr<Altitude> >& Alt,
 const std::vector<boost::shared_ptr<GasAbsorption> >& Gas_absorption,
 const boost::shared_ptr<Constant>& C,
 int Nsub)
: press(Press), temp(Temp), alt(Alt), vmr(Vmr), 
  gas_absorption(Gas_absorption), c(C), nsub(Nsub), 
  cache_tau_gas_stale(true)
{
  Range ra(Range::all());
  press->add_observer(*this);
  temp->add_observer(*this);
  BOOST_FOREACH(boost::shared_ptr<Altitude>& a, alt)
    a->add_observer(*this);
  BOOST_FOREACH(boost::shared_ptr<AbsorberVmr>& a, vmr)
    a->add_observer(*this);
  h2o_index = gas_index("H2O");
  // Right now, we only support H2O as a broadener.
  BOOST_FOREACH(boost::shared_ptr<GasAbsorption>& ga, gas_absorption)
    if(ga->broadener_name() != "" &&
       ga->broadener_name() != "h2o") {
      Exception e;
      e << "Right now, we only support H2O as a broadener. "
	<< "The GasAbsorption has a broadener of \"" 
	<< ga->broadener_name() << "\"";
      throw e;
    }
}

void AbsorberAbsco::notify_add(StateVector& Sv)
{
  BOOST_FOREACH(boost::shared_ptr<AbsorberVmr>& a, vmr)
    Sv.add_observer(*a);
}

void AbsorberAbsco::notify_remove(StateVector& Sv)
{
  BOOST_FOREACH(boost::shared_ptr<AbsorberVmr>& a, vmr)
    Sv.remove_observer(*a);
}

// See base class for description
ArrayAd<double, 2> 
AbsorberAbsco::optical_depth_each_layer(double wn, int spec_index) const
{
  FunctionTimer ft(timer.function_timer());
  return tau_gas_der(wn,spec_index);
}

Array<double, 2> 
AbsorberAbsco::optical_depth_each_layer_nder(double wn, int spec_index) const
{
  FunctionTimer ft(timer.function_timer());
  return tau_gas_nder(wn,spec_index);
}

// See base class for description
boost::shared_ptr<Absorber> AbsorberAbsco::clone() const
{
  boost::shared_ptr<Pressure> pressure_clone = press->clone();
  boost::shared_ptr<Temperature> temperature_clone = 
    temp->clone(pressure_clone);
  std::vector<boost::shared_ptr<Altitude> > alt_clone;
  BOOST_FOREACH(const boost::shared_ptr<Altitude>& a, alt)
    alt_clone.push_back(a->clone(pressure_clone, temperature_clone));
  return clone(pressure_clone, temperature_clone, alt_clone);
}

// See base class for description
boost::shared_ptr<Absorber> AbsorberAbsco::clone
(const boost::shared_ptr<Pressure>& Press,
 const boost::shared_ptr<Temperature>& Temp,
 const std::vector<boost::shared_ptr<Altitude> >& Alt) const
{
  std::vector<boost::shared_ptr<AbsorberVmr> > vmr_clone;
  BOOST_FOREACH(const boost::shared_ptr<AbsorberVmr>& a, vmr)
    vmr_clone.push_back(a->clone(Press));
  boost::shared_ptr<Absorber> res
    (new AbsorberAbsco(vmr_clone, Press, Temp, Alt, 
		       gas_absorption, c, nsub));
  return res;
}

//-----------------------------------------------------------------------
/// Returns specific humidity by layer
//-----------------------------------------------------------------------

ArrayAdWithUnit<double, 1> AbsorberAbsco::specific_humidity_layer() const
{
  ArrayAdWithUnit<double, 1> pgrid(press->pressure_grid());
  ArrayAdWithUnit<double, 1> spec_hum;
  spec_hum.units = units::dimensionless;
  Array<AutoDerivative<double>, 1> spec_hum_t(pgrid.rows() - 1);
  DoubleWithUnit epsilon = c->molar_weight_water() / c->molar_weight_dry_air();
  for(int i = 0; i < spec_hum_t.rows(); ++i) {
    AutoDerivativeWithUnit<double> dp = pgrid(i + 1) - pgrid(i);
    AutoDerivativeWithUnit<double> play = pgrid(i) + dp / 2;
    AutoDerivativeWithUnit<double> h2o_lay
      ((h2o_index < 0 ? 0 : vmr[h2o_index]->volume_mixing_ratio(play.value)),
       units::dimensionless);
    spec_hum_t(i) = (h2o_lay / (1.0 / epsilon + h2o_lay)).
      convert(spec_hum.units).value;
  }
  spec_hum.value.reference(ArrayAd<double, 1>(spec_hum_t));
  return spec_hum;
}

//-----------------------------------------------------------------------
/// This is a helper function to compute the common part of the
/// dry air mass and wet air mass routines.
/// It returns dry air molecular density (per square meter) in the case
/// where there is no water vapor.
/// Returns units of molecules m^-2
//-----------------------------------------------------------------------

ArrayAdWithUnit<double, 1> AbsorberAbsco::dry_air_molecular_density_layer() const
{
  ArrayAdWithUnit<double, 1> pgrid(press->pressure_grid());
  // We pick the first spectrometer to calculate gravity out. For
  // GOSAT this doesn't matter since all the spectrometers look at the
  // same point. For a different instrument, this might matter and
  // we'll need to revisit this.
  int spec_index = 0;
  ArrayAdWithUnit<double, 1> mol_dens;
  mol_dens.units = Unit("m^-2");
  Array<AutoDerivative<double>, 1> mol_dens_t(pgrid.rows() - 1);
  for(int i = 0; i < mol_dens_t.rows(); ++i) {
    AutoDerivativeWithUnit<double> dp = pgrid(i + 1) - pgrid(i);
    AutoDerivativeWithUnit<double> play(pgrid(i) + dp / 2);
    AutoDerivativeWithUnit<double> grav_lay = alt[spec_index]->gravity(play);
    mol_dens_t(i) = (dp / (c->molar_weight_dry_air() * grav_lay) * 
		     c->avogadro_constant()).convert(mol_dens.units).value;
  }
  mol_dens.value.reference(ArrayAd<double, 1>(mol_dens_t));
  return mol_dens;
}

//-----------------------------------------------------------------------
/// This is the dry air column thickness by layer. This is the size
/// of pressure_grid() - 1.
//-----------------------------------------------------------------------

ArrayAdWithUnit<double, 1> AbsorberAbsco::dry_air_column_thickness_layer() const
{
  ArrayAdWithUnit<double, 1> spec_hum(specific_humidity_layer());
  ArrayAdWithUnit<double, 1> mol_dens(dry_air_molecular_density_layer());
  ArrayAdWithUnit<double, 1> dry_am;
  dry_am.units = mol_dens.units;
  Array<AutoDerivative<double>, 1> dry_am_t(mol_dens.rows());
  for(int i = 0; i < dry_am_t.rows(); ++i)
    dry_am_t(i) = (mol_dens(i) * (1 - spec_hum(i))).convert(dry_am.units).value;
  dry_am.value.reference(ArrayAd<double, 1>(dry_am_t));
  return dry_am;
}

//-----------------------------------------------------------------------
/// This is the wet air column thickness by layer. This is the size
/// of pressure_grid() - 1.
//-----------------------------------------------------------------------

ArrayAdWithUnit<double, 1> 
AbsorberAbsco::wet_air_column_thickness_layer() const
{
  ArrayAdWithUnit<double, 1> spec_hum(specific_humidity_layer());
  ArrayAdWithUnit<double, 1> mol_dens(dry_air_molecular_density_layer());
  DoubleWithUnit epsilon = c->molar_weight_water() / c->molar_weight_dry_air();
  ArrayAdWithUnit<double, 1> wet_am;
  wet_am.units = mol_dens.units;
  Array<AutoDerivative<double>, 1> wet_am_t(mol_dens.rows());
  for(int i = 0; i < wet_am_t.rows(); ++i)
    wet_am_t(i) = 
      (mol_dens(i) * (1 + spec_hum(i) * (1 - epsilon) / epsilon)).
      convert(wet_am.units).value;
  wet_am.value.reference(ArrayAd<double, 1>(wet_am_t));
  return wet_am;
}

//-----------------------------------------------------------------------
/// This is the pressure weighting function by layer. This is the size
/// of pressure_grid() - 1.
//-----------------------------------------------------------------------

ArrayAd<double, 1> AbsorberAbsco::pressure_weighting_function_layer() const
{
  Array<AutoDerivative<double>, 1> dry_am(dry_air_column_thickness_layer().
					  value.to_array());
  Array<AutoDerivative<double>, 1> pwf(dry_am.rows());
  pwf = dry_am / sum(dry_am);
  return ArrayAd<double, 1>(pwf);
}

//-----------------------------------------------------------------------
/// This is the pressure weighting function by grid level. This is
/// calculated so that:
/// XCO2 = (co2 on grid levels) dot (press_wf_lev)
//-----------------------------------------------------------------------

ArrayAd<double, 1> AbsorberAbsco::pressure_weighting_function_grid() const
{
  // Note assumption here that pressure varies linearly over layer
  ArrayAd<double, 1> pwlay = pressure_weighting_function_layer();
  int nactive = press->number_level();
  Array<AutoDerivative<double>, 1> pwlev(nactive);
  pwlev(0) = pwlay(0) / 2;
  for(int i = 1; i < nactive - 1; ++i)
    pwlev(i) = (pwlay(i-1) + pwlay(i)) / 2;
  pwlev(nactive - 1) = pwlay(nactive - 2) / 2;
  return ArrayAd<double, 1>(pwlev);
}

//-----------------------------------------------------------------------
/// This is the column thickness of a gas by layer. This is the size
/// of pressure_grid() - 1.
//-----------------------------------------------------------------------

ArrayAdWithUnit<double, 1> 
AbsorberAbsco::gas_column_thickness_layer(const std::string& Gas_name) const
{
  ArrayAdWithUnit<double, 1> pgrid(press->pressure_grid());
  ArrayAdWithUnit<double, 1> dry_am = dry_air_column_thickness_layer();
  ArrayAdWithUnit<double, 1> gas_thickness;
  gas_thickness.units = dry_am.units;
  Array<AutoDerivative<double>, 1> gas_thickness_t(dry_am.rows());
  for(int i = 0; i < gas_thickness_t.rows(); ++i) {
    AutoDerivativeWithUnit<double> dp = pgrid(i + 1) - pgrid(i);
    AutoDerivativeWithUnit<double> play = pgrid(i) + dp / 2;
    gas_thickness_t(i) = 
      (dry_am(i) * absorber_vmr(Gas_name)->volume_mixing_ratio(play.convert(units::Pa).value)).
      convert(gas_thickness.units).value;
  }
  gas_thickness.value.reference(ArrayAd<double, 1>(gas_thickness_t));
  return gas_thickness;
}

//-----------------------------------------------------------------------
/// This is the total column thickness of a gas.
//-----------------------------------------------------------------------

AutoDerivativeWithUnit<double> 
AbsorberAbsco::gas_total_column_thickness(const std::string& Gas_name) const
{
  ArrayAdWithUnit<double, 1> gas_thickness = 
    gas_column_thickness_layer(Gas_name);
  return AutoDerivativeWithUnit<double>(sum(gas_thickness.value.to_array()), 
					gas_thickness.units);
}

// See base class for description
AutoDerivative<double> AbsorberAbsco::xgas(const std::string& Gas_name) const
{ 
  ArrayAd<double, 1> pwf(pressure_weighting_function_grid());
  ArrayAdWithUnit<double, 1> pgrid = press->pressure_grid();
  AutoDerivative<double> res;
  res = pwf(0) * absorber_vmr(Gas_name)->volume_mixing_ratio(pgrid(0).convert(units::Pa).value);
  for(int i = 1; i < pwf.rows(); ++i)
    res += pwf(i) * absorber_vmr(Gas_name)->volume_mixing_ratio(pgrid(i).convert(units::Pa).value);
  return res;
}

//-----------------------------------------------------------------------
/// Returns the simple average volume mixing ratio for a gas.
/// The units returned are mole / mole, same as for vmr.
//-----------------------------------------------------------------------

AutoDerivative<double> AbsorberAbsco::average_vmr(const std::string& Gas_name) const
{
  ArrayAdWithUnit<double, 1> pgrid(press->pressure_grid());
  AutoDerivative<double> avg_vmr = absorber_vmr(Gas_name)->volume_mixing_ratio(pgrid(0).convert(units::Pa).value);
  for(int i = 1; i < press->number_level(); ++i) {
    avg_vmr += absorber_vmr(Gas_name)->volume_mixing_ratio(pgrid(i).convert(units::Pa).value);
  }
  return avg_vmr / press->number_level();
}

void AbsorberAbsco::print(std::ostream& Os) const 
{ 
  OstreamPad opad(Os, "    ");
  Os << "AbsorberAbsco:\n";
  for(int i = 0; i < (int) gas_absorption.size(); ++i) {
    Os << "  Gas Absorber[" << i << "]:\n";
    opad << *gas_absorption[i] << "\n";
    opad.strict_sync();
    Os << "  Absorber VMR[" << i << "]:\n";
    opad << *vmr[i] << "\n";
    opad.strict_sync();
  }
}

//----------------------------------------------------------------
/// Return the gravity we use for each sublayer. This is meant for
/// diagnostic purposes.
//----------------------------------------------------------------

ArrayAdWithUnit<double, 1> AbsorberAbsco::gravity_sublayer(int Spec_index) const
{
  fill_tau_gas_cache(); 
  range_check(Spec_index, 0, number_spectrometer());
  blitz::Array<AutoDerivative<double>, 1> gsub(psub.rows());
  for(int i = 0; i < gsub.rows(); ++i)
    gsub(i) = gravity_func(Spec_index, psub(i)).convert(Unit("m/s^2")).value;
  return ArrayAdWithUnit<double, 1>(ArrayAd<double, 1>(gsub),
				    Unit("m/s^2"));
}

//----------------------------------------------------------------
/// Return the volume mixing ratio we use for each sublayer. This
/// is meant for diagnostic purposes.
//----------------------------------------------------------------

ArrayAdWithUnit<double, 1> AbsorberAbsco::vmr_sublayer(const std::string& Gas_name) const
{     
  fill_tau_gas_cache(); 
  int j = gas_index(Gas_name);
  if(j < 0) {
    Exception err;
    err << "Gas named " << Gas_name << " not present in vmr list";
    throw err;
  }
  Array<AutoDerivative<double>, 1> vmrsub(psub.rows());
  for(int i = 0; i < psub.rows(); ++i)
    vmrsub(i) = vmr_func(j, psub(i)).value;
  return ArrayAdWithUnit<double, 1>(ArrayAd<double, 1>(vmrsub),
				    units::dimensionless);
}

// See base class for description
boost::shared_ptr<AbsorberVmr> AbsorberAbsco::absorber_vmr(const std::string& Gas_name) const
{ 
  int i = gas_index(Gas_name);
  if(i < 0) {
    Exception err;
    err << "Gas named " << Gas_name << " not present in vmr list";
    throw err;
  }
  return vmr[i];
}

//----------------------------------------------------------------
/// Return GasAbsorption as a pointer.
//----------------------------------------------------------------

boost::shared_ptr<GasAbsorption> AbsorberAbsco::gas_absorption_ptr
(const std::string& Gas_name) const
{
  int i = gas_index(Gas_name);
  if(i < 0) {
    Exception err;
    err << "Gas named " << Gas_name << " not present in vmr list";
    throw err;
  }
  return gas_absorption[i];
}

//----------------------------------------------------------------
/// Return the temperature we use for each sublayer. This is meant for
/// diagnostic purposes.
//----------------------------------------------------------------

ArrayAdWithUnit<double, 1> AbsorberAbsco::temperature_sublayer() const
{ 
  fill_tau_gas_cache(); 
  Array<AutoDerivative<double>, 1> tsub_t(psub.rows());
  for(int i = 0; i < psub.rows(); ++i)
    tsub_t(i) = temp_func(psub(i)).convert(units::K).value;
  return ArrayAdWithUnit<double, 1>(ArrayAd<double, 1>(tsub_t),
				    units::K);
}

//----------------------------------------------------------------
/// Return the H2O volume mixing ratio we use for each sublayer. This
/// is meant for diagnostic purposes.
//----------------------------------------------------------------

ArrayAdWithUnit<double, 1> AbsorberAbsco::h2o_vmr_sublayer() const
{ 
  fill_tau_gas_cache(); 
  Array<AutoDerivative<double>, 1> h2o_vmr_t(psub.rows());
  for(int i = 0; i < psub.rows(); ++i)
    h2o_vmr_t(i) = h2o_vmr_func(psub(i)).value;
  return ArrayAdWithUnit<double, 1>(ArrayAd<double, 1>(h2o_vmr_t),
				    units::dimensionless);
}




 




