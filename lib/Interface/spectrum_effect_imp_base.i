// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "spectrum_effect_imp_base.h"
#include "instrument.h"
#include "spectrum_sampling.h"
#include "forward_model_spectral_grid.h"
%}
%base_import(sub_state_vector_array)
%base_import(spectrum_effect)

%fp_shared_ptr(FullPhysics::SpectrumEffectImpBase);

namespace FullPhysics {

// Allow this classes to be derived from in Python.
%feature("director") SpectrumEffectImpBase;

// Note, at least for SWIG 2.0.4 a class that is derived from in python 
// needs to declare every virtual function that can be called on it, even 
// if all that happens is the base class to a director is called. This is 
// because this class is used to create the SwigDirector class, and this 
// class needs each of the member functions to direct things properly. It 
// is *not* necessary to add these function to the underlying
// C++, only that you declare them here.
//
// If you miss something, then you will get something like a recursion
// error in python when a virtual function is used that isn't explicitly
// listed here.
//
// This seems like a bug in 2.0.4, if SWIG needs all the member functions
// it should know to make them itself. So perhaps a future version of SWIG
// won't have this same constraint. But for now, this is required.

class SpectrumEffectImpBase: public SubStateVectorArray<SpectrumEffect> {
public:
  // From SpectrumEffectImpBase
  virtual ~SpectrumEffectImpBase();
  virtual boost::shared_ptr<SpectrumEffect> clone() const = 0;

  // From SpectrumEffect
  virtual void apply_effect(Spectrum& Spec, const ForwardModelSpectralGrid& Forward_model_grid) const = 0;

  %sub_state_virtual_func(SpectrumEffect);
protected:
  SpectrumEffectImpBase();
  SpectrumEffectImpBase(const blitz::Array<double, 1>& Coeff, 
			const blitz::Array<bool, 1>& Used_flag);
};
}

