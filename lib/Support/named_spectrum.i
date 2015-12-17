// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include <std_vector.i>
%include "common.i"

%{
#include "named_spectrum.h"
#include "observer.h"
%}

%base_import(spectrum)
%import "observer.i"

%fp_shared_ptr(FullPhysics::NamedSpectrum)

%fp_shared_ptr(FullPhysics::Observable<FullPhysics::NamedSpectrum>)
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::NamedSpectrum>)
%fp_shared_ptr(FullPhysics::Observable<std::vector<boost::shared_ptr<FullPhysics::NamedSpectrum> > >)
%fp_shared_ptr(FullPhysics::Observable<boost::shared_ptr<FullPhysics::NamedSpectrum> >)
%template(ObservableNamedSpectrum) FullPhysics::Observable<FullPhysics::NamedSpectrum>;
%template(ObservablePtrNamedSpectrum) FullPhysics::Observable<boost::shared_ptr<FullPhysics::NamedSpectrum> >;
%template(ObserverNamedSpectrum) FullPhysics::Observer<FullPhysics::NamedSpectrum>;

namespace FullPhysics {
class NamedSpectrum: public Spectrum {
public:
  NamedSpectrum(const SpectralDomain& Spec_domain, 
                const SpectralRange& Spec_range, const std::string& Name,
                int Index);

  NamedSpectrum(const Spectrum& Spec, const std::string& Name, int Index);
  %python_attribute(name, virtual const std::string&);
  %python_attribute(index, virtual int);
  std::string print_to_string() const;
};
}
%template(vector_named_spectrum) std::vector<boost::shared_ptr<FullPhysics::NamedSpectrum> >;
%fp_shared_ptr(FullPhysics::Observable<std::vector<boost::shared_ptr<FullPhysics::NamedSpectrum> > >)
%fp_shared_ptr(FullPhysics::Observer<std::vector<boost::shared_ptr<FullPhysics::NamedSpectrum> > >)
%template(ObservableStokesUpdate) FullPhysics::Observable<std::vector<boost::shared_ptr<FullPhysics::NamedSpectrum> > >;
// We make this a director, so we can have Python code that is an
// observer of a stokes upate
%feature("director") FullPhysics::Observer<std::vector<boost::shared_ptr<FullPhysics::NamedSpectrum> > >;
%template(ObserverStokesUpdate) FullPhysics::Observer<std::vector<boost::shared_ptr<FullPhysics::NamedSpectrum> > >;
