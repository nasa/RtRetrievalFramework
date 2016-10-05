// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "observer.h"
%}
%base_import(generic_object)

namespace FullPhysics {
template<class T> class Observable;
template<class T> class Observer : public GenericObject {
public:
  Observer();
  virtual ~Observer();
  virtual void notify_update(const T& Observed_object);
  virtual void notify_add(T& Observed_object);
  virtual void notify_remove(T& Observed_object);
};

template<class T> class Observable : public GenericObject {
public:
  virtual ~Observable();
  void add_observer_and_keep_reference(boost::shared_ptr<Observer<T> >& Obs);
  virtual void add_observer(Observer<T>& Obs) = 0;
  virtual void remove_observer(Observer<T>& Obs) = 0;
protected:
  void notify_update_do(const T& Self);
  void add_observer_do(Observer<T>& Obs, T& t);
  void remove_observer_do(Observer<T>& Obs, T& t);
  void clean_dead_ptr();
};  
}
