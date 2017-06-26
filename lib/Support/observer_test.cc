#include "observer.h"
#include "unit_test_support.h"

using namespace FullPhysics;

class ObservableTest : public Observable<ObservableTest> {
public:
  ObservableTest() :v_(0) {}
  virtual ~ObservableTest() {}
  int v() const { return v_;}
  void v(int newv) { v_ = newv; notify_update_do(*this);}
  virtual void add_observer(Observer<ObservableTest>& Obs) 
  { add_observer_do(Obs, *this);}
  virtual void remove_observer(Observer<ObservableTest>& Obs) 
  { remove_observer_do(Obs, *this);}
private:
  int v_;
};

class ObserverTest : public Observer<ObservableTest> {
public:
  ObserverTest() : myv_(-1) {}
  virtual ~ObserverTest() {}
  int myv() { return myv_; }
  virtual void notify_update(const ObservableTest& T)
  {
    myv_ = T.v();
  }
private:
  int myv_;
};

BOOST_FIXTURE_TEST_SUITE(observer, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  ObservableTest t;
  ObserverTest t2;
  t.add_observer(t2);
  BOOST_CHECK_EQUAL(t2.myv(), -1);
  t.v(10);
  BOOST_CHECK_EQUAL(t2.myv(), 10);
}

BOOST_AUTO_TEST_CASE(proper_cleanup)
{
  ObservableTest t;
  boost::shared_ptr<ObserverTest> t2(new ObserverTest);
  t.add_observer(*t2);
  BOOST_CHECK_EQUAL(t2->myv(), -1);
  t.v(10);
  BOOST_CHECK_EQUAL(t2->myv(), 10);
  t2.reset();
  t.v(10);
}

BOOST_AUTO_TEST_SUITE_END()
