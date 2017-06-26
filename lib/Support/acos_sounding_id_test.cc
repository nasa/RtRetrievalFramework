#include "unit_test_support.h"
#include "acos_sounding_id.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(acos_sounding_id, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  HdfFile h(test_data_dir() + "l1b.h5");
  AcosSoundingId sid_p(h, "20090725020316", AcosSoundingId::P_SOUNDING);
  BOOST_CHECK_EQUAL(sid_p.frame_number(), 151);
  BOOST_CHECK_EQUAL(sid_p.sounding_id(), (int64_t) 20090725020316);
  BOOST_CHECK_EQUAL(sid_p.sounding_number(), 0);
  AcosSoundingId sid_s(h, "20090725020316", AcosSoundingId::S_SOUNDING);
  BOOST_CHECK_EQUAL(sid_s.frame_number(), 151);
  BOOST_CHECK_EQUAL(sid_s.sounding_id(), (int64_t) 20090725020316);
  BOOST_CHECK_EQUAL(sid_s.sounding_number(), 1);
}

BOOST_AUTO_TEST_CASE(create)
{ 
  HdfFile h(test_data_dir() + "l1b.h5");
  std::vector<boost::shared_ptr<HdfSoundingId> > res;
  res = AcosSoundingId::create(h, "20090725020316");
  BOOST_CHECK_EQUAL((int) res.size(), 2);
  BOOST_CHECK_EQUAL(res[0]->frame_number(), 151);
  BOOST_CHECK_EQUAL(res[0]->sounding_id(), (int64_t) 20090725020316);
  BOOST_CHECK_EQUAL(res[0]->sounding_number(), 0);
  BOOST_CHECK_EQUAL(res[1]->frame_number(), 151);
  BOOST_CHECK_EQUAL(res[1]->sounding_id(), (int64_t) 20090725020316);
  BOOST_CHECK_EQUAL(res[1]->sounding_number(), 1);
  res = AcosSoundingId::create(h, "20090725020316P");
  BOOST_CHECK_EQUAL((int) res.size(), 1);
  BOOST_CHECK_EQUAL(res[0]->frame_number(), 151);
  BOOST_CHECK_EQUAL(res[0]->sounding_id(), (int64_t) 20090725020316);
  BOOST_CHECK_EQUAL(res[0]->sounding_number(), 0);
  res = AcosSoundingId::create(h, "20090725020316S");
  BOOST_CHECK_EQUAL((int) res.size(), 1);
  BOOST_CHECK_EQUAL(res[0]->frame_number(), 151);
  BOOST_CHECK_EQUAL(res[0]->sounding_id(), (int64_t) 20090725020316);
  BOOST_CHECK_EQUAL(res[0]->sounding_number(), 1);
}

BOOST_AUTO_TEST_SUITE_END()
