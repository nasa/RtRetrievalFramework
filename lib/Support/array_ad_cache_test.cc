#include "unit_test_support.h"
#include "array_ad_cache.h"
#include <cstdlib>

using namespace FullPhysics;
using namespace blitz;

class ArrayAdCacheFixture: public GlobalFixture {
  public:
    ArrayAdCacheFixture() 
    {
      std::srand(8675309);
      int num_check = 5;
      for(int x = 0; x < num_check; x++) {
        double rand_val = ((double)rand()/(double)RAND_MAX);

        boost::shared_ptr<ArrayAd<double,1> > ad_ptr(new ArrayAd<double, 1>());
        ad_ptr->resize(1, 0);
        ad_ptr->value()(0) = rand_val;
 
        test_keys.push_back((double) x / (double) num_check);
        test_values.push_back(ad_ptr);
      }
      
    }

    std::vector<double> test_keys;
    std::vector<boost::shared_ptr<ArrayAd<double,1> > > test_values;
  private:
};

BOOST_FIXTURE_TEST_SUITE(array_ad_cache, ArrayAdCacheFixture)

BOOST_AUTO_TEST_CASE(map_cache)
{
  ArrayAdMapCache<double, double, 1> map_cache;
  // Insert value
  for(int tst_idx = 0; tst_idx < (int) test_values.size(); tst_idx++)
  {
    map_cache.insert(test_keys[tst_idx], *test_values[tst_idx]);
  }
  // Test retrieving them
  for(int tst_idx = 0; tst_idx < (int) test_values.size(); tst_idx++)
  {
    BOOST_CHECK_EQUAL(map_cache.is_valid(test_keys[tst_idx]), true); 
    BOOST_CHECK_CLOSE(test_values[tst_idx]->value()(0), map_cache[test_keys[tst_idx]].value()(0), 1e-15);
  }
  // Check we error when cache is cleared
  map_cache.clear();
  try {
    map_cache[test_keys[0]];
    BOOST_ERROR("Should have thrown exception");
  } catch (const Exception& e) {
    BOOST_CHECK(true);
  }
  // Check that we don't just return some random value
  map_cache.insert(test_keys[0], *test_values[0]);
  BOOST_CHECK_EQUAL(map_cache.is_valid(test_keys[1]), false); 
  try {
    map_cache[test_keys[1]];
    BOOST_ERROR("Should have thrown exception");
  } catch (const Exception& e) {
    BOOST_CHECK(true);
  }
  // Check that erase works
  map_cache.insert(test_keys[1], *test_values[1]);
  map_cache.erase(test_keys[0]);
  BOOST_CHECK_EQUAL(map_cache.is_valid(test_keys[0]), false); 
  BOOST_CHECK_EQUAL(map_cache.is_valid(test_keys[1]), true); 

}

BOOST_AUTO_TEST_CASE(vector_cache)
{
  ArrayAdVectorCache<double, double, 1> vec_cache;
  // Insert value
  for(int tst_idx = 0; tst_idx < (int) test_values.size(); tst_idx++)
  {
    vec_cache.insert(test_keys[tst_idx], *test_values[tst_idx]);
  }
  // Test retrieving them
  for(int tst_idx = 0; tst_idx < (int) test_values.size(); tst_idx++)
  {
    BOOST_CHECK_EQUAL(vec_cache.is_valid(test_keys[tst_idx]), true); 
    BOOST_CHECK_CLOSE(test_values[tst_idx]->value()(0), vec_cache[test_keys[tst_idx]].value()(0), 1e-15);
  }
  // Check that the cached data is the same size and contents as what we expect
  std::vector<std::pair<double, ArrayAd<double, 1> > > cached_data = vec_cache.cached_data();
  BOOST_CHECK_EQUAL(test_values.size(), cached_data.size());
  for(int tst_idx = 0; tst_idx < (int) test_values.size(); tst_idx++)
  {
    BOOST_CHECK_CLOSE(test_keys[tst_idx], cached_data[tst_idx].first, 1e-15);
    BOOST_CHECK_CLOSE(test_values[tst_idx]->value()(0), cached_data[tst_idx].second.value()(0), 1e-15);
  }
  // Check we error when cache is cleared
  vec_cache.clear();
  BOOST_CHECK_EQUAL((int) vec_cache.cached_data().size(), 0);
  BOOST_CHECK_EQUAL(vec_cache.is_valid(test_keys[0]), false); 
  try {
    vec_cache[test_keys[0]];
    BOOST_ERROR("Should have thrown exception");
  } catch (const Exception& e) {
    BOOST_CHECK(true);
  }
  // Check that we don't just return some random value
  vec_cache.insert(test_keys[0], *test_values[0]);
  BOOST_CHECK_EQUAL(vec_cache.is_valid(test_keys[1]), false); 
  try {
    vec_cache[test_keys[1]];
    BOOST_ERROR("Should have thrown exception");
  } catch (const Exception& e) {
    BOOST_CHECK(true);
  }

  // Check that erase works
  vec_cache.insert(test_keys[1], *test_values[1]);
  vec_cache.erase(test_keys[0]);
  BOOST_CHECK_EQUAL(vec_cache.is_valid(test_keys[0]), false); 
  BOOST_CHECK_EQUAL(vec_cache.is_valid(test_keys[1]), true); 

  // Check that total clear works
  vec_cache.clear(true);
  BOOST_CHECK_EQUAL((int) vec_cache.cached_data().size(), 0);
  BOOST_CHECK_EQUAL(vec_cache.is_valid(test_keys[0]), false); 
 
}

BOOST_AUTO_TEST_SUITE_END()
