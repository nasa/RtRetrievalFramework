#include "unit_test_support.h"
#include "output.h"

using namespace FullPhysics;

class OutputTest : public Output {
public:
  std::map<std::string, int> val;
protected:
  void start_write()
  { }
  void end_write() { }
  void write_data(const std::string& Dataset_name, int Val)
  { val[Dataset_name] = Val;}
  void write_data(const std::string& Dataset_name, const std::string& Val)
  { }
  void write_data(const std::string& Dataset_name, const char* Val)
  { }
  void write_data(const std::string& Dataset_name, int64_t Val)
  { }
  void write_data(const std::string& Dataset_name, double Val)
  { }
  void write_data(const std::string& Dataset_name, 
		  const blitz::Array<int, 1>& Val)
  { }
  void write_data(const std::string& Dataset_name, 
		  const blitz::Array<std::string, 1>& Val)
  { }
  void write_data(const std::string& Dataset_name, 
		  const blitz::Array<const char*, 1>& Val)
  { }
  void write_data(const std::string& Dataset_name, 
		  const blitz::Array<double, 1>& Val)
  { }
  void write_data(const std::string& Dataset_name, 
		  const blitz::Array<int, 2>& Val)
  { }
  void write_data(const std::string& Dataset_name, 
		  const blitz::Array<std::string, 2>& Val)
  { }
  void write_data(const std::string& Dataset_name, 
		  const blitz::Array<const char*, 2>& Val)
  { }
  void write_data(const std::string& Dataset_name, 
		  const blitz::Array<double, 2>& Val)
  { }
  void write_data(const std::string& Dataset_name, 
		  const blitz::Array<int, 3>& Val)
  { }
  void write_data(const std::string& Dataset_name, 
		  const blitz::Array<std::string, 3>& Val)
  { }
  void write_data(const std::string& Dataset_name, 
		  const blitz::Array<const char*, 3>& Val)
  { }
  void write_data(const std::string& Dataset_name, 
		  const blitz::Array<double, 3>& Val)
  { }
};

class CalcData {
public:
  CalcData(int D) : d(D) {}
  int d_val() const {return d;}
private:
  int d;
};

BOOST_FIXTURE_TEST_SUITE(output, GlobalFixture)

BOOST_AUTO_TEST_CASE(read)
{
  OutputTest ot;
  add_file_to_cleanup("out.h5");
  boost::shared_ptr<CalcData> c1(new CalcData(1));
  boost::shared_ptr<CalcData> c2(new CalcData(2));
  ot.register_data_source("C1", &CalcData::d_val, c1);
  ot.register_data_source("C2", &CalcData::d_val, c2);
  ot.write();
  BOOST_CHECK_EQUAL(ot.val["C1"], 1);
  BOOST_CHECK_EQUAL(ot.val["C2"], 2);
}

BOOST_AUTO_TEST_SUITE_END()
