#ifndef LIDORT_FIXTURE_H
#define LIDORT_FIXTURE_H

#include "configuration_fixture.h"
#include "lidort_rt.h"

namespace FullPhysics {

/****************************************************************//**
 Base fixture for testing LidortRtDriver
*******************************************************************/

class LidortDriverCommonFixture {
public:
  LidortDriverCommonFixture();
  virtual ~LidortDriverCommonFixture() {};

  /// Lidort driver to use in testing
  boost::shared_ptr<LidortRtDriver> lidort_driver;

  /// solar zenith, view zenith, and azimuth angles used in low_rt and
  /// high_rt.
  blitz::Array<double, 1> sza, zen, azm;
  bool pure_nadir;

  /// Wavenumber array used for testing
  blitz::Array<double, 1> wn_arr;
};

/****************************************************************//**
 Fixture for LidortRtDriver lambertian mode
*******************************************************************/

class LidortDriverLambertianFixture : public LidortDriverCommonFixture {
public:
  LidortDriverLambertianFixture();
  virtual ~LidortDriverLambertianFixture() {};
};


/****************************************************************//**
 Fixture for LidortRtDriver coxmunk mode
*******************************************************************/

class LidortDriverCoxmunkFixture : public LidortDriverCommonFixture {
public:
  LidortDriverCoxmunkFixture();
  virtual ~LidortDriverCoxmunkFixture() {};
};

/****************************************************************//**
 Base fixture for testing LidortRt modes
*******************************************************************/

class LidortRtCommonFixture {
public:
  LidortRtCommonFixture();
  virtual ~LidortRtCommonFixture() {};

  /// Lidort driver to use in testing
  boost::shared_ptr<LidortRt> lidort_rt;

  /// solar zenith, view zenith, and azimuth angles used in low_rt and
  /// high_rt.
  blitz::Array<double, 1> sza, zen, azm;
  bool pure_nadir;

  /// Wavenumber array used for testing
  blitz::Array<double, 1> wn_arr;

  boost::shared_ptr<StokesCoefficient> stokes_coefs;

  std::string pp_expected_filename, ps_expected_filename, pp_and_ss_expected_filename;
};

/****************************************************************//**
  This is a test fixture that creates a Atmosphere, Statevector, and
  LidortRt based on the ConfigurationFixture.
*******************************************************************/
class LidortLambertianFixture : public LidortRtCommonFixture, 
				public ConfigurationFixture {
public:
  LidortLambertianFixture(const std::string& Config_file = "config.lua");
  virtual ~LidortLambertianFixture() {};
};

/****************************************************************//**
  This is a test fixture that creates a Atmosphere, Statevector, and
  LidortRt based on the ConfigurationCoxmunkFixture. Same as 
  LidortLambertianFixture but using coxmunk surface type
*******************************************************************/
class LidortCoxmunkFixture : public LidortRtCommonFixture, 
			     public ConfigurationCoxmunkFixture {
public:
  LidortCoxmunkFixture(const std::string& Config_file = "config_coxmunk.lua");
  virtual ~LidortCoxmunkFixture() {};
};

/****************************************************************//**
  This is a test fixture that creates a Atmosphere, Statevector, and
  LidortRt based on the ConfigurationCoxmunkPlusLambertianFixture. Same as 
  LidortLambertianFixture but using coxmunk + lambertian surface type
*******************************************************************/
class LidortCoxmunkPlusLambertianFixture : public LidortRtCommonFixture, 
                                           public ConfigurationCoxmunkPlusLambertianFixture {
public:
  LidortCoxmunkPlusLambertianFixture();
  virtual ~LidortCoxmunkPlusLambertianFixture() {};
};

/*******************************************************************/

class LidortLowHighCommon {
 public:
   LidortLowHighCommon();
   virtual ~LidortLowHighCommon() {};
   
  /// Low stream version of lidort
  boost::shared_ptr<LidortRt> low_rt;

  /// High stream version of lidort
  boost::shared_ptr<LidortRt> high_rt;

  /// solar zenith, view zenith, and azimuth angles used in low_rt and
  /// high_rt.
  blitz::Array<double, 1> sza, zen, azm;
  bool pure_nadir;

  // Used for converting to final radiance value 
  boost::shared_ptr<StokesCoefficient> stokes_coefs;

  // Use these for all fixtures
  int nstreams_low, nstreams_high, nmoms_low, nmoms_high;
};

/****************************************************************//**
  This is a test fixture that creates a Atmosphere, Statevector, and
  LidortLowHighDriver based on the ConfigurationFixture.
*******************************************************************/
 class LidortLowHighLambertianFixture : public LidortLowHighCommon,
                                        public ConfigurationFixture
  {
public:
  LidortLowHighLambertianFixture();
  virtual ~LidortLowHighLambertianFixture() {}

};

/****************************************************************//**
  This is a test fixture that creates a Atmosphere, Statevector, and
  LidortLowHighDriver based on the ConfigurationFixture. This is very
  similar to LidortLowHighFixture, but we have the polarization turned off
  for this.
*******************************************************************/
class LidortLowHighFixtureNoPolarization : public LidortLowHighCommon,
                                           public ConfigurationFixture {
public:
  LidortLowHighFixtureNoPolarization();
  virtual ~LidortLowHighFixtureNoPolarization() {}
};

/****************************************************************//**
  This is a test fixture that creates a Atmosphere, Statevector, and
  LidortLowHighDriver based on the ConfigurationCoxmunkFixture. Same as 
  LidortLowHighLambertianFixture but using coxmunk surface type
*******************************************************************/
class LidortLowHighCoxmunkFixture : public LidortLowHighCommon,
                                    public ConfigurationCoxmunkFixture {
public:
  LidortLowHighCoxmunkFixture();
  virtual ~LidortLowHighCoxmunkFixture() {}
};

/****************************************************************//**
  This is a test fixture that creates a Atmosphere, Statevector, and
  LidortRt based on the ConfigurationCoxmunkPlusLambertianFixture. Same as 
  LidortLambertianFixture but using coxmunk surface type with
  a lambertian component
*******************************************************************/
class LidortLowHighCoxmunkPlusLambertianFixture : public LidortLowHighCommon,
                                                  public ConfigurationCoxmunkPlusLambertianFixture {
public:
  LidortLowHighCoxmunkPlusLambertianFixture();
  virtual ~LidortLowHighCoxmunkPlusLambertianFixture() {}
};

}
#endif

