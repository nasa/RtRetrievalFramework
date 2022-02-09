#include "register_lua.h"
#include <iostream>
#include <sstream>
#include "unit.h"

using namespace FullPhysics;

//-----------------------------------------------------------------------
/// Here are some registrations that do not really belong anywhere else.
//-----------------------------------------------------------------------

std::string string_vector_tostring(const std::vector<std::string>& Vec)
{
  std::ostringstream os;
  os << "std::vector<std::string>:" << std::endl;
  for(int i = 0; i < (int) Vec.size(); i++)
    os << "["<<i<<"] = " << Vec[i] << std::endl;
  return os.str();
}

std::string double_vector_tostring(const std::vector<double>& Vec)
{
  std::ostringstream os;
  os << "std::vector<std::double>:" << std::endl;
  for(int i = 0; i < (int) Vec.size(); i++)
    os << "["<<i<<"] = " << Vec[i] << std::endl;
  return os.str();
}

std::string int_vector_tostring(const std::vector<int>& Vec)
{
  std::ostringstream os;
  os << "std::vector<std::int>:" << std::endl;
  for(int i = 0; i < (int) Vec.size(); i++)
    os << "["<<i<<"] = " << Vec[i] << std::endl;
  return os.str();
}

std::string string_value(const std::vector<std::string>& Vec, int index) {
  return Vec.at(index);
}

double double_value(const std::vector<double>& Vec, int index) {
  return Vec.at(index);
}

int int_value(const std::vector<int>& Vec, int index) {
  return Vec.at(index);
}

// typedef to distinguish between copying value or moving value (C++11) push_back prototoypes 
typedef void(std::vector<std::string>::*pbt1)(const std::vector<std::string>::value_type&);
typedef void(std::vector<double>::*pbt2)(const std::vector<double>::value_type&);
typedef void(std::vector<int>::*pbt3)(const std::vector<int>::value_type&);

REGISTER_LUA_CLASS_NAME(std::vector<std::string>, VectorString)
.def(luabind::constructor<>())
.def("size", &std::vector<std::string>::size)
.def("push_back", ((pbt1) &std::vector<std::string>::push_back))
.def("value", &string_value)
.def("__tostring", &string_vector_tostring)
REGISTER_LUA_END()

REGISTER_LUA_CLASS_NAME(std::vector<double>, VectorDouble)
.def(luabind::constructor<>())
.def("size", &std::vector<double>::size)
.def("push_back", ((pbt2) &std::vector<double>::push_back))
.def("value", &double_value)
.def("__tostring", &double_vector_tostring)
REGISTER_LUA_END()

REGISTER_LUA_CLASS_NAME(std::vector<int>, VectorInt)
.def(luabind::constructor<>())
.def("size", &std::vector<int>::size)
.def("push_back", ((pbt3) &std::vector<int>::push_back))
.def("value", &int_value)
.def("__tostring", &int_vector_tostring)
REGISTER_LUA_END()

//-----------------------------------------------------------------------
/// Declare Unit registration here since there are template
/// expansion difficulties doing it inside the Unit .cc file itself
//-----------------------------------------------------------------------

REGISTER_LUA_CLASS(Unit)
.def(luabind::constructor<const std::string&>())
.property("name", 
	  (const std::string&(Unit::*)() const) &Unit::name,
	  (void(Unit::*)(const std::string&)) &Unit::name)
.def("conversion_to_si", &Unit::conversion_to_si)
.def("is_commensurate", &Unit::is_commensurate)
REGISTER_LUA_END()

//-----------------------------------------------------------------------
/// Place all the classes that have registered with Lua into the given
/// instance of Lua.
//-----------------------------------------------------------------------

int RegisterLua::add_file_and_line(lua_State* L)
{
   lua_Debug d;
   lua_getstack(L, 1, &d);
   lua_getinfo(L, "Sln", &d);
   std::string err = lua_tostring(L, -1);
   lua_pop(L, 1);
   std::stringstream msg;
   msg << d.short_src << ":" << d.currentline;

   if (d.name != 0)
   {
      msg << "(" << d.namewhat << " " << d.name << ")";
   }
   msg << " " << err;
   lua_pushstring(L, msg.str().c_str());
   return 1;
}

namespace FullPhysics {
void RegisterLua::register_lua(lua_State* ls)
{
  // There order is not important here, *except* that you need to make
  // sure that if the Lua interface for class B refers to class A
  // (e.g., is derived from it, has a function listed for Lua that
  // uses A), then B appears after A in the list.
  REGISTER_LUA_LIST(Unit);
  REGISTER_LUA_LIST(Range);
  REGISTER_LUA_LIST(Blitz_double_array_1d);
  REGISTER_LUA_LIST(Blitz_double_array_2d);
  REGISTER_LUA_LIST(Blitz_double_array_3d);
  REGISTER_LUA_LIST(Blitz_double_array_4d);
  REGISTER_LUA_LIST(Blitz_double_array_5d);
  REGISTER_LUA_LIST(Blitz_int_array_1d);
  REGISTER_LUA_LIST(Blitz_int_array_2d);
  REGISTER_LUA_LIST(Blitz_int_array_3d);
  REGISTER_LUA_LIST(Blitz_bool_array_1d);
  REGISTER_LUA_LIST(Blitz_bool_array_2d);
  REGISTER_LUA_LIST(Blitz_bool_array_3d);
  REGISTER_LUA_LIST(VectorString);
  REGISTER_LUA_LIST(VectorDouble);
  REGISTER_LUA_LIST(VectorInt);
  REGISTER_LUA_LIST(LogImp);
  REGISTER_LUA_LIST(Logger);
  REGISTER_LUA_LIST(Constant);
  REGISTER_LUA_LIST(DefaultConstant);
  REGISTER_LUA_LIST(OldConstant);
  REGISTER_LUA_LIST(LuaCallback);
  REGISTER_LUA_LIST(DoubleWithUnit);
  REGISTER_LUA_LIST(ArrayWithUnit_1d);
  REGISTER_LUA_LIST(ArrayWithUnit_2d);
  REGISTER_LUA_LIST(ArrayWithUnit_3d);
  REGISTER_LUA_LIST(Pressure);
  REGISTER_LUA_LIST(StokesCoefficient);
  REGISTER_LUA_LIST(Temperature);
  REGISTER_LUA_LIST(Ground);
  REGISTER_LUA_LIST(Absorber);
  REGISTER_LUA_LIST(RtAtmosphere);
  REGISTER_LUA_LIST(AerosolProperty);
  REGISTER_LUA_LIST(VectorAerosolProperty);
  REGISTER_LUA_LIST(AerosolExtinction);
  REGISTER_LUA_LIST(VectorAerosolExtinction);
  REGISTER_LUA_LIST(Aerosol);
  REGISTER_LUA_LIST(AerosolOptical);
  REGISTER_LUA_LIST(MerraAerosol);
  REGISTER_LUA_LIST(AerosolMetPrior);
  REGISTER_LUA_LIST(Altitude);
  REGISTER_LUA_LIST(VectorAltitude);
  REGISTER_LUA_LIST(AbsorberVmr);
  REGISTER_LUA_LIST(VectorAbsorberVmr);
  REGISTER_LUA_LIST(GasAbsorption);
  REGISTER_LUA_LIST(VectorGasAbsorption);
  REGISTER_LUA_LIST(Absco);
  REGISTER_LUA_LIST(RadiativeTransfer);
  REGISTER_LUA_LIST(SpectralBound);
  REGISTER_LUA_LIST(SpectralWindow);
  REGISTER_LUA_LIST(NoiseModel);
  REGISTER_LUA_LIST(HdfSoundingId);
  REGISTER_LUA_LIST(Level1b);
  REGISTER_LUA_LIST(VectorLevel1b);
  REGISTER_LUA_LIST(Ils);
  REGISTER_LUA_LIST(VectorIls);
  REGISTER_LUA_LIST(Instrument);
  REGISTER_LUA_LIST(IlsFunction);
  REGISTER_LUA_LIST(InstrumentCorrection);
  REGISTER_LUA_LIST(VectorInstrumentCorrection);
  REGISTER_LUA_LIST(VectorVectorInstrumentCorrection);
  REGISTER_LUA_LIST(SpectrumEffect);
  REGISTER_LUA_LIST(VectorSpectrumEffect);
  REGISTER_LUA_LIST(VectorVectorSpectrumEffect);
  REGISTER_LUA_LIST(SpectrumEffectImpBase);
  REGISTER_LUA_LIST(FluorescenceEffect);
  REGISTER_LUA_LIST(Dispersion);
  REGISTER_LUA_LIST(VectorDispersion);
  REGISTER_LUA_LIST(IlsTableLinear);
  REGISTER_LUA_LIST(IlsTableLog);
  REGISTER_LUA_LIST(DispersionPolynomial);
  REGISTER_LUA_LIST(IlsConvolution);
  REGISTER_LUA_LIST(IlsInstrument);
  REGISTER_LUA_LIST(IlsFts);
  REGISTER_LUA_LIST(SpectrumSampling);
  REGISTER_LUA_LIST(ForwardModel);
  REGISTER_LUA_LIST(ConvergenceCheck);
  REGISTER_LUA_LIST(SolarModel);
  REGISTER_LUA_LIST(SolarAbsorptionSpectrum);
  REGISTER_LUA_LIST(SolarContinuumSpectrum);
  REGISTER_LUA_LIST(SolarDopplerShift);
  REGISTER_LUA_LIST(SolarAbsorptionAndContinuum);
  REGISTER_LUA_LIST(SolarDopplerShiftPolynomial);
  REGISTER_LUA_LIST(SolarDopplerShiftL1b);
  REGISTER_LUA_LIST(SolarAbsorptionOcoFile);
  REGISTER_LUA_LIST(SolarAbsorptionTable);
  REGISTER_LUA_LIST(SolarContinuumPolynomial);
  REGISTER_LUA_LIST(SolarContinuumTable);
  REGISTER_LUA_LIST(CostFunction);
  REGISTER_LUA_LIST(RegisterOutputBase);
  REGISTER_LUA_LIST(VectorRegisterOutput);
  REGISTER_LUA_LIST(SpectralRange);
  REGISTER_LUA_LIST(SpectralDomain);
  REGISTER_LUA_LIST(Spectrum);
  REGISTER_LUA_LIST(FpLogger);
  REGISTER_LUA_LIST(HeritageFile);
  REGISTER_LUA_LIST(AcosSoundingId);
  REGISTER_LUA_LIST(OcoSoundingId);
  REGISTER_LUA_LIST(UqSoundingId);
  REGISTER_LUA_LIST(VectorHdfSoundingId);
  REGISTER_LUA_LIST(AerosolAodOutput);
  REGISTER_LUA_LIST(AerosolParamOutput);
  REGISTER_LUA_LIST(AerosolConsolidatedOutput);
  REGISTER_LUA_LIST(Level1bHdf);
  REGISTER_LUA_LIST(Level1bHeritage);
  REGISTER_LUA_LIST(Level1bAcos);
  REGISTER_LUA_LIST(Level1bAverage);
  REGISTER_LUA_LIST(Level1bScaleRadiance);
  REGISTER_LUA_LIST(FtsRunLog);
  REGISTER_LUA_LIST(FtsRunLogRecord);
  REGISTER_LUA_LIST(FtsRunLogVector);
  REGISTER_LUA_LIST(FtsRunLogOutput);
  REGISTER_LUA_LIST(Level1bFts);
  REGISTER_LUA_LIST(Level1bOco);
  REGISTER_LUA_LIST(Level1bUq);
  REGISTER_LUA_LIST(Meteorology);
  REGISTER_LUA_LIST(AcosMetFile);
  REGISTER_LUA_LIST(OcoMetFile);
  REGISTER_LUA_LIST(UqEcmwf);
  REGISTER_LUA_LIST(OcoSimMetEcmwf);
  REGISTER_LUA_LIST(HdfFile);
  REGISTER_LUA_LIST(HdfConstant);
  REGISTER_LUA_LIST(PressureLevelInput);
  REGISTER_LUA_LIST(PressureFixedLevel);
  REGISTER_LUA_LIST(PressureFixedLevelOutput);
  REGISTER_LUA_LIST(PressureOutput);
  REGISTER_LUA_LIST(PressureSigma);
  REGISTER_LUA_LIST(AltitudeOutput);
  REGISTER_LUA_LIST(StokesCoefficientConstant);
  REGISTER_LUA_LIST(StokesCoefficientFraction);
  REGISTER_LUA_LIST(TemperatureFixedLevel);
  REGISTER_LUA_LIST(TemperatureMet);
  REGISTER_LUA_LIST(TemperatureMetOutput);
  REGISTER_LUA_LIST(TemperatureLevel);
  REGISTER_LUA_LIST(TemperatureLevelOffset);
  REGISTER_LUA_LIST(TemperatureLevelOffsetOutput);
  REGISTER_LUA_LIST(TemperatureFixedLevelOutput);
  REGISTER_LUA_LIST(RelativeHumidity);
  REGISTER_LUA_LIST(InitialGuess);
  REGISTER_LUA_LIST(InitialGuessBuilder);
  REGISTER_LUA_LIST(CompositeInitialGuess);
  REGISTER_LUA_LIST(InitialGuessValue);
  REGISTER_LUA_LIST(AerosolPropertyHdf);
  REGISTER_LUA_LIST(AerosolPropertyRhHdf);
  REGISTER_LUA_LIST(AerosolExtinctionLinear);
  REGISTER_LUA_LIST(AerosolExtinctionLog);
  REGISTER_LUA_LIST(AerosolShapeGaussian);
  REGISTER_LUA_LIST(GroundBrdfWeight);
  REGISTER_LUA_LIST(GroundBrdfWeightOutput);
  REGISTER_LUA_LIST(GroundLambertian);
  REGISTER_LUA_LIST(GroundLambertianOutput);
  REGISTER_LUA_LIST(GroundCoxmunk);
  REGISTER_LUA_LIST(GroundCoxmunkOutput);
  REGISTER_LUA_LIST(GroundCoxmunkPlusLambertian);
  REGISTER_LUA_LIST(GroundCoxmunkPlusLambertianOutput);
  REGISTER_LUA_LIST(GroundCoxmunkScaled);
  REGISTER_LUA_LIST(GroundCoxmunkScaledOutput);
  REGISTER_LUA_LIST(GroundBrdfVeg);
  REGISTER_LUA_LIST(GroundBrdfSoil);
  REGISTER_LUA_LIST(GroundBrdfOutput);
  REGISTER_LUA_LIST(AltitudeHydrostatic);
  REGISTER_LUA_LIST(AbsorberVmrFixedLevel);
  REGISTER_LUA_LIST(AbsorberVmrFixedLevelOutput);
  REGISTER_LUA_LIST(AbsorberVmrFixedLevelScaled);
  REGISTER_LUA_LIST(AbsorberVmrFixedLevelScaledOutput);
  REGISTER_LUA_LIST(AbsorberVmrMet);
  REGISTER_LUA_LIST(AbsorberVmrMetOutput);
  REGISTER_LUA_LIST(AbsorberVmrLevel);
  REGISTER_LUA_LIST(AbsorberVmrLevelOutput);
  REGISTER_LUA_LIST(AbsorberVmrLogLevel);
  REGISTER_LUA_LIST(AbsorberVmrLogLevelOutput);
  REGISTER_LUA_LIST(AbsorberVmrLevelScaled);
  REGISTER_LUA_LIST(AbsorberVmrLevelScaledOutput);
  REGISTER_LUA_LIST(AbscoHdf);
  REGISTER_LUA_LIST(AbsorberAbsco);
  REGISTER_LUA_LIST(AbsorberAbscoOutput);
  REGISTER_LUA_LIST(AtmosphereOco);
  REGISTER_LUA_LIST(LidortRt);
  REGISTER_LUA_LIST(SpectralWindowRange);
  REGISTER_LUA_LIST(LRadRt);
  REGISTER_LUA_LIST(LsiRt);
  REGISTER_LUA_LIST(TwostreamRt);
  REGISTER_LUA_LIST(HresWrapper);
  REGISTER_LUA_LIST(UplookingRaytracing);
  REGISTER_LUA_LIST(ChapmanBoaRT);
  REGISTER_LUA_LIST(PrecomputedNoiseModel);
  REGISTER_LUA_LIST(GosatNoiseModel);
  REGISTER_LUA_LIST(OcoNoiseModel);
  REGISTER_LUA_LIST(UqNoiseModel);
  REGISTER_LUA_LIST(BadSampleNoiseModel);
  REGISTER_LUA_LIST(SpectrallyResolvedNoise);
  REGISTER_LUA_LIST(SpectrumSamplingFixedSpacing);
  REGISTER_LUA_LIST(StateVector);
  REGISTER_LUA_LIST(OcoForwardModel);
  REGISTER_LUA_LIST(ForwardModelCostFunction);
  REGISTER_LUA_LIST(ConnorConvergence);
  REGISTER_LUA_LIST(ConnorSolver);
  REGISTER_LUA_LIST(ErrorAnalysis);
  REGISTER_LUA_LIST(StateVectorOutput);
  REGISTER_LUA_LIST(Level1bOutput);
  REGISTER_LUA_LIST(DispersionPolynomialOutput);
  REGISTER_LUA_LIST(IlsTableLinearOutput);
  REGISTER_LUA_LIST(IlsTableLogOutput);
  REGISTER_LUA_LIST(StokesCoefficientFractionOutput);
  REGISTER_LUA_LIST(ErrorAnalysisOutput);
  REGISTER_LUA_LIST(ForwardModelOutput);
  REGISTER_LUA_LIST(OcoForwardModelOutput);
  REGISTER_LUA_LIST(ConnorConvergenceOutput);
  REGISTER_LUA_LIST(ConnorSolverOutput);
  REGISTER_LUA_LIST(ZeroOffsetWaveform);
  REGISTER_LUA_LIST(ZeroOffsetWaveformOutput);
  REGISTER_LUA_LIST(EmpiricalOrthogonalFunction);
  REGISTER_LUA_LIST(EmpiricalOrthogonalFunctionOutput);
  REGISTER_LUA_LIST(RadianceScalingSvFit);
  REGISTER_LUA_LIST(RadianceScalingLinearFit);
  REGISTER_LUA_LIST(NonuniformSpectrumSampling);
  REGISTER_LUA_LIST(TcconApriori);
  REGISTER_LUA_LIST(CO2ProfilePrior);
  REGISTER_LUA_LIST(GasVmrApriori);
  REGISTER_LUA_LIST(GasVmrAprioriOutput);
  REGISTER_LUA_LIST(OcoSimApriori);
  REGISTER_LUA_LIST(DispersionFit);
  REGISTER_LUA_LIST(DispersionFitOutput);
  REGISTER_LUA_LIST(SolverIterationLog);
  REGISTER_LUA_LIST(SourceFilesOutput);
  REGISTER_LUA_LIST(MetPassThroughOutput);
  REGISTER_LUA_LIST(HighResSpectrumOutput);
  REGISTER_LUA_LIST(MaxAPosteriori);
  REGISTER_LUA_LIST(MaxLikelihood);
  REGISTER_LUA_LIST(MaxAPosterioriOCO);
  REGISTER_LUA_LIST(MaxLikelihoodOCO);
  REGISTER_LUA_LIST(CostFunc);
  REGISTER_LUA_LIST(CostFuncDiff);
  REGISTER_LUA_LIST(NLLSProblem);
  REGISTER_LUA_LIST(NLLSProblemScaled);
  REGISTER_LUA_LIST(NLLSMaxLikelihood);
  REGISTER_LUA_LIST(NLLSMaxAPosteriori);
  REGISTER_LUA_LIST(IterativeSolver);
  REGISTER_LUA_LIST(IterativeSolverDer);
  REGISTER_LUA_LIST(NLLSSolver);
  REGISTER_LUA_LIST(NLLSSolverGSL);
  REGISTER_LUA_LIST(NLLSSolverGSLLMDER);
  REGISTER_LUA_LIST(NLLSSolverGSLLMSDER);
  REGISTER_LUA_LIST(ConnorSolverMAP);
  REGISTER_LUA_LIST(CostMinimizer);
  REGISTER_LUA_LIST(CostMinimizerGSL);
  REGISTER_LUA_LIST(MaxAPosterioriOutput);
  REGISTER_LUA_LIST(InstrumentDoppler);
  REGISTER_LUA_LIST(SolarAbsorptionGfitFile);
  REGISTER_LUA_LIST(FluorescenceEffectOutput);
  REGISTER_LUA_LIST(RadianceScalingOutput);
}
}

extern "C" {
  int luaopen_libfull_physics(lua_State* ls)
  {
    luabind::open(ls);
    RegisterLua::register_lua(ls);
    luabind::set_pcall_callback(RegisterLua::add_file_and_line);
    return 1;
  }
}
