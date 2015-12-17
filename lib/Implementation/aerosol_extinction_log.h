#ifndef AEROSOL_EXTINCTION_LOG_H
#define AEROSOL_EXTINCTION_LOG_H
#include "aerosol_extinction_imp_base.h"
#include "pressure.h"

namespace FullPhysics {
/****************************************************************//**
  This class maps the state vector to the aerosol extinction on each
  level.

  This implementation just gets the log of the extinction coefficient
  for each level from the state vector.
*******************************************************************/
class AerosolExtinctionLog : public AerosolExtinctionImpBase {
public:
//-----------------------------------------------------------------------
/// Constructor.
/// \param Press The pressure to use
/// \param Flag Boolean flag indicating which levels are to be set by
///   the state vector. A value of false means the level is held fixed
///   when the state vector changes.
/// \param Aext The aerosol extinction value.
/// \param Aerosol_name The name of the aerosol. This is used to
///   generate the state vector name metadata, so it should be
///   whatever is convenient.
//-----------------------------------------------------------------------

  AerosolExtinctionLog(const boost::shared_ptr<Pressure>& Press,
			  const blitz::Array<bool, 1>& Flag, 
			  const blitz::Array<double, 1>& Aext,
			  const std::string& Aerosol_name)
    : AerosolExtinctionImpBase(Aerosol_name, Aext, Flag, Press) {}
  virtual ~AerosolExtinctionLog() {}
  virtual boost::shared_ptr<AerosolExtinction> clone() const
  { return clone(press->clone()); }
  virtual boost::shared_ptr<AerosolExtinction> clone
  (const boost::shared_ptr<Pressure>& P) const;
  virtual std::string state_vector_name_i(int i) const
  { return "Aerosol " + aerosol_name() + " Log(extinction) for Press Lvl " +
      boost::lexical_cast<std::string>(i + 1); }
  virtual std::string model_short_name() const { return "profile_log"; }
  virtual void print(std::ostream& Os) const;
protected:
  virtual void calc_aerosol_extinction() const;
};
}
#endif
