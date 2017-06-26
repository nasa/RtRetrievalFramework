#ifndef LEVEL_1B_UQ_H
#define LEVEL_1B_UQ_H
#include "level_1b_oco.h"

namespace FullPhysics {

/****************************************************************//**
  This reads a Uncertainty Quantification L1B file.
*******************************************************************/
class Level1bUq: public Level1bOco {
public:

    Level1bUq(const boost::shared_ptr<HdfFile>& Hfile,
              const boost::shared_ptr<HdfSoundingId>& Sounding_id);

    void set_radiance(int Spec_index, boost::shared_ptr<SpectralRange>& Rad);

protected:
    virtual SpectralRange radiance_no_uncertainty(int Spec_index) const;

private:
    void initialize();

    std::vector<boost::shared_ptr<SpectralRange> > radiance;
};
}

#endif
