#ifdef PAIR_CLASS
PairStyle(sph/heatconduction/phasechange,PairSPHHeatConductionPhaseChange);
#else

#ifndef LMP_PAIR_SPH_HEATCONDUCTION_PHASECHANGE_H
#define LMP_PAIR_SPH_HEATCONDUCTION_PHASECHANGE_H

#include "pair.h"

namespace LAMMPS_NS {

class PairSPHHeatConductionPhaseChange : public Pair {
 public:
  PairSPHHeatConductionPhaseChange(class LAMMPS *);
  ~PairSPHHeatConductionPhaseChange() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  void init_style() override;

 protected:
  double **cut, **alpha, **tc;
  int **fixflag;
  void allocate();
};

}

#endif
#endif
