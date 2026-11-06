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
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

 protected:
  double **cut, **alpha, **tc;
  int **fixflag, **latent, anylatent = 0, nmax;
  double *dlat;
  void allocate();
};

}

#endif
#endif
