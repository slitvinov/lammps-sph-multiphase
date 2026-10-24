#ifdef PAIR_CLASS
PairStyle(sph/heatconduction/multiphase,PairSPHHeatConductionMultiphase);
#else

#ifndef LMP_PAIR_SPH_HEATCONDUCTION_MULTIPHASE_H
#define LMP_PAIR_SPH_HEATCONDUCTION_MULTIPHASE_H

#include "pair.h"

namespace LAMMPS_NS {

class PairSPHHeatConductionMultiphase : public Pair {
 public:
  PairSPHHeatConductionMultiphase(class LAMMPS *);
  ~PairSPHHeatConductionMultiphase() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  void init_style() override;

 protected:
  double **cut, **alpha;
  void allocate();
};

}

#endif
#endif
