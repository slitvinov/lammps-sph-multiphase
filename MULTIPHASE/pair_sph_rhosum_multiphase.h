#ifdef PAIR_CLASS
PairStyle(sph/rhosum/multiphase,PairSPHRhoSumMultiphase);
#else

#ifndef LMP_PAIR_SPH_RHOSUM_MULTIPHASE_H
#define LMP_PAIR_SPH_RHOSUM_MULTIPHASE_H

#include "pair.h"

namespace LAMMPS_NS {

class PairSPHRhoSumMultiphase : public Pair {
 public:
  PairSPHRhoSumMultiphase(class LAMMPS *);
  ~PairSPHRhoSumMultiphase() override;
  void init_style() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

 protected:
  double **cut;
  int nstep;
  void allocate();
};

}

#endif
#endif
