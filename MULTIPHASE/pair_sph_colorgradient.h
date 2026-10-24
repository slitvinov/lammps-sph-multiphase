#ifdef PAIR_CLASS
PairStyle(sph/colorgradient,PairSPHColorGradient);
#else

#ifndef LMP_PAIR_SPH_COLORGRADIENT_H
#define LMP_PAIR_SPH_COLORGRADIENT_H

#include "pair.h"

namespace LAMMPS_NS {

class PairSPHColorGradient : public Pair {
 public:
  PairSPHColorGradient(class LAMMPS *);
  ~PairSPHColorGradient() override;
  void init_style() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

 protected:
  double **cut, **alpha;
  int nstep;
  void allocate();
};

}

#endif
#endif
