#ifdef PAIR_CLASS
PairStyle(sph/surfacetension,PairSPHSurfaceTension);
#else

#ifndef LMP_PAIR_SPH_SURFACETENSION_H
#define LMP_PAIR_SPH_SURFACETENSION_H

#include "pair.h"

namespace LAMMPS_NS {

class PairSPHSurfaceTension : public Pair {
 public:
  PairSPHSurfaceTension(class LAMMPS *);
  ~PairSPHSurfaceTension() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  void init_style() override;

 protected:
  double **cut;
  void allocate();
};

}

#endif
#endif
