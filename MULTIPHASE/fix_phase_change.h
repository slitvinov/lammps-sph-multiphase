#ifdef FIX_CLASS
FixStyle(phase_change,FixPhaseChange);
#else

#ifndef LMP_FIX_PHASE_CHANGE_H
#define LMP_FIX_PHASE_CHANGE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPhaseChange : public Fix {
 public:
  FixPhaseChange(class LAMMPS *, int, char **);
  ~FixPhaseChange() override;
  int setmask() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void pre_exchange() override;
  void write_restart(FILE *) override;
  void restart(char *) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

 private:
  int from_type, to_type, nfreq, seed, maxattempt, nfirst;
  bool energy_chance_flag;
  int convert_flag;
  double convert_cv;
  double Tc, Tt, Hwv, dr, to_mass, cutoff, change_chance, phase_change_rate;
  char *idregion;
  double *dmass;
  class NeighList *list;
  class RanPark *random;

  bool owns(double *, double *, double *);
  bool insert_one_atom(double *, double *, double *);
  void create_newpos(double *, double *, double, double *);
  void create_newpos_simple(double *, double, double *);
  bool isfromphasearound(int);
};

}

#endif
#endif
