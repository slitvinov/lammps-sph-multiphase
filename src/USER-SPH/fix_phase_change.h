/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(phase_change,FixPhaseChange)

#else

#ifndef LMP_FIX_PHASECHANGE_H
#define LMP_FIX_PHASECHANGE_H

#include "stdio.h"
#include "fix.h"

namespace LAMMPS_NS {

class FixPhaseChange : public Fix {
 public:
  FixPhaseChange(class LAMMPS *, int, char **);
  ~FixPhaseChange();
  int setmask();
  void init();
  void init_list(int, class NeighList *);
  void pre_exchange();
  void write_restart(FILE *);
  void restart(char *);
  
 private:
  int from_type,to_type,nfreq,seed;
  int iregion,maxattempt,scaleflag;
  char *idregion;
  double xlo,xhi,ylo,yhi,zlo,zhi;
  double cutoff;
  int nfirst;
  class NeighList *list;

  // critical temperature
  double Tc;
  // temperature which triggers phase change
  double Tt;
  // latent heat of evaporation (in energy per mass)
  double Hwv;
  // how far a new atom should be placed
  double dr;
  // mass of a particle of a new phase
  double to_mass;
  /// probability of phase change [0.0:1.0]
  double change_chance;
  
  /// if true we define probability of phase change based on energy
  bool energy_chance_flag;
  double phase_change_rate;

  // mass to distribute
  double* dmass;
  void unpack_reverse_comm(int n, int *list, double *buf);
  int pack_reverse_comm(int n, int first, double *buf);

  class RanPark *random;

  void options(int, char **);
  // insert one atom, return true if sucsesfull
  bool insert_one_atom(double* coord, double* sublo, double* subhi);
  void create_newpos(double* xone, double* cgone, double delta, double* coord);
  void create_newpos_simple(double* xone, double delta, double* coord);
  bool isfromphasearound(int i);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Must specify a region in fix deposit

The region keyword must be specified with this fix.

E: Fix deposit region does not support a bounding box

Not all regions represent bounded volumes.  You cannot use
such a region with the fix deposit command.

E: Fix deposit region cannot be dynamic

Only static regions can be used with fix deposit.

E: Deposition region extends outside simulation box

Self-explanatory.

E: Use of fix deposit with undefined lattice

Must use lattice command with compute fix deposit command if units
option is set to lattice.

E: Region ID for fix deposit does not exist

Self-explanatory.

W: Particle deposition was unsuccessful

The fix deposit command was not able to insert as many atoms as
needed.  The requested volume fraction may be too high, or other atoms
may be in the insertion region.

*/
