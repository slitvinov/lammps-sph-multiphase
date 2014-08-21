/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "string.h"
#include "stdlib.h"
#include "fix_setmesode.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

FixSetMesodE::FixSetMesodE(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal fix setmesode command");

  vector_flag = 1;
  size_vector = 1;
  global_freq = 1;
  extvector = 1;

  xstr = NULL;

  if (strstr(arg[3],"v_") == arg[3]) {
    int n = strlen(&arg[3][2]) + 1;
    xstr = new char[n];
    strcpy(xstr,&arg[3][2]);
  } else if (strcmp(arg[3],"NULL") == 0) {
    xstyle = NONE;
  } else {
    xvalue = atof(arg[3]);
    xstyle = CONSTANT;
  }

  // optional args

  iregion = -1;
  idregion = NULL;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix setmesode command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR,"Region ID for fix setmesode does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix setmesode command");
  }

  force_flag = 0;
  deoriginal = 0.0;

  maxatom = 0;
  sde = NULL;
}

/* ---------------------------------------------------------------------- */

FixSetMesodE::~FixSetMesodE()
{
  delete [] xstr;
  delete [] idregion;
  memory->destroy(sde);
}

/* ---------------------------------------------------------------------- */

int FixSetMesodE::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSetMesodE::init()
{
  // check variables

  if (xstr) {
    xvar = input->variable->find(xstr);
    if (xvar < 0)
      error->all(FLERR,"Variable name for fix setmesode does not exist");
    if (input->variable->equalstyle(xvar)) xstyle = EQUAL;
    else if (input->variable->atomstyle(xvar)) xstyle = ATOM;
    else error->all(FLERR,"Variable for fix setmesode is invalid style");
  }

  // set index and check validity of region

  if (iregion >= 0) {
    iregion = domain->find_region(idregion);
    if (iregion == -1)
      error->all(FLERR,"Region ID for fix setmesode does not exist");
  }

  if (xstyle == ATOM )
    varflag = ATOM;
  else if (xstyle == EQUAL)
    varflag = EQUAL;
  else varflag = CONSTANT;

  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  // cannot use non-zero forces for a minimization since no energy is integrated
  // use fix addforce instead

  int flag = 0;
  if (update->whichflag == 2) {
    if (xstyle == EQUAL || xstyle == ATOM) flag = 1;
    if (xstyle == CONSTANT && xvalue != 0.0) flag = 1;
  }
  if (flag)
    error->all(FLERR,"Cannot use non-zero forces in an energy minimization");
}

/* ---------------------------------------------------------------------- */

void FixSetMesodE::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else
    for (int ilevel = 0; ilevel < nlevels_respa; ilevel++) {
      ((Respa *) update->integrate)->copy_flevel_f(ilevel);
      post_force_respa(vflag,ilevel,0);
      ((Respa *) update->integrate)->copy_f_flevel(ilevel);
    }
}

/* ---------------------------------------------------------------------- */

void FixSetMesodE::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSetMesodE::post_force(int vflag)
{
  double **x = atom->x;
  double *de = atom->de;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // reallocate sde array if necessary

  if (varflag == ATOM && nlocal > maxatom) {
    maxatom = atom->nmax;
    memory->destroy(sde);
    memory->create(sde,maxatom,"setmesode:sde");
  }

  deoriginal = 0.0;
  force_flag = 0;

  if (varflag == CONSTANT) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (iregion >= 0 &&
            !domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2]))
          continue;

        deoriginal += de[i];
        if (xstyle) de[i] = xvalue;
      }

  // variable force, wrap with clear/add

  } else {

    modify->clearstep_compute();

    if (xstyle == EQUAL) xvalue = input->variable->compute_equal(xvar);
    else if (xstyle == ATOM && sde)
      input->variable->compute_atom(xvar,igroup,&sde[0],1,0);

    modify->addstep_compute(update->ntimestep + 1);

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (iregion >= 0 &&
            !domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2]))
          continue;

        deoriginal += de[i];
        if (xstyle == ATOM) de[i] = sde[i];
        else if (xstyle) de[i] = xvalue;
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixSetMesodE::post_force_respa(int vflag, int ilevel, int iloop)
{
  // set force to desired value on outermost level, 0.0 on other levels

  if (ilevel == nlevels_respa-1) post_force(vflag);
  else {
    double *de = atom->de;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (xstyle) de[i] = 0.0;
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixSetMesodE::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   return components of total force on fix group before force was changed
------------------------------------------------------------------------- */

double FixSetMesodE::compute_scalar()
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(&deoriginal,&deoriginal_all,1,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return deoriginal_all;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixSetMesodE::memory_usage()
{
  double bytes = 0.0;
  if (varflag == ATOM) bytes = atom->nmax*3 * sizeof(double);
  return bytes;
}
