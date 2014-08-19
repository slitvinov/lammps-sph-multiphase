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
#include "fix_setmeso.h"
#include "sph_energy_equation.h"
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

FixSetMeso::FixSetMeso(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal fix setmeso command");

  global_freq = 1;
  xstr = NULL;
  rhoflag = eflag = tflag = 0;
  if (strcmp(arg[3],"meso_rho")==0) {
    rhoflag = 1;
  } else if (strcmp(arg[3],"meso_e")==0) {
    eflag = 1;
  } else if (strcmp(arg[3],"meso_t")==0) {
    tflag = 1;
  }
  else {
    error->all(FLERR,"Illegal fix setmeso command, meso_rho or meso_e must be given");
  }

  if (strstr(arg[4],"v_") == arg[4]) {
    int n = strlen(&arg[4][2]) + 1;
    xstr = new char[n];
    strcpy(xstr,&arg[4][2]);
  } else {
    xvalue = atof(arg[4]);
    xstyle = CONSTANT;
  }
  // optional args

  iregion = -1;
  idregion = NULL;
  regionflag = 1;

  int iarg = 5;
  while (iarg < narg) {
    if ( (strcmp(arg[iarg],"region") == 0) || (strcmp(arg[iarg],"noregion") == 0) ) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix setmesode command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR,"Region ID for fix setmesode does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      if (strcmp(arg[iarg],"noregion") == 0) {
	regionflag = 0;
      } 
      iarg += 2;
    } else error->all(FLERR,"Illegal fix setmesode command");
  }

  force_flag = 0;
  mesovarorg = 0.0;

  maxatom = 0;
  smesovar = NULL;
}

/* ---------------------------------------------------------------------- */

FixSetMeso::~FixSetMeso()
{
  delete [] xstr;
  delete [] idregion;
  memory->destroy(smesovar);
}

/* ---------------------------------------------------------------------- */

int FixSetMeso::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSetMeso::init()
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

void FixSetMeso::setup(int vflag)
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

void FixSetMeso::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSetMeso::post_force(int vflag)
{
  double **x = atom->x;
  double *cv = atom->cv;
  double *rho = atom->rho;
  double *mesovar;
  if (rhoflag) {
    mesovar = atom->rho;
  } else if (eflag) {
    mesovar = atom->e;
  } else if (tflag) {
    // also use energy
    mesovar = atom->e;
  } else {
    error->all(FLERR,"Illegal fix setmeso command");
  }

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // reallocate smesovar array if necessary

  if (varflag == ATOM && nlocal > maxatom) {
    maxatom = atom->nmax;
    memory->destroy(smesovar);
    memory->create(smesovar,maxatom,"setmesode:smesovar");
  }

  mesovarorg = 0.0;
  force_flag = 0;

  if (varflag == CONSTANT) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	if (iregion >= 0) {
	  if ( (regionflag) && 
	      !domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2]))
          continue;
	  /// applay fix outside of the region
	  if ( (!regionflag) && 
	      domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2]))
          continue;
	}

        mesovarorg += mesovar[i];
        if (xstyle) {
	  if (tflag) {
	    // xvalue is temperature and we set energy
	    mesovar[i] = sph_t2energy(xvalue,cv[i]);
	  } else {
	    mesovar[i] = xvalue;
	  }
	}
      }

  // variable force, wrap with clear/add

  } else {

    modify->clearstep_compute();

    if (xstyle == EQUAL) xvalue = input->variable->compute_equal(xvar);
    else if (xstyle == ATOM && smesovar)
      input->variable->compute_atom(xvar,igroup,&smesovar[0],1,0);

    modify->addstep_compute(update->ntimestep + 1);

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (iregion >= 0 &&
            !domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2]))
          continue;

        mesovarorg += mesovar[i];
        if (xstyle == ATOM) {
	  if (tflag) {
	    mesovar[i] = sph_t2energy(smesovar[i],cv[i]);
	  } else {
	    mesovar[i] = smesovar[i];	    
	  }
	}
        else if (xstyle) {
	  if (tflag) {
	    // xvalue is temperature and we set energy
	    mesovar[i] = sph_t2energy(xvalue,cv[i]);
	  } else {
	    mesovar[i] = xvalue;
	  }
	}
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixSetMeso::post_force_respa(int vflag, int ilevel, int iloop)
{
  // set force to desired value on outermost level, 0.0 on other levels

  if (ilevel == nlevels_respa-1) post_force(vflag);
  else {
    
    double *mesovar;
    if (rhoflag) {
      mesovar = atom->rho;
    } else if (eflag) {
      mesovar = atom->e;
    } else {
      error->all(FLERR,"Illegal fix setmeso command");
    }
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (xstyle) mesovar[i] = 0.0;
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixSetMeso::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   return components of total force on fix group before force was changed
------------------------------------------------------------------------- */

double FixSetMeso::compute_scalar()
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(&mesovarorg,&mesovarorg_all,1,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return mesovarorg_all;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixSetMeso::memory_usage()
{
  double bytes = 0.0;
  if (varflag == ATOM) bytes = atom->nmax*3 * sizeof(double);
  return bytes;
}
