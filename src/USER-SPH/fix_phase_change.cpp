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
#include "assert.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_phase_change.h"
#include "sph_kernel_quintic.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "pair.h"
#include "random_park.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "comm.h"
#include "domain.h"
#include "lattice.h"
#include "region.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPhaseChange::FixPhaseChange(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  int nnarg = 13;
  if (narg < nnarg) error->all(FLERR,"Illegal fix phase_change command");

  restart_global = 1;
  time_depend = 1;

  // required args
  int m = 3;
  Tc = atof(arg[m++]);
  Tt = atof(arg[m++]);
  cp = atof(arg[m++]);
  dr = atof(arg[m++]);
  to_mass = atof(arg[m++]);
  cutoff = atof(arg[m++]);
  from_type = atoi(arg[m++]);
  to_type = atoi(arg[m++]);
  nfreq = atoi(arg[m++]);
  seed = atoi(arg[m++]);
  assert(m==nnarg);

  if (seed <= 0) error->all(FLERR,"Illegal value for seed");

  iregion = -1;
  idregion = NULL;
  maxattempt = 10;
  scaleflag = 1;

  // read options from end of input line

  options(narg-nnarg,&arg[nnarg]);

  // error checks on region and its extent being inside simulation box

  if (iregion == -1) error->all(FLERR,"Must specify a region in fix phase_change");
  if (domain->regions[iregion]->bboxflag == 0)
    error->all(FLERR,"Fix phase_change region does not support a bounding box");
  if (domain->regions[iregion]->dynamic_check())
    error->all(FLERR,"Fix phase_change region cannot be dynamic");

  xlo = domain->regions[iregion]->extent_xlo;
  xhi = domain->regions[iregion]->extent_xhi;
  ylo = domain->regions[iregion]->extent_ylo;
  yhi = domain->regions[iregion]->extent_yhi;
  zlo = domain->regions[iregion]->extent_zlo;
  zhi = domain->regions[iregion]->extent_zhi;

  if (domain->triclinic == 0) {
    if (xlo < domain->boxlo[0] || xhi > domain->boxhi[0] ||
        ylo < domain->boxlo[1] || yhi > domain->boxhi[1] ||
        zlo < domain->boxlo[2] || zhi > domain->boxhi[2])
      error->all(FLERR,"Phase change region extends outside simulation box");
  } else {
    if (xlo < domain->boxlo_bound[0] || xhi > domain->boxhi_bound[0] ||
        ylo < domain->boxlo_bound[1] || yhi > domain->boxhi_bound[1] ||
        zlo < domain->boxlo_bound[2] || zhi > domain->boxhi_bound[2])
      error->all(FLERR,"Phase change region extends outside simulation box");
  }

  // setup scaling

  if (scaleflag && domain->lattice == NULL)
    error->all(FLERR,"Use of fix phase_change with undefined lattice");

  double xscale,yscale,zscale;
  if (scaleflag) {
    xscale = domain->lattice->xlattice;
    yscale = domain->lattice->ylattice;
    zscale = domain->lattice->zlattice;
  }
  else xscale = yscale = zscale = 1.0;

  random = new RanPark(lmp,seed);
  // set up reneighboring

  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;
  nfirst = next_reneighbor;
}

/* ---------------------------------------------------------------------- */

FixPhaseChange::~FixPhaseChange()
{
  delete random;
  delete [] idregion;
}

/* ---------------------------------------------------------------------- */

int FixPhaseChange::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPhaseChange::init()
{
  // set index and check validity of region

  iregion = domain->find_region(idregion);
  if (iregion == -1)
    error->all(FLERR,"Region ID for fix phase_change does not exist");

  // need a full neighbor list, built whenever re-neighboring occurs

  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

void FixPhaseChange::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ----------------------------------------------------------------------
   perform phase change
------------------------------------------------------------------------- */

void FixPhaseChange::pre_exchange()
{
  // just return if should not be called on this timestep

  if (next_reneighbor != update->ntimestep) return;

  double *sublo,*subhi;
  if (domain->triclinic == 0) {
    sublo = domain->sublo;
    subhi = domain->subhi;
  } else {
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  }

  // attempt an insertion until successful


    // choose random position for new atom within region
  int nins = 0;
  int nlocal = atom->nlocal;
  int inum = list->inum;
  int* ilist = list->ilist;
  int* numneigh = list->numneigh;
  double **x = atom->x;
  double **v = atom->v;
  double **cg = atom->colorgradient;
  double *rmass = atom->rmass;
  double *e   = atom->e;
  double *de   = atom->de;
  int *type = atom->type;
  
  /// TODO: clear de
  int nall;
  if (force->newton) nall = atom->nlocal + atom->nghost;
  else nall = atom->nlocal;
  for (int i = 0; i < nall; i++) {
    de[i] = 0.0;
  }

  /// TODO: how to distribute to ghosts?
  for (int i = 0; i < nlocal; i++) {
    // clear de, I will use it to transfer energy
    double abscgi = sqrt(cg[i][0]*cg[i][0] +
			 cg[i][1]*cg[i][1] +
			 cg[i][2]*cg[i][2]);
    if ( (abscgi>1e-20) && (e[i]>Tt) && (type[i] == to_type) ) {
      double coord[3];
      // place an atom in the directions of color gradient
      double eij[3];
      // eij[0] = random->uniform() - 0.5;
      // eij[1] = random->uniform() - 0.5;
      // eij[2] = random->uniform() - 0.5;
      eij[0] = -cg[i][1];
      eij[1] = cg[i][0];
      eij[2] = 0.0;
      double eijabs = sqrt(eij[0]*eij[0] + eij[1]*eij[1] + eij[2]*eij[2]);
      coord[0] = x[i][0] + eij[0]*dr/eijabs;
      coord[1] = x[i][1] + eij[1]*dr/eijabs;
      coord[2] = x[i][2] + eij[2]*dr/eijabs;
      bool ok = insert_one_atom(coord, sublo, subhi);
      if (ok) {
	/// distribute energy to neighboring particles
	/// latent heat + heat particle i + heat a new particle
	double energy_to_dist = cp + (Tc - e[i] + Tc)*rmass[i]/to_mass;
	nins++;
	// look for the neighbors of the type from_type
	// and extract energy from all of them
	double xtmp = x[i][0];
	double ytmp = x[i][1];
	double ztmp = x[i][2];
	int** firstneigh = list->firstneigh;
	int jnum = numneigh[i];
	int* jlist = firstneigh[i];
	// collect
	double wtotal = 0.0;
	  /// TODO: make it run in parallel
	for (int jj = 0; jj < jnum; jj++) {
	  int j = jlist[jj];
	  j &= NEIGHMASK;
	  if (type[j]==from_type) {
	    double delx = xtmp - x[j][0];
	    double dely = ytmp - x[j][1];
	    double delz = ztmp - x[j][2];
	    double rsq = delx * delx + dely * dely + delz * delz;
	    double wfd;
	    if (domain->dimension == 3) {
	      wfd = sph_dw_quintic2d(sqrt(rsq)*cutoff);
	    } else {
	      wfd = sph_dw_quintic3d(sqrt(rsq)*cutoff);
	    }
	    wtotal+=wfd;
	  }
	}

	// distribute
	for (int jj = 0; jj < jnum; jj++) {
	  int j = jlist[jj];
	  j &= NEIGHMASK;
	  if (type[j]==from_type) {
	    double delx = xtmp - x[j][0];
	    double dely = ytmp - x[j][1];
	    double delz = ztmp - x[j][2];
	    double rsq = delx * delx + dely * dely + delz * delz;
	    double wfd;
	    if (domain->dimension == 3) {
	      wfd = sph_dw_quintic2d(sqrt(rsq)*cutoff);
	    } else {
	      wfd = sph_dw_quintic3d(sqrt(rsq)*cutoff);
	    }
	    de[j] -= energy_to_dist*wfd/wtotal;
	  }
	}
	
	e[i] = Tc;
	// for a new atom
	rmass[atom->nlocal-1] = to_mass;
	e[atom->nlocal-1] = Tc;
	de[atom->nlocal-1] = 0.0;
	v[atom->nlocal-1][0] = v[i][0];
	v[atom->nlocal-1][1] = v[i][1];
	v[atom->nlocal-1][2] = v[i][2];
      }
    }
  }

  /// find a total number of inserted atoms
  int ninsall;
  MPI_Allreduce(&nins,&ninsall,1,MPI_INT,MPI_SUM,world);
  next_reneighbor += nfreq;
  
  // reset global natoms
  // set tag # of new particle beyond all previous atoms
  // if global map exists, reset it now instead of waiting for comm
  // since deleting atoms messes up ghosts
  int success = 1;
  if (success) {
    atom->natoms += ninsall;
    if (atom->tag_enable) {
      atom->tag_extend();
      if (atom->map_style) {
        atom->nghost = 0;
        atom->map_init();
        atom->map_set();
      }
    }
  }

  MPI_Allreduce(&nins,&ninsall,1,MPI_INT,MPI_SUM,world);
  comm->reverse_comm();
  for (int i = 0; i < nlocal; i++) {
    e[i] += de[i];
  }
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of input line
------------------------------------------------------------------------- */

void FixPhaseChange::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR,"Illegal fix indent command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix phase_change command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR,"Region ID for fix phase_change does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"attempt") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix phase_change command");
      maxattempt = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix phase_change command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all(FLERR,"Illegal fix phase_change command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix phase_change command");
  }
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixPhaseChange::write_restart(FILE *fp)
{
  int n = 0;
  double list[4];
  list[n++] = next_reneighbor;

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixPhaseChange::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;

  seed = static_cast<int> (list[n++]);
  nfirst = static_cast<int> (list[n++]);
  next_reneighbor = static_cast<int> (list[n++]);

  random->reset(seed);
}

bool FixPhaseChange::insert_one_atom(double* coord, double* sublo, double* subhi)
{
  int flagall, flag;
  double lamda[3];
  double *newcoord;

  int i, j;
  int nfix = modify->nfix;
  Fix **fix = modify->fix;

  int nlocal = atom->nlocal;
  int success = 1;

  // check if new atom is in my sub-box or above it if I'm highest proc
  // if so, add to my list via create_atom()
  // initialize info about the atoms
  // set group mask to "all" plus fix group
  if (domain->triclinic) {
    domain->x2lamda(coord,lamda);
    newcoord = lamda;
  } else newcoord = coord;

  flag = 0;
  if (newcoord[0] >= sublo[0] && newcoord[0] < subhi[0] &&
      newcoord[1] >= sublo[1] && newcoord[1] < subhi[1] &&
      newcoord[2] >= sublo[2] && newcoord[2] < subhi[2]) flag = 1;
  else if (domain->dimension == 3 && newcoord[2] >= domain->boxhi[2] &&
	   comm->myloc[2] == comm->procgrid[2]-1 &&
	   newcoord[0] >= sublo[0] && newcoord[0] < subhi[0] &&
	   newcoord[1] >= sublo[1] && newcoord[1] < subhi[1]) flag = 1;
  else if (domain->dimension == 2 && newcoord[1] >= domain->boxhi[1] &&
	   comm->myloc[1] == comm->procgrid[1]-1 &&
	   newcoord[0] >= sublo[0] && newcoord[0] < subhi[0]) flag = 1;

  if (flag) {
    atom->avec->create_atom(to_type,coord);
    int m = atom->nlocal - 1;
    atom->type[m] = to_type;
    atom->mask[m] = 1 | groupbit;
    for (j = 0; j < nfix; j++)
      if (fix[j]->create_attribute) fix[j]->set_arrays(m);
  }
  if (flag) {
    return true;
  } else {
    return false;
  }
}

