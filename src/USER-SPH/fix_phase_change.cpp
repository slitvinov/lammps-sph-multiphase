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
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_phase_change.h"
#include "sph_kernel_quintic.h"
#include "sph_energy_equation.h"
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

#define CG_SMALL 1.0e-20

/* ---------------------------------------------------------------------- */

FixPhaseChange::FixPhaseChange(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  // communicate energy change due to phase change
  comm_reverse = 1;
  int nnarg = 14;
  if (narg < nnarg) error->all(FLERR,"Illegal fix phase_change command");

  restart_global = 1;
  time_depend = 1;

  // required args
  int m = 3;
  Tc = atof(arg[m++]);
  Tt = atof(arg[m++]);
  Hwv = atof(arg[m++]);
  dr = atof(arg[m++]);
  to_mass = atof(arg[m++]);
  cutoff = atof(arg[m++]);
  from_type = atoi(arg[m++]);
  to_type = atoi(arg[m++]);
  nfreq = atoi(arg[m++]);
  seed = atoi(arg[m++]);
  if (seed <= 0) error->all(FLERR,"Illegal value for seed");
  // chance is based on energy
  if (strcmp(arg[m++],"ENERGY") == 0) {
    energy_chance_flag = true;
    phase_change_rate = atof(arg[m++]);
    nnarg = 15;
  } else {
    energy_chance_flag = false;
    change_chance = atof(arg[m-1]);
    if (change_chance < 0) error->all(FLERR,"Illegal value for change_chance");
  }

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

void FixPhaseChange::init_list(int, NeighList *ptr)
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

  int nins = 0;
  int nlocal = atom->nlocal;
  int* numneigh = list->numneigh;
  double **x = atom->x;
  double **v = atom->v;
  double **vest = atom->vest;
  double **cg = atom->colorgradient;
  double *rmass = atom->rmass;
  double *rho = atom->rho;
  double *cv = atom->cv;
  double *e   = atom->e;
  dmass   = atom->de;
  int *type = atom->type;
  
  int nall;
  if (force->newton) nall = atom->nlocal + atom->nghost;
  else nall = atom->nlocal;
  for (int i = 0; i < nall; i++) {
    dmass[i] = 0.0;
  }

  for (int i = 0; i < nlocal; i++) {
    double Ti = sph_energy2t(e[i], cv[i]);
    bool isphasechange;
    if ( (Ti<Tc) || (type[i] != to_type) ) {
      isphasechange = false;
    } else if (energy_chance_flag) {
      double threshold  = (e[i] - sph_t2energy(Tc, cv[i]))/Hwv*update->dt*phase_change_rate;
      isphasechange = (random->uniform()<threshold) && isfromphasearound(i);
    } else {
      isphasechange = (random->uniform()<change_chance) && (Ti>Tt) && isfromphasearound(i);
    }
    if  (isphasechange)  {
      double coord[3];
      bool ok;
      double delta = dr;
      int natempt = 0;
      do { 
	create_newpos(x[i], cg[i], delta, coord);
	ok = insert_one_atom(coord, sublo, subhi);
	// reduce dr
	delta = 0.75*delta;
	natempt++;
      } while (!ok && natempt<10);
      
      if (!ok) {
	double delta = dr;
	natempt = 0;
	do { 
	  create_newpos_simple(x[i], delta, coord);
	  ok = insert_one_atom(coord, sublo, subhi);
	  // reduce dr
	  delta = 0.75*delta;
	  natempt++;
	} while (!ok && natempt<10);
      }

      if (ok) {
	nins++;
	/// we should take mass from neighboring from_type atoms
	double xtmp = x[i][0];
	double ytmp = x[i][1];
	double ztmp = x[i][2];
	int** firstneigh = list->firstneigh;
	int jnum = numneigh[i];
	int* jlist = firstneigh[i];
	// collect wights
	double wtotal = 0.0;
	for (int jj = 0; jj < jnum; jj++) {
	  int j = jlist[jj];
	  j &= NEIGHMASK;
	  if  ( (type[j]==from_type) && (rmass[j]>0.5*to_mass) ) {
	    double delx = xtmp - x[j][0];
	    double dely = ytmp - x[j][1];
	    double delz = ztmp - x[j][2];
	    double rsq = delx * delx + dely * dely + delz * delz;
	    double wfd;
	    if (domain->dimension == 3) {
	      wfd = sph_kernel_quintic3d(sqrt(rsq)*cutoff);
	    } else {
	      wfd = sph_kernel_quintic2d(sqrt(rsq)*cutoff);
	    }
	    wtotal+=wfd;
	  }
	}

	// take mass
	// and keep energy and momentum we are taking
	double dmom[3];
	double dmomest[3];
	dmom[0]=0.0; dmom[1]=0.0; dmom[2]=0.0;
	dmomest[0]=0.0; dmomest[1]=0.0; dmomest[2]=0.0;
	double denergy = 0.0;
	for (int jj = 0; jj < jnum; jj++) {
	  int j = jlist[jj];
	  j &= NEIGHMASK;
	  if  ( (type[j]==from_type) && (rmass[j]>0.5*to_mass) ) {
	    double delx = xtmp - x[j][0];
	    double dely = ytmp - x[j][1];
	    double delz = ztmp - x[j][2];
	    double rsq = delx * delx + dely * dely + delz * delz;
	    double wfd;
	    if (domain->dimension == 3) {
	      wfd = sph_kernel_quintic3d(sqrt(rsq)*cutoff);
	    } else {
	      wfd = sph_kernel_quintic2d(sqrt(rsq)*cutoff);
	    }
	    double dmass_aux = to_mass*wfd/wtotal;
 	    dmass[j] += dmass_aux;
	    denergy += e[j]*dmass_aux;
	    dmom[0] += v[j][0]*dmass_aux;
	    dmom[1] += v[j][1]*dmass_aux;
	    dmom[2] += v[j][2]*dmass_aux;

	    dmomest[0] += vest[j][0]*dmass_aux;
	    dmomest[1] += vest[j][1]*dmass_aux;
	    dmomest[2] += vest[j][2]*dmass_aux;
	  }
	}
	
 	// for a new atom
	int m = atom->nlocal-1;
	rmass[m] = to_mass;
	rho[m] = rho[i];
	cv[m] = cv[i];
	// TODO: think about a better momentum conservation
	v[m][0] = dmom[0]/to_mass;
	v[m][1] = dmom[1]/to_mass;
	v[m][2] = dmom[2]/to_mass;
	vest[m][0] = dmomest[0]/to_mass;
	vest[m][1] = dmomest[1]/to_mass;
	vest[m][2] = dmomest[2]/to_mass;

	// conserve energy
	double energy_aux = 0.5*(e[i] - Hwv);
	e[i] = energy_aux;
	e[m] = energy_aux;
      }
    }
  }

  /// substract mass and update energy
  comm->reverse_comm_fix(this);
  for (int i = 0; i < nlocal; i++) {
    double mold = rmass[i];
    rmass[i] -= dmass[i];
    // renormalize energy
    e[i] = e[i]*mold/rmass[i];
    dmass[i] = 0;
  }
  
  // reset global natoms
  // set tag # of new particle beyond all previous atoms
  // if global map exists, reset it now instead of waiting for comm
  // since deleting atoms messes up ghosts
  next_reneighbor += nfreq;
  int ninsall;
  MPI_Allreduce(&nins,&ninsall,1,MPI_INT,MPI_SUM,world);
  if (ninsall>0) {
    atom->natoms += ninsall;
    if (atom->tag_enable) {
      atom->tag_extend();
    }
    atom->nghost = 0;
    if (atom->map_style) {
      atom->map_init();
      atom->map_set();
    }
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
      else if (strcmp(arg[iarg+1],"lattice") == 0) 
	error->all(FLERR,"Illegal fix phase_change command: 'units lattice' "
		   "is not implemented");
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
  int flag;
  double lamda[3];
  double *newcoord;

  int nfix = modify->nfix;
  Fix **fix = modify->fix;

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
    for (int j = 0; j < nfix; j++)
      if (fix[j]->create_attribute) fix[j]->set_arrays(m);
  }
  if (flag) {
    return true;
  } else {
    return false;
  }
}

void FixPhaseChange::create_newpos_simple(double* xone, double delta, double* coord) {
  coord[0] = xone[0] + (random->uniform() - 0.5)*delta;
  coord[1] = xone[1] + (random->uniform() - 0.5)*delta;
  coord[2] = xone[2] + (random->uniform() - 0.5)*delta;
}

void FixPhaseChange::create_newpos(double* xone, double* cgone, double delta, double* coord) {
  double eij[3];
  if (domain->dimension==3) {
    double b1[3];
    b1[0] =  -cgone[1] ;
    b1[1] =  cgone[0] ;
    b1[2] =  0 ;
    double b1abs = sqrt(b1[0]*b1[0] + b1[1]*b1[1] + b1[2]*b1[2]);
    if (b1abs>CG_SMALL) {
      b1[0] = b1[0]/b1abs;     b1[1] = b1[1]/b1abs;     b1[2] = b1[2]/b1abs;
    }

    double b2[3];
    b2[0] =  -cgone[0]*cgone[1]*cgone[2]/(pow(cgone[1],2)+pow(cgone[0],2)) ;
    b2[1] =  -cgone[2]*pow(cgone[1],2)/(pow(cgone[1],2)+pow(cgone[0],2)) ;
    b2[2] =  cgone[1];
    double b2abs = sqrt(b2[0]*b2[0] + b2[1]*b2[1] + b2[2]*b2[2]);
    if (b1abs>CG_SMALL) {
      b2[0] = b2[0]/b2abs;     b2[1] = b2[1]/b2abs;     b2[2] = b2[2]/b2abs;
    }

    double atmp = random->uniform() - 0.5;
    double btmp = random->uniform() - 0.5;
    eij[0] = atmp*b1[0] + btmp*b2[0];
    eij[1] = atmp*b1[1] + btmp*b2[1];
    eij[2] = atmp*b1[2] + btmp*b2[2];

  } else {
    // TODO: find direction cheaper
    double atmp = random->uniform();
    if (atmp>0.5) atmp=1; else atmp=-1;
    eij[0] = -atmp*cgone[1];
    eij[1] = atmp*cgone[0];
    eij[2] = 0.0;
  }
  double eijabs = sqrt(eij[0]*eij[0] + eij[1]*eij[1] + eij[2]*eij[2]);
  /// TODO: add scale
  coord[0] = xone[0] + eij[0]*delta/eijabs;
  coord[1] = xone[1] + eij[1]*delta/eijabs;
  coord[2] = xone[2] + eij[2]*delta/eijabs;
}

int FixPhaseChange::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = dmass[i];
  }
  return comm_reverse;
}

/* ---------------------------------------------------------------------- */

void FixPhaseChange::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    dmass[j] += buf[m++];
  }
}

bool FixPhaseChange::isfromphasearound(int i) {
  int* numneigh = list->numneigh;
  double **x = atom->x;
  int *type = atom->type;
  double cutoff2 = cutoff*cutoff;
  double xtmp = x[i][0];
  double ytmp = x[i][1];
  double ztmp = x[i][2];
  int** firstneigh = list->firstneigh;
  int jnum = numneigh[i];
  int* jlist = firstneigh[i];
  for (int jj = 0; jj < jnum; jj++) {
    int j = jlist[jj];
    j &= NEIGHMASK;
    if (type[j]==from_type) {
	double delx = xtmp - x[j][0];
	double dely = ytmp - x[j][1];
	double delz = ztmp - x[j][2];
	double rsq = delx * delx + dely * dely + delz * delz;
	if (rsq<=cutoff2) return true;
    }
  }
  return  false;
}
