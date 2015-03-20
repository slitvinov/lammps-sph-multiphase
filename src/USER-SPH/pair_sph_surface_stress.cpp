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
#include "pair_sph_surface_stress.h"
#include "sph_kernel_quintic.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include "neighbor.h"
#include "update.h"
#include "domain.h"
#include <iostream>
#include <assert.h>
#include <vector>

using namespace LAMMPS_NS;
typedef std::vector<double> vd;
typedef std::vector<vd> vdd;

#define EPSILON 1.0e-15

/* ---------------------------------------------------------------------- */

PairSPHSurfaceStress::PairSPHSurfaceStress(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;

  // set comm size needed by this Pair

  comm_forward = 3;
  first = 1;
}

/* ---------------------------------------------------------------------- */

PairSPHSurfaceStress::~PairSPHSurfaceStress() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(alpha);
  }
}

/* ----------------------------------------------------------------------
 init specific to this pair style
 ------------------------------------------------------------------------- */

void PairSPHSurfaceStress::init_style() {
  // need a full neighbor list
  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ---------------------------------------------------------------------- */

void PairSPHSurfaceStress::compute(int eflag, int vflag) {
  int i, j, ii, jj, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz;
  double rsq, h, ih;
  int *jlist;
  
  const int ndim = domain->dimension;
  double eij[ndim];
  // neighbor list variables
  int inum, *ilist, *numneigh, **firstneigh;

  if (eflag || vflag)
    ev_setup(eflag, vflag);
  else
    evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double *rho = atom->rho;
  int *type = atom->type;
  double *rmass = atom->rmass;
  double **surface_stress = atom->surface_stress;

  // check consistency of pair coefficients
  if (first) {
    for (i = 1; i <= atom->ntypes; i++) {
      for (j = 1; i <= atom->ntypes; i++) {
        if (cutsq[i][j] > 0.0) {
          if (!setflag[i][i] || !setflag[j][j]) {
            if (comm->me == 0) {
              printf(
		     "SPH particle types %d and %d interact, but not all of their single particle properties are set.\n",
		     i, j);
            }
          }
        }
      }
    }
    first = 0;
  }

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // we use a full neighborlist here
  if (nstep != 0) {
    if ((update->ntimestep % nstep) == 0) {
      for (ii = 0; ii < inum; ii++) {
	// 3D array of color gradients cg[atomtype][dimenshion]
	const size_t nt = atom->ntypes + 1;
	std::vector<vd> cg(nt, vd(3, 0.0));
      
        i = ilist[ii];
        xtmp = x[i][0];
        ytmp = x[i][1];
        ztmp = x[i][2];
        itype = type[i];
        jlist = firstneigh[i];
        jnum = numneigh[i];
	double sigmai = rho[i]/rmass[i];
        for (jj = 0; jj < jnum; jj++) {
          j = jlist[jj];
          j &= NEIGHMASK;

          jtype = type[j];
	  if (jtype==itype) continue;
	  
          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          rsq = delx * delx + dely * dely + delz * delz;

	  if (rsq < cutsq[itype][jtype]) {
	    double r = sqrt(rsq);
	    eij[0]= delx/r;
	    eij[1]= dely/r;
	    eij[2]= delz/r;

	    h = cut[itype][jtype];
	    ih = 1.0 / h;

	    double wfd;
	    if (domain->dimension == 3) {
	      // Quintic spline
	      wfd = sph_dw_quintic3d(r*ih);
	      wfd = wfd * ih * ih * ih * ih;
	    } else {
	      wfd = sph_dw_quintic2d(r*ih);
	      wfd = wfd * ih * ih * ih;
	    }
	    double sigmaj = rho[j]/rmass[j];
	    double dphi = wfd*sigmai/sigmaj/sigmaj;

	    cg[jtype][0] += dphi*eij[0];
	    cg[jtype][1] += dphi*eij[1];
	    cg[jtype][2] += dphi*eij[2];
          } // rsq < cutsq
        } // jj loop
	for (size_t d=0; d<6; d++)
	  surface_stress[i][d] = 0.0;

	for (jtype=1; jtype<=atom->ntypes; jtype++) {
	  if (domain->dimension == 3) {
	    const double al = alpha[itype][jtype];
	    if (al<EPSILON) continue;
	    const vd&    c = cg[jtype];
	    const double c_abs=pow(pow(c[0],2)+pow(c[1],2)+pow(c[2],2),1.0/2.0);
	    if (c_abs<EPSILON) continue;
	    surface_stress[i][0]+=(-3.0*pow(c[0],2)*al+al*pow(c_abs,2))/3.0/c_abs;
	    surface_stress[i][1]+=(-3.0*pow(c[1],2)*al+al*pow(c_abs,2))/3.0/c_abs;
	    surface_stress[i][2]+=(-3.0*pow(c[2],2)*al+al*pow(c_abs,2))/3.0/c_abs;
	    surface_stress[i][3]+=-c[0]*c[1]*al/c_abs;
	    surface_stress[i][4]+=-c[0]*c[2]*al/c_abs;
	    surface_stress[i][5]+=-c[1]*c[2]*al/c_abs;
	  } else {
	    const double al = alpha[itype][jtype];
	    if (al<EPSILON) continue;
	    const vd&    c = cg[jtype];	    
	    const double c_abs=pow(pow(c[0],2)+pow(c[1],2),1.0/2.0);
	    if (c_abs<EPSILON) continue;
	    surface_stress[i][0]+=(-pow(c[0],2)+pow(c[1],2)+pow(c[2],2))*al/2.0/c_abs;
	    surface_stress[i][1]+=(pow(c[0],2)-pow(c[1],2)+pow(c[2],2))*al/2.0/c_abs;
	    surface_stress[i][3]+=-c[0]*c[1]*al/c_abs;
	  }
	}
	
      } // ii loop
      // communicate surface_stress
      comm->forward_comm_pair(this);
    }
  }
}

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairSPHSurfaceStress::allocate() {
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");

  memory->create(cut, n + 1, n + 1, "pair:cut");
  memory->create(alpha, n + 1, n + 1, "pair:alpha");
}

/* ----------------------------------------------------------------------
 global settings
 ------------------------------------------------------------------------- */

void PairSPHSurfaceStress::settings(int narg, char **arg) {
  if (narg != 1)
    error->all(FLERR,
        "Illegal number of setting arguments for pair_style sph/surface_stress");
  nstep = force->inumeric(FLERR, arg[0]);
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairSPHSurfaceStress::coeff(int narg, char **arg) {
  if (narg != 4)
    error->all(FLERR,"Incorrect number of args for sph/surface_stress coefficients");
  if (!allocated)
    allocate();

  int ilo, ihi, jlo, jhi;
  force->bounds(arg[0], atom->ntypes, ilo, ihi);
  force->bounds(arg[1], atom->ntypes, jlo, jhi);

  double cut_one = force->numeric(FLERR, arg[2]);
  double alpha_one = force->numeric(FLERR, arg[3]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut[i][j] = cut_one;
      alpha[i][j] = alpha_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0)
    error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
 init for one type pair i,j and corresponding j,i
 ------------------------------------------------------------------------- */

double PairSPHSurfaceStress::init_one(int i, int j) {
  if (setflag[i][j] == 0) {
    error->all(FLERR,"All pair sph/surface_stress coeffs are not set");
  }

  cut[j][i] = cut[i][j];
  alpha[j][i] = alpha[i][j];

  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

double PairSPHSurfaceStress::single(int i, int j, int itype, int jtype, double rsq,
    double factor_coul, double factor_lj, double &fforce) {
  fforce = 0.0;

  return 0.0;
}

/* ---------------------------------------------------------------------- */

int PairSPHSurfaceStress::pack_comm(int n, int *list, double *buf, int pbc_flag,
    int *pbc) {
  int i, j, m;
  double **surface_stress = atom->surface_stress;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = surface_stress[j][0];
    buf[m++] = surface_stress[j][1];
    buf[m++] = surface_stress[j][2];
    buf[m++] = surface_stress[j][3];
    buf[m++] = surface_stress[j][4];
    buf[m++] = surface_stress[j][5];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairSPHSurfaceStress::unpack_comm(int n, int first, double *buf) {
  int i, m, last;
  double **surface_stress = atom->surface_stress;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    surface_stress[i][0] = buf[m++];
    surface_stress[i][1] = buf[m++];
    surface_stress[i][2] = buf[m++];
    surface_stress[i][3] = buf[m++];
    surface_stress[i][4] = buf[m++];
    surface_stress[i][5] = buf[m++];
}

