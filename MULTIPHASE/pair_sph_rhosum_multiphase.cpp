#include "pair_sph_rhosum_multiphase.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "sph_multiphase.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace SPHMultiphase;

PairSPHRhoSumMultiphase::PairSPHRhoSumMultiphase(LAMMPS *lmp) : Pair(lmp), cut(nullptr)
{
  restartinfo = 0;
  comm_forward = 1;
}

PairSPHRhoSumMultiphase::~PairSPHRhoSumMultiphase()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
  }
}

void PairSPHRhoSumMultiphase::init_style()
{
  require_rmass(atom, error);
  neighbor->add_request(this, NeighConst::REQ_FULL);
}

void PairSPHRhoSumMultiphase::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);
  if (nstep == 0 || update->ntimestep % nstep != 0) {
    comm->forward_comm(this);
    return;
  }

  double **x = atom->x;
  double *rho = atom->rho;
  double *rmass = atom->rmass;
  int *type = atom->type;
  const int dim = domain->dimension;
  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    int itype = type[i];
    double h = cut[itype][itype];
    rho[i] = w(dim, 0.0) / (dim == 3 ? h * h * h : h * h);
  }

  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    int itype = type[i];
    int *jlist = firstneigh[i];
    int jnum = numneigh[i];
    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj] & NEIGHMASK;
      int jtype = type[j];
      double delx = x[i][0] - x[j][0];
      double dely = x[i][1] - x[j][1];
      double delz = x[i][2] - x[j][2];
      double rsq = delx * delx + dely * dely + delz * delz;
      if (rsq >= cutsq[itype][jtype]) continue;
      double ih = 1.0 / cut[itype][jtype];
      rho[i] += w(dim, sqrt(rsq) * ih) * (dim == 3 ? ih * ih * ih : ih * ih);
    }
    rho[i] *= rmass[i];
  }

  comm->forward_comm(this);
}

void PairSPHRhoSumMultiphase::allocate()
{
  allocated = 1;
  int n = atom->ntypes;
  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++) setflag[i][j] = 0;
  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
  memory->create(cut, n + 1, n + 1, "pair:cut");
}

void PairSPHRhoSumMultiphase::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR, "Illegal pair_style sph/rhosum/multiphase command");
  nstep = utils::inumeric(FLERR, arg[0], false, lmp);
}

void PairSPHRhoSumMultiphase::coeff(int narg, char **arg)
{
  if (narg != 3)
    error->all(FLERR, "Incorrect args for pair_style sph/rhosum/multiphase coefficients");
  if (!allocated) allocate();
  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);
  double cut_one = utils::numeric(FLERR, arg[2], false, lmp);
  int count = 0;
  for (int i = ilo; i <= ihi; i++)
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

double PairSPHRhoSumMultiphase::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR, "All pair sph/rhosum/multiphase coeffs are not set");
  cut[j][i] = cut[i][j];
  return cut[i][j];
}

int PairSPHRhoSumMultiphase::pack_forward_comm(int n, int *list, double *buf, int, int *)
{
  double *rho = atom->rho;
  for (int i = 0; i < n; i++) buf[i] = rho[list[i]];
  return n;
}

void PairSPHRhoSumMultiphase::unpack_forward_comm(int n, int first, double *buf)
{
  double *rho = atom->rho;
  for (int i = 0; i < n; i++) rho[first + i] = buf[i];
}
