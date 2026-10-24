#include "pair_sph_colorgradient.h"

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

PairSPHColorGradient::PairSPHColorGradient(LAMMPS *lmp) :
    Pair(lmp), cut(nullptr), alpha(nullptr)
{
  restartinfo = 0;
  comm_forward = 3;
}

PairSPHColorGradient::~PairSPHColorGradient()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(alpha);
  }
}

void PairSPHColorGradient::init_style()
{
  require_rmass(atom, error);
  colorgradient(atom, error);
  neighbor->add_request(this, NeighConst::REQ_FULL);
}

void PairSPHColorGradient::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);
  if (nstep == 0 || update->ntimestep % nstep != 0) {
    comm->forward_comm(this);
    return;
  }

  double **x = atom->x;
  double *rho = atom->rho;
  double *rmass = atom->rmass;
  double **cg = colorgradient(atom, error);
  int *type = atom->type;
  const int dim = domain->dimension;
  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    cg[i][0] = cg[i][1] = cg[i][2] = 0.0;
  }

  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    int itype = type[i];
    double sigmai = rho[i] / rmass[i];
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
      double r = sqrt(rsq);
      double ih = 1.0 / cut[itype][jtype];
      double wfd = dw(dim, r * ih) * (dim == 3 ? ih * ih * ih * ih : ih * ih * ih);
      double sigmaj = rho[j] / rmass[j];
      double dphi = -wfd * alpha[itype][jtype] / (sigmaj * sigmaj) * sigmai;
      cg[i][0] += dphi * delx / r;
      cg[i][1] += dphi * dely / r;
      if (dim == 3) cg[i][2] += dphi * delz / r;
    }
  }

  comm->forward_comm(this);
}

void PairSPHColorGradient::allocate()
{
  allocated = 1;
  int n = atom->ntypes;
  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++) setflag[i][j] = 0;
  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
  memory->create(cut, n + 1, n + 1, "pair:cut");
  memory->create(alpha, n + 1, n + 1, "pair:alpha");
}

void PairSPHColorGradient::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR, "Illegal pair_style sph/colorgradient command");
  nstep = utils::inumeric(FLERR, arg[0], false, lmp);
}

void PairSPHColorGradient::coeff(int narg, char **arg)
{
  if (narg != 4) error->all(FLERR, "Incorrect args for pair_style sph/colorgradient coefficients");
  if (!allocated) allocate();
  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);
  double cut_one = utils::numeric(FLERR, arg[2], false, lmp);
  double alpha_one = utils::numeric(FLERR, arg[3], false, lmp);
  int count = 0;
  for (int i = ilo; i <= ihi; i++)
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      cut[i][j] = cut_one;
      alpha[i][j] = alpha_one;
      setflag[i][j] = 1;
      count++;
    }
  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

double PairSPHColorGradient::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR, "All pair sph/colorgradient coeffs are not set");
  cut[j][i] = cut[i][j];
  alpha[j][i] = alpha[i][j];
  return cut[i][j];
}

int PairSPHColorGradient::pack_forward_comm(int n, int *list, double *buf, int, int *)
{
  double **cg = colorgradient(atom, error);
  int m = 0;
  for (int i = 0; i < n; i++) {
    int j = list[i];
    buf[m++] = cg[j][0];
    buf[m++] = cg[j][1];
    buf[m++] = cg[j][2];
  }
  return m;
}

void PairSPHColorGradient::unpack_forward_comm(int n, int first, double *buf)
{
  double **cg = colorgradient(atom, error);
  int m = 0;
  for (int i = first; i < first + n; i++) {
    cg[i][0] = buf[m++];
    cg[i][1] = buf[m++];
    cg[i][2] = buf[m++];
  }
}
