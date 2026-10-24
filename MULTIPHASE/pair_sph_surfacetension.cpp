#include "pair_sph_surfacetension.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "sph_multiphase.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace SPHMultiphase;

static constexpr double EPSILON = 1.0e-12;

static void stress(int dim, const double *e, const double *g, double *s)
{
  s[0] = s[1] = s[2] = 0.0;
  if (dim == 2) {
    double a = sqrt(g[0] * g[0] + g[1] * g[1]);
    if (a <= EPSILON) return;
    double t = (g[0] * g[0] + g[1] * g[1]) / 2;
    s[0] = (e[0] * (t - g[0] * g[0]) - g[0] * e[1] * g[1]) / a;
    s[1] = (e[1] * (t - g[1] * g[1]) - e[0] * g[0] * g[1]) / a;
  } else {
    double a = sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
    if (a <= EPSILON) return;
    double t = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]) / 3;
    s[0] = (e[0] * (t - g[0] * g[0]) - g[0] * e[2] * g[2] - g[0] * e[1] * g[1]) / a;
    s[1] = (e[1] * (t - g[1] * g[1]) - g[1] * e[2] * g[2] - e[0] * g[0] * g[1]) / a;
    s[2] = (e[2] * (t - g[2] * g[2]) - e[1] * g[1] * g[2] - e[0] * g[0] * g[2]) / a;
  }
}

PairSPHSurfaceTension::PairSPHSurfaceTension(LAMMPS *lmp) : Pair(lmp), cut(nullptr)
{
  restartinfo = 0;
}

PairSPHSurfaceTension::~PairSPHSurfaceTension()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
  }
}

void PairSPHSurfaceTension::init_style()
{
  require_rmass(atom, error);
  colorgradient(atom, error);
  Pair::init_style();
}

void PairSPHSurfaceTension::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *rho = atom->rho;
  double **cg = colorgradient(atom, error);
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  const int dim = domain->dimension;
  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    int itype = type[i];
    double Vi = rmass[i] / rho[i];
    int *jlist = firstneigh[i];
    int jnum = numneigh[i];
    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj] & NEIGHMASK;
      int jtype = type[j];
      double del[3] = {x[i][0] - x[j][0], x[i][1] - x[j][1], x[i][2] - x[j][2]};
      double rsq = del[0] * del[0] + del[1] * del[1] + del[2] * del[2];
      if (rsq >= cutsq[itype][jtype]) continue;
      double r = sqrt(rsq);
      double ih = 1.0 / cut[itype][jtype];
      double wfd = dw(dim, r * ih) * (dim == 3 ? ih * ih * ih * ih : ih * ih * ih);
      double e[3] = {del[0] / r, del[1] / r, dim == 3 ? del[2] / r : 0.0};
      double si[3], sj[3];
      stress(dim, e, cg[i], si);
      stress(dim, e, cg[j], sj);
      double Vj = rmass[j] / rho[j];
      double fs[3];
      for (int d = 0; d < 3; d++) fs[d] = (si[d] * Vi * Vi + sj[d] * Vj * Vj) * wfd;
      for (int d = 0; d < dim; d++) f[i][d] += fs[d];
      if (newton_pair || j < nlocal)
        for (int d = 0; d < dim; d++) f[j][d] -= fs[d];
    }
  }
}

void PairSPHSurfaceTension::allocate()
{
  allocated = 1;
  int n = atom->ntypes;
  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++) setflag[i][j] = 0;
  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
  memory->create(cut, n + 1, n + 1, "pair:cut");
}

void PairSPHSurfaceTension::settings(int narg, char **)
{
  if (narg != 0) error->all(FLERR, "Illegal pair_style sph/surfacetension command");
}

void PairSPHSurfaceTension::coeff(int narg, char **arg)
{
  if (narg != 3)
    error->all(FLERR, "Incorrect args for pair_style sph/surfacetension coefficients");
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

double PairSPHSurfaceTension::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR, "All pair sph/surfacetension coeffs are not set");
  cut[j][i] = cut[i][j];
  return cut[i][j];
}
