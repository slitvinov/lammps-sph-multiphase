#include "compute_sph_de_atom.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;


ComputeSPHDeAtom::ComputeSPHDeAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3)
    error->all(FLERR,"Number of arguments for compute sph/de/atom command != 3");
  if (atom->esph_flag != 1)
    error->all(FLERR,"Compute sph/de/atom requires atom attribute energy, e.g. in atom_style sph");

  peratom_flag = 1;
  size_peratom_cols = 0;

  nmax = 0;
  devector = nullptr;
}


ComputeSPHDeAtom::~ComputeSPHDeAtom()
{
  memory->sfree(devector);
}


void ComputeSPHDeAtom::init()
{
  if ((comm->me == 0) && (modify->get_compute_by_style("^sph/de/atom").size() > 1))
    error->warning(FLERR, "More than one compute {}", style);
}


void ComputeSPHDeAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;


  if (atom->nmax > nmax) {
    memory->sfree(devector);
    nmax = atom->nmax;
    devector = (double *) memory->smalloc(nmax*sizeof(double),"devector/atom:devector");
    vector_atom = devector;
  }

  double *desph = atom->desph;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      devector[i] = desph[i];
    } else {
      devector[i] = 0.0;
    }
  }
}


double ComputeSPHDeAtom::memory_usage()
{
  double bytes = (double)nmax * sizeof(double);
  return bytes;
}
