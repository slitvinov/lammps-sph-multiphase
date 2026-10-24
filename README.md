# Multiphase SPH for LAMMPS

Multiphase SPH (color gradient, surface tension, phase change) for
current LAMMPS. Tested with LAMMPS 30 Sep 2026.

## Build with make

```
cp MULTIPHASE/*.cpp MULTIPHASE/*.h lammps/src/
cd lammps/src
make yes-SPH
make -j mpi
```

## Build with CMake

```
cp MULTIPHASE/*.cpp MULTIPHASE/*.h lammps/src/
mkdir lammps/build && cd lammps/build
cmake ../cmake -D PKG_SPH=on
cmake --build . -j
```

## Input

```
atom_style  sph
create_box  2 box
mass        * 1.0
fix         mp all property/atom rmass d2_cg 3 ghost yes
```

## Styles

- pair: `sph/rhosum/multiphase`, `sph/colorgradient`, `sph/taitwater/multiphase`,
  `sph/surfacetension`, `sph/heatconduction/multiphase`, `sph/heatconduction/phasechange`
- fix: `phase_change`
- compute: `sph/de/atom`

## Examples

```
cd examples/two_atoms
maxima -b colorgradient.mac
lmp_mpi -in colorgradient.lmp

cd examples/bubble_growth
LMP=lmp_mpi bash run.sh
```
