#!/bin/bash
#SBATCH --job-name=fortran_ifx
#SBATCH --output=fortran_ifx.%j.out
#SBATCH --partition=broadwell
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --mem=2G

module load oneapi/2022.3.1/compiler

export IFX_FLAGS="-O3 -xHost"
export FORTRAN_SOURCES="subgrid_boundary_initializer.f90 \
  macrogrid_boundary_initializer.f90 \
  subgrid_solver.f90 \
  macrogrid_solver.f90 \
  main.f90"

mkdir -p build
ifx $IFX_FLAGS -module build $FORTRAN_SOURCES -o build/main

./build/main
