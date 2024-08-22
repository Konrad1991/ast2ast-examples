#!/bin/bash


quarto render ./r2sundials/r2sundials.qmd
mv  ./r2sundials/r2sundials.html ./docs/r2sundials.html.html
sudo rm -r ./docs/r2sundials_files/
sudo mv ./r2sundials/r2sundials_files ./docs/r2sundials_files/

quarto render ./ParticleSwarmOptimization/pso.qmd
mv  ./ParticleSwarmOptimization/pso.html ./docs/pso.html
sudo rm -r ./docs/pso_files/
sudo mv ./ParticleSwarmOptimization/pso_files ./docs/pso_files/

quarto render ./paropt/paropt.qmd
mv  ./paropt/paropt.html ./docs/paropt.html
sudo rm -r ./docs/paropt_files/
sudo mv ./paropt/paropt_files ./docs/paropt_files/

