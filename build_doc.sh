#!/bin/bash


quarto render ./r2sundials/r2sundials.qmd

mv  ./r2sundials/r2sundials.html ./docs/index.html
sudo rm -r ./docs/r2sundials_files/
sudo mv ./r2sundials/r2sundials_files ./docs/r2sundials_files/


