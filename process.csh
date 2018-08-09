#!/bin/csh -f
echo 'Data processing: flux.F90'
echo 'Compling the program'
pgf90 -o flux flux.F90
echo 'running the program'
flux

echo 'Data processing is over.'
rm -f flux
