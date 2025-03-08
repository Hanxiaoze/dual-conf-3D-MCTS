#!/bin/bash



#$AMBERHOME/bin/pmemd.cuda -O -i min1.in -o min1.out -p nsp5_8UTE.top -c nsp5_8UTE.crd -r nsp5_8UTE_min1.rst -ref nsp5_8UTE.crd&&

#$AMBERHOME/bin/pmemd.cuda -O -i min2.in -o min2.out -p nsp5_8UTE.top -c nsp5_8UTE_min1.rst -r nsp5_8UTE_min2.rst&&

#$AMBERHOME/bin/pmemd.cuda -O -i heat.in -o heat.out -p nsp5_8UTE.top -c nsp5_8UTE_min2.rst -r nsp5_8UTE_heat.rst -x nsp5_8UTE_heat.nc -ref nsp5_8UTE_min2.rst&&

#$AMBERHOME/bin/pmemd.cuda -O -i equil.in -o equil.out -p nsp5_8UTE.top -c nsp5_8UTE_heat.rst -r nsp5_8UTE_equil.rst -x nsp5_8UTE_equil.nc -ref nsp5_8UTE_heat.rst&&

#mpirun -np 36 $AMBERHOME/bin/sander.MPI -O -i equil.in -o equil.out -p nsp5_8UTE.top -c nsp5_8UTE_heat.rst -r nsp5_8UTE_equil.rst -x nsp5_8UTE_equil.nc -ref nsp5_8UTE_heat.rst&&

#$AMBERHOME/bin/pmemd.cuda -O -i prod.in -o prod1.out -p nsp5_8UTE.top -c nsp5_8UTE_equil.rst -r nsp5_8UTE_prod1.rst -x nsp5_8UTE_prod1.nc&&

#$AMBERHOME/bin/pmemd.cuda -O -i prod.in -o prod2.out -p nsp5_8UTE.top -c nsp5_8UTE_prod1.rst -r nsp5_8UTE_prod2.rst -x nsp5_8UTE_prod2.nc&&

#$AMBERHOME/bin/pmemd.cuda -O -i prod.in -o prod3.out -p nsp5_8UTE.top -c nsp5_8UTE_prod2.rst -r nsp5_8UTE_prod3.rst -x nsp5_8UTE_prod3.nc&&

#$AMBERHOME/bin/pmemd.cuda -O -i prod.in -o prod4.out -p nsp5_8UTE.top -c nsp5_8UTE_prod3.rst -r nsp5_8UTE_prod4.rst -x nsp5_8UTE_prod4.nc&&

#$AMBERHOME/bin/pmemd.cuda -O -i prod.in -o prod5.out -p nsp5_8UTE.top -c nsp5_8UTE_prod4.rst -r nsp5_8UTE_prod5.rst -x nsp5_8UTE_prod5.nc


$AMBERHOME/bin/pmemd.cuda -O -i prod.in -o prod6.out -p nsp5_8UTE.top -c nsp5_8UTE_prod5.rst -r nsp5_8UTE_prod6.rst -x nsp5_8UTE_prod6.nc&&

$AMBERHOME/bin/pmemd.cuda -O -i prod.in -o prod7.out -p nsp5_8UTE.top -c nsp5_8UTE_prod6.rst -r nsp5_8UTE_prod7.rst -x nsp5_8UTE_prod7.nc&&

$AMBERHOME/bin/pmemd.cuda -O -i prod.in -o prod8.out -p nsp5_8UTE.top -c nsp5_8UTE_prod7.rst -r nsp5_8UTE_prod8.rst -x nsp5_8UTE_prod8.nc&&

$AMBERHOME/bin/pmemd.cuda -O -i prod.in -o prod9.out -p nsp5_8UTE.top -c nsp5_8UTE_prod8.rst -r nsp5_8UTE_prod9.rst -x nsp5_8UTE_prod9.nc&&

$AMBERHOME/bin/pmemd.cuda -O -i prod.in -o prod10.out -p nsp5_8UTE.top -c nsp5_8UTE_prod9.rst -r nsp5_8UTE_prod10.rst -x nsp5_8UTE_prod10.nc
