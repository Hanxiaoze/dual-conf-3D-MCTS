#!/bin/bash
cpptraj << EOF
parm nsp5_8UTE.top
parmstrip !(:1-302&!@H=)
parmwrite out nsp5_8UTE_nowater.prmtop
go
EOF
cpptraj << EOF
parm nsp5_8UTE.top
trajin nsp5_8UTE_prod1.nc
trajin nsp5_8UTE_prod2.nc
trajin nsp5_8UTE_prod3.nc
trajin nsp5_8UTE_prod4.nc
trajin nsp5_8UTE_prod5.nc
strip !(:1-302&!@H=)
trajout nsp5_8UTE_nowater.ncdf
run
EOF
cpptraj << EOF
parm nsp5_8UTE_nowater.prmtop
trajin nsp5_8UTE_nowater.ncdf
average average.rst7 rst7
run
EOF
cpptraj << EOF
parm nsp5_8UTE_nowater.prmtop
trajin nsp5_8UTE_nowater.ncdf
autoimage
reference average.rst7
rms reference
cluster c1 \
 kmeans clusters 10 randompoint maxit 500 \
 rms :22-27,41-54,138-145,164-173&!@H= \
 sieve 10 random \
 out cnumvtime.dat \
 summary summary.dat \
 info info.dat \
 cpopvtime cpopvtime.agr normframe \
 repout rep repfmt pdb \
 singlerepout singlerep.nc singlerepfmt netcdf \
 avgout avg avgfmt pdb
run
EOF
