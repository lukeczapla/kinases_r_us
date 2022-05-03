#!/bin/csh
#
#

set init = Abl_kinase
set equi_prefix = equilibration
set prod_prefix = production
set prod_step   = /data/aiml/czaplal/snapshots/prod2

set platform    = OpenCL

# Equilibration
#set input_param = "-t toppar.str -p ${init}.psf -c ${init}.crd -b sysinfo.dat"
#python -u openmm_run.py -i ${equi_prefix}.inp ${input_param} -orst ${equi_prefix}.rst -odcd ${equi_prefix}.dcdi --platform ${platform} > ${equi_prefix}.out

# Production
# The checkpoint file (.chk) can't be used on a different machine (and platform) environment.
# So please make sure if you are using the same OpenCL or CUDA version while doing additional
# production steps with the checkpoint file.  The .rst file has no such issues.

set cnt = 1
set cntmax = 300

while ( ${cnt} <= ${cntmax} )
    @ pcnt = ${cnt} - 1
    set istep = ${prod_step}_${cnt}
    set pstep = ${prod_step}_${pcnt}
    cp ${istep}.out ${istep}.out.bak
    cp ${istep}.dcd ${istep}.dcd.bak
    set input_param = "-t toppar.str -p ${init}.psf -c ${init}.crd -irst ${pstep}.rst"
    python -u openmm_run.py -i ${prod_prefix}.inp ${input_param} -orst ${istep}.rst -odcd ${istep}.dcd --platform ${platform} > ${istep}.out
    @ cnt += 1
end


