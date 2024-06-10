#!/bin/bash

#set -x

date

#-- system related dir, from module env or manually set
#MPIDIR=/share/apps/gnu-4.8.5/mpich-3.3
MPIDIR=$MPI_ROOT

MAINCONF=SeisFD3D.conf

EXEC_DIR=/share/home/litong/spherpml/FD3DspherEw/bin

#-- generate grid
echo "+ generate grid ..."
${EXEC_DIR}/seis3d_grid $MAINCONF

if [ $? -ne 0 ]; then
    printf "\ngrid generation fail! stop!\n"
    exit 1
fi

#-- generate metric
echo "+ calculate grid metric..."
${EXEC_DIR}/seis3d_metric $MAINCONF

if [ $? -ne 0 ]; then
    printf "\ncoordinate mapping calculation fail! stop!\n"
    exit 1
fi

#-- generate media
echo "+ evaluate media ..."
${EXEC_DIR}/seis3d_media $MAINCONF
if [ $? -ne 0 ]; then
    printf "\ndiscretize media fail! stop!\n"
    exit 1
fi

#-- generate source
echo "+ discrete source ..."
${EXEC_DIR}/seis3d_source $MAINCONF
if [ $? -ne 0 ]; then
    printf "\nlocate source fail! stop!\n"
    exit 1
fi

#-- generate station
echo "+ discrete station ..."
${EXEC_DIR}/seis3d_station $MAINCONF
if [ $? -ne 0 ]; then
    printf "\nlocate station fail! stop!\n"
    exit 1
fi

date

# vim:ts=4:sw=4:nu:et:ai:
