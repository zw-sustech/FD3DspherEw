#!/bin/sh

#$ -cwd
#$ -m beas
#$ -j y
#$ -S /bin/sh

#$ -N fd3dspher
#$ -pe mpich 25
#$ -v MPI_HOME=/opt/mpich/p4-intel,COMMD_PORT

# $NSLOTS
# $TMPDIR/machines

echo "Got $NSLOTS slots."
PATH=$TMPDIR:$PATH

echo "begin simulation, please go to bed ..."
$MPI_HOME/bin/mpirun -np $NSLOTS -machinefile $TMPDIR/machines ./bin/seis3d_wave
#nohup mpirun -np 6 -nolocal -machinefile hostlist ./bin/seis3d_wave &

