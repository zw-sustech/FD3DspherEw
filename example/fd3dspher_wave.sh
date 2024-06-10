#!/bin/bash

#set -x

date


#-- system related dir, from module env or manually set
MPIDIR=/share/apps/mpich/3.4.1_intel_2019.5

#-- program related dir
EXE_DIR=/share/home/litong/spherpml/FD3DspherEw/bin
EXEC_WAVE=${EXE_DIR}/seis3d_wave_mpi
echo "EXEC_WAVE=$EXEC_WAVE"

PROJDIR=`pwd`
PAR_FILE=SeisFD3D.conf
EVTNM=spher
echo "PROJDIR=${PROJDIR}"
echo "EVTNM=${EVTNM}"


#RUN_SCRIPT_FILE=${PROJDIR}/runscript.sh
RUN_SCRIPT_FILE=${PROJDIR}/runscript.lsf
echo "RUN_SCRIPT_FILE=${RUN_SCRIPT_FILE}"

#-- create dir
#mkdir -p $PROJDIR
#mkdir -p $OUTPUT_DIR
#mkdir -p $GRID_DIR
#mkdir -p $MEDIA_DIR

#-- tmp dir should be in local, here is for test
#TMP_DIR=${PROJDIR}/scratch
#mkdir -p ${TMP_DIR}

#----------------------------------------------------------------------
#-- grid and mpi configurations
#----------------------------------------------------------------------

#-- total x mpi procs
NPROCS_X=5
#-- total y mpi procs
NPROCS_Y=8
#-- total mpi procs
NPROCS=$(( NPROCS_X*NPROCS_Y ))

#-------------------------------------------------------------------------------
#-- generate run script
#-------------------------------------------------------------------------------

#-- for lsf
create_script_lsf()
{
cat << ieof > ${RUN_SCRIPT_FILE}
#!/bin/bash

#-- Job Name
#BSUB -J ${EVTNM} 
#-- Queue name
#--  mars: normal, large;  earth: short, large, long
#BSUB -q normal
#-- requires number of cores (default: 1)
#BSUB -n ${NPROCS}
#-- Other specification
##BSUB -R "hname!=c013"
#-- Merge stderr with stdout, %J is the job-id
#BSUB -o stdout.%J

printf "\nUse $NPROCS CPUs on following nodes:\n"

MPI_CMD="$MPIDIR/bin/mpiexec -np ${NPROCS} ${EXEC_WAVE} ${PAR_FILE}"

printf "\nStart simualtion ...\n";
printf "%s\n\n" "\${MPI_CMD}";

time \${MPI_CMD};
if [ \$? -ne 0 ]; then
    printf "\nSimulation fail! stop!\n"
    exit 1
fi
ieof
}

#-- for pbs
create_script_pbs()
{
cat << ieof > ${RUN_SCRIPT_FILE}
#!/bin/bash

#-- Job Name
#PBS -N ${EVTNM} 
#-- Queue name, should check the naming of the system
#PBS -q normal
#-- requires number of cores
#--  be aware that pbs requires to specify number of nodes, ncpus per node
#PBS -l select=1:ncpus=${NPROCS}
#-- Merge stderr with stdout
#PBS -j oe

printf "\nUse $NPROCS CPUs on following nodes:\n"
printf "%s " \`cat ${PBS_NODEFILE} | sort\`;

MPI_CMD="$MPIDIR/bin/mpiexec -np ${NPROCS} ${EXEC_WAVE} ${PAR_FILE} 110"

printf "\nStart simualtion ...\n";
printf "%s\n\n" "\${MPI_CMD}";

time \${MPI_CMD};
if [ \$? -ne 0 ]; then
    printf "\nSimulation fail! stop!\n"
    exit 1
fi
ieof
}

#-- for directly run
create_script_run()
{
#-- create hostlist
cat << ieof > ${PROJDIR}/hostlist
server1
ieof

#-- create runscript
cat << ieof > ${RUN_SCRIPT_FILE}
#!/bin/bash

printf "\nUse $NPROCS CPUs on following nodes:\n"
printf "%s " \`cat ${PROJDIR}/hostlist | sort\`;

MPI_CMD="$MPIDIR/bin/mpiexec -machinefile ${PROJDIR}/hostlist -np $NPROCS $EXEC_WAVE $PAR_FILE 10"

printf "\nStart simualtion ...\n";
printf "%s\n\n" "'\${MPI_CMD}'";

time \${MPI_CMD};
if [ \$? -ne 0 ]; then
    printf "\nSimulation fail! stop!\n"
    exit 1
fi
ieof

chmod 755 ${RUN_SCRIPT_FILE}
}

#-------------------------------------------------------------------------------
# submit or run
#-------------------------------------------------------------------------------

#-- common steps

#-- run with lsf
echo "sumbit to lsf ..."
create_script_lsf;
bsub < ${RUN_SCRIPT_FILE}

#-- run with pbs
#echo "sumbit to pbs ..."
#create_script_pbs;
#qsub ${RUN_SCRIPT_FILE}

#-- directly run
#echo "start run script ..."
#create_script_run;
#${RUN_SCRIPT_FILE}

date

# vim:ts=4:sw=4:nu:et:ai:
