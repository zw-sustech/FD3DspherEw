#!/bin/bash
# modified from PSSClabs scripts

### Number of nodes n nodes using m cpus Processor Per Node
#PBS -l nodes=12:ppn=2+1:ppn=1
##PBS -l nodes=20:ppn=3

### Example: to request 2 VPs on each of 3 nodes and 1 VPs on 2 more nodes
## #PBS -l nodes=3:ppn=2+2:ppn=1

### Queue name
##PBS -q high
#PBS -q default
### Job name
#PBS -N d30rincfy-670

### Merge stderr with stdout
#PBS -j oe
### Mail to user
#PBS -m eb
### Declare job-non-rerunable
#PBS -r n

PBS_PWD="`pwd`";
cd "${PBS_O_WORKDIR}";
THIS_HOST="`hostname`";
MPICH_ROOT="/opt/mpich/p4-intel";
MPIRUN_BIN="${MPICH_ROOT}/bin/mpirun";
MPIEXEC_BIN="/opt/mpiexec/bin/mpiexec";
CLEANIPCS_BIN="${MPICH_ROOT}/sbin/cleanipcs";
FNM_BIN="./bin/seis3d_wave";

NPROCS="`wc -l < ${PBS_NODEFILE} | tr -d '[:blank:]'`";

hr()
{
	perl -e 'print "\n" . "-"x70 . "\n\n"';
}

print_job_info()
{
	hr;
	printf "Torque Job ID: %s\n" "${PBS_JOBID}";
	printf "\nRunning on host %s @ %s\n" "${THIS_HOST}" "`date`";
	printf "\nStarting directory was %s\n" "${PBS_PWD}";
	printf "Working directory is %s\n" "${PBS_O_WORKDIR}";
	printf "The PWD is %s\n" "`pwd`";
	printf "\nThis job runs on the following processors:\n\n\t";
	printf "%s " `cat ${PBS_NODEFILE} | sort`;
	printf "\n\n";
	printf "This job has allocated %s nodes/processors.\n" "${NPROCS}";
	hr;
}

clean_ipcs()
{
	for NODE in `cat ${PBS_NODEFILE} | sort -u`; do
		ssh ${NODE} ${CLEANIPCS_BIN};
	done;
}

run_mpirun()
{
	clean_ipcs;
	MPI_CMD="${MPIRUN_BIN} -nolocal -np ${NPROCS} -machinefile ${PBS_NODEFILE} ${FNM_BIN}";
    printf "begin simulation, please go to bed ...\n";
	printf "%s\n\n" "${MPI_CMD}";
	time ${MPI_CMD};
	#sleep 10;
}

run_mpiexec()
{
	#clean_ipcs;
	#MPIEXEC_CMD="${MPIEXEC_BIN} --comm p4 -mpich-p4-no-shmem ${FNM_BIN}";
	MPIEXEC_CMD="${MPIEXEC_BIN} -n ${NPROCS} --comm p4 ${FNM_BIN}";
    printf "begin simulation, please go to bed ...\n";
	printf "%s\n\n" "${MPIEXEC_CMD}";
	time ${MPIEXEC_CMD};
	#sleep 10;
}

main()
{
	print_job_info;
	#run_mpirun;
	run_mpiexec;
	hr;
}

time main;

