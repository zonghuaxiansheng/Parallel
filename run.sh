if [ $1 == "sim_bcast" ]
then
	echo "Copy sim_bcast to node3 & node4 ..."
	scp ./sim_bcast pp11@node3:/home/pp11/SA19011136/.experiment/Parallel/
	scp ./sim_bcast pp11@node4:/home/pp11/SA19011136/.experiment/Parallel/
	echo "Run sim_bcast ..."
	mpirun -np 12 -f mpi_config ./sim_bcast
elif [ $1 == "sim_alltoall" ]
then
	mpirun -np 4 ./sim_alltoall
elif [ $1 == "sim_pipeline" ]
then
	mpirun -np 3 ./sim_pipeline
else
	echo "Unrecognized option: "$1
fi
