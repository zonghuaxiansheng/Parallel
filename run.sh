if [ $1 == "sim_bcast" ]
then
	echo "Copy sim_bcast to node3 & node4 ..."
	scp ./sim_bcast pp11@node3:/home/pp11/SA19011136/.experiment/Parallel/
	scp ./sim_bcast pp11@node4:/home/pp11/SA19011136/.experiment/Parallel/
	echo "Run sim_bcast ..."
	mpirun -np 12 -f mpi_config ./sim_bcast
elif [ $1 == "sim_lusplit" ]
then
	echo "Copy sim_lusplit to node3 & node4 ..."
	scp ./$1 pp11@node3:/home/pp11/SA19011136/.experiment/Parallel/
	scp ./$1 pp11@node4:/home/pp11/SA19011136/.experiment/Parallel/
	scp ./Matrix_A pp11@node3:/home/pp11/SA19011136/.experiment/Parallel/
	scp ./Matrix_L pp11@node3:/home/pp11/SA19011136/.experiment/Parallel/
	scp ./Matrix_U pp11@node3:/home/pp11/SA19011136/.experiment/Parallel/
	scp ./Matrix_A pp11@node4:/home/pp11/SA19011136/.experiment/Parallel/
	scp ./Matrix_L pp11@node4:/home/pp11/SA19011136/.experiment/Parallel/
	scp ./Matrix_U pp11@node4:/home/pp11/SA19011136/.experiment/Parallel/
	echo "Run sim_lusplit ..."
	mpirun -np 12 -f mpi_config ./$1
elif [ $1 == "sim_dishsum" ]
then
	echo "Copy sim_dishsum to node3 & node4 ..."
	scp ./$1 pp11@node3:/home/pp11/SA19011136/.experiment/Parallel/
	scp ./$1 pp11@node4:/home/pp11/SA19011136/.experiment/Parallel/
	echo "Run sim_dishsum ..."
	for ((i = 2; i < 256; i *= 2));
	do
		mpirun -np $i -f mpi_config ./$1 2>> dishsum_parallel_out.txt
	done
elif [ $1 == "sim_btreesum" ]
then
	echo "Copy sim_btreesum to node3 & node4 ..."
	scp ./$1 pp11@node3:/home/pp11/SA19011136/.experiment/Parallel/
	scp ./$1 pp11@node4:/home/pp11/SA19011136/.experiment/Parallel/
	echo "Run sim_btreesum ..."
	mpirun -np 32 -f mpi_config ./$1
elif [ $1 == "sim_foxmult" ]
then
	echo "Copy sim_foxmult to node3 & node4 ..."
	scp ./$1 pp11@node3:/home/pp11/SA19011136/.experiment/Parallel/
	scp ./$1 pp11@node4:/home/pp11/SA19011136/.experiment/Parallel/
	echo "Run sim_foxmult ..."
	mpirun -np 36 -f mpi_config ./$1 2>> fox_parallel_out
elif [ $1 == "sim_pserver" ]
then
	echo "Copy sim_pserver to node3 & node4 ..."
	scp ./$1 pp11@node3:/home/pp11/SA19011136/.experiment/Parallel/
	scp ./$1 pp11@node4:/home/pp11/SA19011136/.experiment/Parallel/
	echo "Run sim_pserver ..."
	mpirun -np 32 -f mpi_config ./$1
elif [ $1 == "sim_mc_single" ]
then
	# echo "Copy sim_mc_single to node3 & node4 ..."
	# scp ./$1 pp11@node3:/home/pp11/SA19011136/.experiment/Parallel/
	# scp ./$1 pp11@node4:/home/pp11/SA19011136/.experiment/Parallel/
	echo "Run sim_mc_single ..."
	mpirun -np 1 ./$1 2>> mc_parallel_out
elif [ $1 == "sim_mc_parallel" ]
then
	echo "Copy sim_mc_parallel to node3 & node4 ..."
	scp ./$1 pp11@node3:/home/pp11/SA19011136/.experiment/Parallel/
	scp ./$1 pp11@node4:/home/pp11/SA19011136/.experiment/Parallel/
	echo "Run sim_mc_parallel ..."
	mpirun -np 160 -f mpi_config ./$1 2>> mc_parallel_out
elif [ $1 == "sim_alltoall" ]
then
	mpirun -np 4 ./sim_alltoall
elif [ $1 == "sim_pipeline" ]
then
	mpirun -np 3 ./sim_pipeline
else
	echo "Unrecognized option: "$1
fi
