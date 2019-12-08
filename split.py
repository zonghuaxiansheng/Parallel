import os
import sys

f = open("./fox_parallel_out", "r")
# text = f.readline()
cnt = 0
mpi_list = []
for line in f:
    if cnt == 0:
        print(line)
        out = line.split("-MPI-Fox-Avg-Time: ")
        print(out)
        mpi_list.append((out[0], out[1]))
        cnt = 1
    else:
        out = line.split("-Normal-Fox-Time: ")
        mpi_list.append((out[0], out[1]))
        cnt = 0
of = open("./fox_parallel_out.txt", "w")
for item in mpi_list:
    index, time = item
    of.write(index + "\t" + time)
of.close()
